library(tidyverse)
set.seed(2023)

# change log
# s = c(40, 200, 400)
# beta01
# Y1 = Y0

# 1 setups ---------------------------------------------------------------------

setups0 <- expand_grid(
  n          = 300, # c(300, 600)
  s          = c(40, 200, 400),
  strata     = c("none", "several large coarse", "many small coarse", 
                 "hybrid coarse", "triplet"),
  propensity = c("eq", "uneq"),
  p          = 400,
  SNR        = 5,
  X_cov      = c("toeplitz")
) %>% 
  filter(!(strata == "none" & propensity == "uneq")) %>% 
  mutate(case = row_number(), .before = everything())
  

# 2 generate data --------------------------------------------------------------

# .x <- 2
simdata <- map(1:nrow(setups0), ~ {
  # par
  par <- setups0 %>% slice(.x) 
  attach(par)
  # B & n1
  if (strata == "none") {
    M <- 1
    B <- rep(1, n)
    n1 <- n * 0.5
  } else if (strata == "several large coarse") {
    M <- 2
    nm <- n / M
    B <- rep(1:M, each = nm)
    if (propensity == "eq") {
      n1 <- round(rep(nm, M) * 0.5)
    } else {
      ps <- seq(0.3, 0.7, length.out = M)
      n1 <- round(rep(nm, M) * ps)
    }
  } else if (strata == "many small coarse") {
    nm <- 10
    M <- n / nm
    B <- rep(1:M, each = nm)
    if (propensity == "eq") {
      n1 <- round(rep(nm, M) * 0.5)
    } else {
      ps <- seq(0.3, 0.7, length.out = M)
      n1 <- round(rep(nm, M) * ps)
    }
  } else if (strata == "hybrid coarse") {
    n_S <- n / 3
    nm_S <- 10
    M_S <- n_S / nm_S
    B_S <- rep(1:M_S, each = nm_S)
    n_L <- n * 2 / 3
    M_L <- 2
    nm_L <- n_L / M_L
    B_L <- rep(1:M_L, each = nm_L) + M_S
    B <- c(B_S, B_L)
    M <- M_S + M_L
    if (propensity == "eq") {
      n1 <- round(c(rep(nm_S, M_S), rep(nm_L, M_L)) * 0.5)
    } else {
      ps <- seq(0.3, 0.7, length.out = M)
      n1 <- round(c(rep(nm_S, M_S), rep(nm_L, M_L)) * ps)
    }
  } else if (strata == "triplet") {
    nm <- 3
    M <- n / nm
    B <- rep(1:M, each = nm)
    if (propensity == "eq") {
      n1 <- round(rep(nm, M) * 0.5) # all are 2
    } else {
      ps <- seq(0.3, 0.7, length.out = M)
      n1 <- round(rep(nm, M) * ps) # half of 1, half of 2
    }
  }
  # X
  if (X_cov == "highly corr") {
    corr_value <- 0.8
    cov_matrix <- matrix(corr_value, nrow = p, ncol = p)
    diag(cov_matrix) <- 1  # Set diagonal elements to 1
    X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = cov_matrix)
  } else if (X_cov == "toeplitz") {
    X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = toeplitz(0.6^(0:(p - 1))))
  }
  colnames(X) <- paste0("X", 1:ncol(X))
  # beta
  #beta01 <- map(1:2, ~ c(runif(s), rep(0, p - s)))
  beta01 <- map(1:2, ~ c(runif(s, -0.1, 0.1), rep(0, p - s)))
  beta01[[1]] <- beta01[[2]]
  for (trt in c(0, 1)) {
    beta <- beta01[[trt + 1]]
    X_sc <- tibble(B, as_tibble(X)) %>% 
      group_by(B) %>% 
      mutate(across(everything(), function(x) {x - mean(x)})) %>% 
      ungroup() %>% 
      select(-B) %>% 
      as.matrix()
    fx <- X %*% beta - 2 * X_sc %*% beta
    fit <- lm(fx ~ X[,1:s])
    varfx <- var(fit$fitted.values)
    vare <- varfx / SNR
    mu <- (B / M)^(2 * trt + 1)
    Yz <- as.vector(mu + fx + MASS::mvrnorm(n, mu = 0, Sigma = vare, empirical = T))
    if (trt == 0) {Y0 <- Yz} else {Y1 <- Yz}
  }
  Y1 <- Y0
  detach(par)
  print(.x)
  return(lst(B, n1, X, Y1, Y0))
})


# 3 summary & save -------------------------------------------------------------

data_sum <- imap_dfr(simdata, ~ {
  attach(.x)
  out <- tibble(
    case = .y,
    tau = mean(Y1) - mean(Y0),
    varY1 = round(var(Y1), 1), 
    varY0 = round(var(Y0), 1), 
    R2_1_20 = round(summary(lm(Y1 ~ X[,1:20]))$r.squared, 2),
    R2_0_20 = round(summary(lm(Y0 ~ X[,1:20]))$r.squared, 2)
  )
  detach(.x)
  out
})

setups <- setups0 %>% left_join(data_sum, by = "case")

save(setups, simdata, file = "data/sim_dense.RData")
save(setups, file = "data/setups_dense.RData")






