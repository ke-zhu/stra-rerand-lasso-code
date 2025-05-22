library(tidyverse)
library(parallel)
library(tictoc)
library(glmnet)
library(glmnetUtils)
source("stra.R")
set.seed(2023)

job_id <- as.numeric(commandArgs(T)[1])
print(str_glue("job_id is {job_id}"))

n_rep <- 1000
n_cores <- 48

load("data/sim_dense.RData")
attach(simdata[[job_id]]) # B, n1, X, Y1, Y0
XK <- X[,1:5]
#XS <- X[,1:10]



# 0 functions -------------------------------------------------------------

stra_adj_output_lambda <- function(Y, Z, B, X, base, 
                                   sel = "1se", s_res = 1, alpha = 0.05,
                                   return_Sz = FALSE) {
  n <- length(Y)
  variable_name <- colnames(X)
  colnames(X) <- paste0("X", 1:ncol(X))
  # transformation
  YXtrans <- tibble(Y, as_tibble(X), B, Z) %>% 
    # compute weights
    group_by(B) %>% 
    mutate(
      nm = n(),
      type_m = ifelse(sum(Z == 0) > 1 & sum(Z == 1) > 1, "coarse", "fine")
    ) %>% 
    group_by(B, Z) %>% 
    mutate(
      nmz = n(),
      emz = nmz / nm,
      wF = nm^2 / (nmz * (nm - 1) * emz),
      wX = 1 / (1 - emz),
      wY = 1 - emz
    ) %>% 
    # fine transformation for all [m]
    group_by(B) %>% 
    mutate(
      across(starts_with("X"), ~ {
        sqrt(wF * wX) * (.x - mean(.x))
      }),
      Y = sqrt(wF * wY) * Y
    ) %>%
    ungroup()
  Yt <- YXtrans$Y
  Xt <- YXtrans %>% select(starts_with("X")) %>% as.matrix
  
  map(c(0, 1), ~ {
    fit <- glmnet(Xt[Z == .x,], Yt[Z == .x], intercept = T)
    fit$lambda
  })
}


stra_adj_1se_lambda <- function(Y, Z, B, X, base, 
                     sel = "1se", s_res = 1, alpha = 0.05,
                     return_Sz = FALSE) {
  n <- length(Y)
  variable_name <- colnames(X)
  colnames(X) <- paste0("X", 1:ncol(X))
  # transformation
  YXtrans <- tibble(Y, as_tibble(X), B, Z) %>% 
    # compute weights
    group_by(B) %>% 
    mutate(
      nm = n(),
      type_m = ifelse(sum(Z == 0) > 1 & sum(Z == 1) > 1, "coarse", "fine")
    ) %>% 
    group_by(B, Z) %>% 
    mutate(
      nmz = n(),
      emz = nmz / nm,
      wF = nm^2 / (nmz * (nm - 1) * emz),
      wX = 1 / (1 - emz),
      wY = 1 - emz
    ) %>% 
    # fine transformation for all [m]
    group_by(B) %>% 
    mutate(
      across(starts_with("X"), ~ {
        sqrt(wF * wX) * (.x - mean(.x))
      }),
      Y = sqrt(wF * wY) * Y
    ) %>%
    ungroup()
  Yt <- YXtrans$Y
  Xt <- YXtrans %>% select(starts_with("X")) %>% as.matrix
  # fit
  if (base == "ols") {
    bhat <- map(c(0, 1), ~ {
      lm(Yt[Z == .x] ~ Xt[Z == .x,]) %>% coef() %>% .[-1]
    }) %>% reduce(`+`)
  } else if (base == "lasso") {
    fit_res <- map(c(0, 1), ~ {
      fit <- cv.glmnet(Xt[Z == .x,], Yt[Z == .x], intercept = T)
      bhatz <- fit %>% 
        coef(s = max(
          fit[[str_glue("lambda.{sel}")]],
          min(fit$lambda[fit$nzero <= fit$glmnet.fit$nobs * s_res])
        )) %>% 
        as.matrix() %>%
        as_tibble(rownames = "var") %>% 
        filter(str_detect(var, "X")) %>% 
        pull(2)
      Sz <- which(bhatz != 0)
      lambda.1se <- fit$lambda.1se
      lst(bhatz, Sz, lambda.1se)
    })
    bhat <- fit_res[[1]]$bhatz + fit_res[[2]]$bhatz
  } else if (base == "ridge") {
    fit_res <- map(c(0, 1), ~ {
      fit <- cv.glmnet(Xt[Z == .x,], Yt[Z == .x], intercept = T, alpha = 0)
      bhatz <- fit %>% 
        coef(s = fit[[str_glue("lambda.{sel}")]]) %>% 
        as.matrix() %>%
        as_tibble(rownames = "var") %>% 
        filter(str_detect(var, "X")) %>% 
        pull(2)
      Sz <- which(bhatz != 0)
      lst(bhatz, Sz)
    })
    bhat <- fit_res[[1]]$bhatz + fit_res[[2]]$bhatz
  } else if (base == "enet") {
    fit_res <- map(c(0, 1), ~ {
      fit <- cva.glmnet(Xt[Z == .x,], Yt[Z == .x], intercept = T)
      fit_mse <- map_dbl(fit$modlist, function(.m) {
        .m$cvm[which(.m$lambda == .m$lambda.1se)]
      })
      fit_best <- fit$modlist[[which.min(fit_mse)]]
      bhatz <- fit_best %>% 
        coef("lambda.1se") %>% 
        as.matrix() %>%
        as_tibble(rownames = "var") %>% 
        filter(str_detect(var, "X")) %>% 
        pull(2)
      Sz <- which(bhatz != 0)
      lst(bhatz, Sz)
    })
    bhat <- fit_res[[1]]$bhatz + fit_res[[2]]$bhatz
  } else if (base == "lasso+ols") {
    fit_res <- map(c(0, 1), ~ {
      fit <- cv.glmnet(Xt[Z == .x,], Yt[Z == .x], intercept = T)
      bhatz <- fit %>% 
        coef(s = fit[[str_glue("lambda.{sel}")]]) %>% 
        as.matrix() %>%
        as_tibble(rownames = "var") %>% 
        filter(str_detect(var, "X")) %>% 
        pull(2)
      Sz <- which(bhatz != 0)
      if (any(bhatz != 0)) {
        bhatz_ols <- lm(Yt[Z == .x] ~ Xt[Z == .x, Sz]) %>% coef() %>% .[-1]
        bhatz[Sz] <- bhatz_ols
      }
      lst(bhatz, Sz)
    })
    bhat <- fit_res[[1]]$bhatz + fit_res[[2]]$bhatz
  }
  bhat <- ifelse(is.na(bhat), 0, bhat)
  shat <- sum(bhat != 0)
  X_sc <- tibble(B, as_tibble(X)) %>% 
    group_by(B) %>% 
    mutate(across(everything(), function(x) {x - mean(x)})) %>% 
    ungroup() %>% 
    select(-B) %>% 
    as.matrix()
  R <- as.vector(Y - X_sc %*% bhat) # residual
  # inference
  tauhat <- stra_compute_dif(B, Z, R)
  if (base == "ridge" | base == "enet") {
    sighat <- sqrt(stra_estimate_var(R, Z, B))
  } else {
    sighat <- sqrt(stra_estimate_var(R, Z, B) * n / (n - shat))
  }
  qhat <- qnorm(1 - alpha / 2)
  # name of estimator
  if (base == "ols") {
    estimator <- base
  } else {
    if (s_res == 1) {
      estimator <- str_glue("{base}.{sel}")
    } else {
      estimator <- str_glue("{base}.r{sel}/{as.character(1/s_res)}")
    }
  }
  # output
  if (return_Sz) {
    lst(
      inf_res = tibble(
        estimator,
        tauhat,
        ci_l = tauhat - qhat * sighat, 
        ci_u = tauhat + qhat * sighat,
        shat
      ),
      S0 = variable_name[fit_res[[1]]$Sz],
      S1 = variable_name[fit_res[[2]]$Sz]
    )
  } else {
    tibble(
      estimator,
      tauhat,
      ci_l = tauhat - qhat * sighat, 
      ci_u = tauhat + qhat * sighat,
      shat,
      lambda = fit_res[[1]]$lambda.1se
    )
  }
}


stra_adj_multi_lambda <- function(Y, Z, B, X, l_seq = NULL, alpha = 0.05) {
  n <- length(Y)
  variable_name <- colnames(X)
  colnames(X) <- paste0("X", 1:ncol(X))
  # transformation
  YXtrans <- tibble(Y, as_tibble(X), B, Z) %>% 
    # compute weights
    group_by(B) %>% 
    mutate(
      nm = n(),
      type_m = ifelse(sum(Z == 0) > 1 & sum(Z == 1) > 1, "coarse", "fine")
    ) %>% 
    group_by(B, Z) %>% 
    mutate(
      nmz = n(),
      emz = nmz / nm,
      wF = nm^2 / (nmz * (nm - 1) * emz),
      wX = 1 / (1 - emz),
      wY = 1 - emz
    ) %>% 
    # fine transformation for all [m]
    group_by(B) %>% 
    mutate(
      across(starts_with("X"), ~ {
        sqrt(wF * wX) * (.x - mean(.x))
      }),
      Y = sqrt(wF * wY) * Y
    ) %>%
    ungroup()
  Yt <- YXtrans$Y
  Xt <- YXtrans %>% select(starts_with("X")) %>% as.matrix
  X_sc <- tibble(B, as_tibble(X)) %>% 
    group_by(B) %>% 
    mutate(across(everything(), function(x) {x - mean(x)})) %>% 
    ungroup() %>% 
    select(-B) %>% 
    as.matrix()
  # fit
  fit_res <- map(c(0, 1), ~ {
    l_seq_z <- l_seq[[.x + 1]]
    fit <- glmnet(Xt[Z == .x,], Yt[Z == .x], intercept = T,
                  lambda = l_seq_z)
    as.matrix(fit$beta)
  })
  bhat_mat <- fit_res[[1]] + fit_res[[2]]
  qhat <- qnorm(1 - alpha / 2)
  inf_res <- apply(bhat_mat, 2, function(bhat) {
    shat <- sum(bhat != 0)
    R <- as.vector(Y - X_sc %*% bhat) # residual
    tauhat <- stra_compute_dif(B, Z, R)
    sighat <- sqrt(stra_estimate_var(R, Z, B) * n / (n - shat))
    ci_l <- tauhat - qhat * sighat
    ci_u <- tauhat + qhat * sighat
    c(tauhat, ci_l, ci_u, shat)
  }) %>% t()
  colnames(inf_res) <- c("tauhat", "ci_l", "ci_u", "shat")
  # output
  inf_res %>% 
    as_tibble() %>% 
    mutate(estimator = paste0("lasso", row_number()), .before = everything()) %>% 
    mutate(lambda = l_seq[[1]])
}


# 1 design ------------------------------------------------------------------

tic()
Z_design_rep <- list(
  stra_rerand(B, n1, p_a = 1, n_rep = n_rep, parallel = T, n_cores = n_cores),
  stra_rerand(B, n1, XK, p_a = 0.001, n_rep = n_rep, parallel = T, n_cores = n_cores)
)
toc()
save(Z_design_rep, file = str_glue("output/result/lambda_des_c{job_id}.RData"))


# 2 analysis ----------------------------------------------------------------

# Z_rep <- Z_design_rep[[1]]
# Z <- Z_design_rep[[1]][[1]]
# p_a <- 0.001

tic()
inf_tab <- map2_dfr(Z_design_rep, c("no", "yes"), function(Z_rep, rerand) {
  p_a <- ifelse(rerand == "no", 1, 0.001)
  # output lambda seq
  Z_one <- Z_rep[[1]]
  Y_one <- Y1 * Z_one + Y0 * (1 - Z_one)
  l_seq <- stra_adj_output_lambda(Y_one, Z_one, B, X, "lasso")
  # inference
  tm <- mclapply(Z_rep, function(Z) {
    Y = Y1 * Z + Y0 * (1 - Z)
    tm <- rbind(
      stra_adj_1se_lambda(Y, Z, B, X, "lasso"),
      stra_adj_multi_lambda(Y, Z, B, X, l_seq = l_seq)
    )
  }, mc.cores = n_cores) %>% 
    map_dfr(~.) %>% 
    mutate(case = job_id, rerand = rerand, .before = everything())
}) %>% 
  mutate(across(c(rerand, estimator), as_factor))
toc()
save(inf_tab, file = str_glue("output/result/lambda_inf_c{job_id}.RData"))


# 3 summary ----------------------------------------------------------------

tic()
sum_tab <- inf_tab %>% 
  mutate(tau = setups$tau[job_id]) %>% 
  group_by(case, rerand, estimator) %>% 
  nest(data = c(tauhat, ci_l, ci_u, shat, tau, lambda)) %>%
  # 1 compute summary statistics (mean & se) of sim
  mutate(data = map(data, function(data_i) {
    # 1.1 mean of sim
    res_mean <- data_i %>% 
      summarise(
        bias = mean(tauhat - tau),
        sd = sd(tauhat),
        rmse = sqrt(mean((tauhat - tau)^2)),
        cp = mean(ci_l <= tau & tau <= ci_u),
        length = mean(ci_u - ci_l),
        shat = mean(shat),
        lambda = mean(lambda)
      )
    # # 1.2 standard error of sim (by bootstrap)
    # res_se <- map_dfr(1:500, ~ {
    #   data_i %>% 
    #     slice_sample(prop = 1, replace = TRUE) %>% 
    #     summarise(
    #       bias = mean(tauhat - tau),
    #       sd = sd(tauhat),
    #       rmse = sqrt(mean((tauhat - tau)^2)),
    #       cp = mean(ci_l <= tau & tau <= ci_u),
    #       length = mean(ci_u - ci_l),
    #       shat = mean(shat)
    #     )
    # }) %>% 
    #   # compute standard error
    #   summarise(across(everything(), ~ sd(.x), .names = "{.col}_se"))
    # # aggregate
    # bind_cols(res_mean, res_se)
    res_mean
  })) %>% 
  unnest(data) %>% 
  # 2 compare
  group_by(case) %>% 
  mutate(
    `b/s%` = round(abs(bias / sd) * 100),
    `sd%` = round(sd / sd[1] * 100),
    `length%` = round(length / length[1] * 100),
  ) %>% 
  ungroup() %>% 
  # 3 relocate
  # select(
  #   case, rerand, estimator, 
  #   bias, bias_se, `b/s%`, sd, sd_se, `sd%`, rmse, rmse_se,
  #   cp, cp_se, length, length_se, `length%`, shat, shat_se
  # )
  select(
    case, rerand, estimator, 
    bias, `b/s%`, sd, `sd%`, rmse, 
    cp, length, `length%`, shat, lambda
  )
toc()
save(sum_tab, file = str_glue("output/lambda_sum_c{job_id}.RData"))
