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


load("data/sim_varest.RData")
attach(simdata[[job_id]]) # B, n1, X, Y1, Y0
XK <- X[,1:5]
#XS <- X[,1:10]

# 1 design ------------------------------------------------------------------

tic()
Z_design_rep <- list(
  stra_rerand(B, n1, p_a = 1, n_rep = n_rep, parallel = T, n_cores = n_cores),
  stra_rerand(B, n1, XK, p_a = 0.001, n_rep = n_rep, parallel = T, n_cores = n_cores)
)
toc()
save(Z_design_rep, file = str_glue("output/result/varest_des_c{job_id}.RData"))


# 2 analysis ----------------------------------------------------------------

#Z <- Z_design_rep[[1]][[1]]

tic()
inf_tab <- map2_dfr(Z_design_rep, c("no", "yes"), function(Z_rep, rerand) {
  p_a <- ifelse(rerand == "no", 1, 0.001)
  mclapply(Z_rep, function(Z) {
    Y = Y1 * Z + Y0 * (1 - Z)
    n <- length(Y)
    inf1 <- stra_unadj(Y, Z, B, XK, p_a = p_a)
    inf2 <- stra_adj(Y, Z, B, X, "lasso")
    inf3 <- inf2 %>% 
      mutate(
        estimator = "lasso_alt",
        qhat = qnorm(1 - 0.05 / 2),
        sighat = (ci_u - ci_l) / 2 / qhat,
        sighat_alt = sqrt(sighat^2 * (n - shat) / n),
        ci_l = tauhat - qhat * sighat_alt, 
        ci_u = tauhat + qhat * sighat_alt
      ) %>% 
      select(-qhat, -sighat, -sighat_alt)
    rbind(inf1, inf2, inf3)
  }, mc.cores = n_cores) %>% 
    map_dfr(~.) %>% 
    mutate(case = job_id, rerand = rerand, .before = everything())
}) %>% 
  mutate(across(c(rerand, estimator), as_factor))
toc()
save(inf_tab, file = str_glue("output/result/varest_inf_c{job_id}.RData"))


# 3 summary ----------------------------------------------------------------

tic()
sum_tab <- inf_tab %>% 
  mutate(tau = setups$tau[job_id]) %>% 
  group_by(case, rerand, estimator) %>% 
  nest(data = c(tauhat, ci_l, ci_u, shat, tau)) %>%
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
        shat = mean(shat)
      )
    # 1.2 standard error of sim (by bootstrap)
    res_se <- map_dfr(1:500, ~ {
      data_i %>% 
        slice_sample(prop = 1, replace = TRUE) %>% 
        summarise(
          bias = mean(tauhat - tau),
          sd = sd(tauhat),
          rmse = sqrt(mean((tauhat - tau)^2)),
          cp = mean(ci_l <= tau & tau <= ci_u),
          length = mean(ci_u - ci_l),
          shat = mean(shat)
        )
    }) %>% 
      # compute standard error
      summarise(across(everything(), ~ sd(.x), .names = "{.col}_se"))
    # aggregate
    bind_cols(res_mean, res_se)
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
  select(
    case, rerand, estimator, 
    bias, bias_se, `b/s%`, sd, sd_se, `sd%`, rmse, rmse_se,
    cp, cp_se, length, length_se, `length%`, shat, shat_se
  )
toc()
save(sum_tab, file = str_glue("output/varest_sum_c{job_id}.RData"))
