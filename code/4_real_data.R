library(tidyverse)
library(parallel)
library(tictoc)
library(glmnet)
source("stra.R")
set.seed(2022)

test <- F
if (test) {
  job_id <- 1
  n_rep <- 4
  n_cores <- 4
} else {
  job_id <- as.numeric(commandArgs(T)[1])
  n_rep <- 1000
  n_cores <- 24
}

tau_true <- NA
if (job_id == 1) {
  load("data/ok.RData") # B, n1, XK, X, Y1, Y0
  tau_true[job_id] <- mean(Y1 - Y0)
} else {
  load("data/fish.RData") # B, n1, XK, X, Y1, Y0
  tau_true[job_id] <- mean(Y1 - Y0)
}

# 1 design ------------------------------------------------------------------

tic()
Z_design_rep <- list(
  stra_rerand(B, n1, p_a = 1, n_rep = n_rep, parallel = T, n_cores = n_cores),
  stra_rerand(B, n1, XK, p_a = 0.001, n_rep = n_rep, parallel = T, n_cores = n_cores)
)
toc()
save(Z_design_rep, file = str_glue("output/result/rd_des_c{job_id}.RData"))


# 2 analysis ----------------------------------------------------------------

tic()
inf_tab <- map2_dfr(Z_design_rep, c("no", "yes"), function(Z_rep, rerand) {
  p_a <- ifelse(rerand == "no", 1, 0.001)
  mclapply(Z_rep, function(Z) {
    Y = Y1 * Z + Y0 * (1 - Z)
    rbind(
      stra_unadj(Y, Z, B, XK, p_a = p_a),
      stra_adj(Y, Z, B, XK, "ols"),
      stra_adj(Y, Z, B, X, "lasso", "min"),
      stra_adj(Y, Z, B, X, "lasso", "min", s_res = 1/3),
      stra_adj(Y, Z, B, X, "lasso", "min", s_res = 1/10),
      stra_adj(Y, Z, B, X, "lasso", "min", s_res = 1/60),
      stra_adj(Y, Z, B, X, "lasso", "min", s_res = 1/120),
      stra_adj(Y, Z, B, X, "lasso", "1se")
    )
  }, mc.cores = n_cores) %>% 
    map_dfr(~.) %>% 
    mutate(case = job_id, rerand = rerand, .before = everything())
}) %>% 
  mutate(across(c(rerand, estimator), as_factor))
toc()
save(inf_tab, file = str_glue("output/result/rd_inf_c{job_id}.RData"))


# 3 summary ----------------------------------------------------------------

tic()
sum_tab <- inf_tab %>% 
  mutate(tau = tau_true[job_id]) %>% 
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
save(sum_tab, file = str_glue("output/rd_sum_c{job_id}.RData"))


