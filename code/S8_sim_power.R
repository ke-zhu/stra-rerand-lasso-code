library(tidyverse)
library(parallel)
library(tictoc)
library(glmnet)
source("stra.R")
set.seed(2022)

# job_id <- as.numeric(commandArgs(T)[1])
# print(str_glue("job_id is {job_id}"))
# n_rep <- 1000
# n_cores <- 48


if (is.na(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))) {
  job_id <- 2
  n_cores <- 4
  n_rep <- 20
} else {
  job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
  n_rep <- 1000
}

print(str_glue("job_id = {job_id}"))
print(str_glue("n_cores = {n_cores}"))
print(str_glue("n_rep = {n_rep}"))




load("data/sim.RData")
attach(simdata[[job_id]]) # B, n1, X, Y1, Y0
XK <- X[,1:5]
XS <- X[,1:10]

# 1 design ------------------------------------------------------------------

# tic()
# Z_design_rep <- list(
#   stra_rerand(B, n1, p_a = 1, n_rep = n_rep, parallel = T, n_cores = n_cores),
#   stra_rerand(B, n1, XK, p_a = 0.001, n_rep = n_rep, parallel = T, n_cores = n_cores)
# )
# toc()
# save(Z_design_rep, file = str_glue("output/result/power_des_c{job_id}.RData"))

load(str_glue("output/result/power_des_c{job_id}.RData"))

tau_seq <- seq(0, 1.6, length.out = 50)

#tau_seq <- c(0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6)

# tic()
# sum_tab1 <- map_dfr(tau_seq, function(tau_alt) {
#   Y1_alt <- Y0 + tau_alt
#   # 2 analysis ----------------------------------------------------------------
#   inf_tab <- map2_dfr(Z_design_rep, c("no", "yes"), function(Z_rep, rerand) {
#     p_a <- ifelse(rerand == "no", 1, 0.001)
#     mclapply(Z_rep, function(Z) {
#       Y = Y1_alt * Z + Y0 * (1 - Z)
#       rbind(
#         stra_unadj(Y, Z, B, XK, p_a = p_a),
#         stra_adj(Y, Z, B, X, "lasso", "1se")
#       )
#     }, mc.cores = n_cores) %>% 
#       map_dfr(~.) %>% 
#       mutate(case = job_id, rerand = rerand, .before = everything())
#   }) %>% 
#     mutate(across(c(rerand, estimator), as_factor))
#   # 3 summary ----------------------------------------------------------------
#   inf_tab %>% 
#     mutate(tau = tau_alt) %>% 
#     group_by(case, tau, rerand, estimator) %>% 
#     summarise(
#       bias = mean(tauhat - tau),
#       sd = sd(tauhat),
#       rmse = sqrt(mean((tauhat - tau)^2)),
#       cp = mean(ci_l <= tau & tau <= ci_u),
#       length = mean(ci_u - ci_l),
#       shat = mean(shat),
#       power = 1 - mean(ci_l <= 0 & 0 <= ci_u),
#       .groups = "drop"
#     ) %>% 
#     # 2 compare
#     group_by(case) %>% 
#     mutate(
#       `b/s%` = round(abs(bias / sd) * 100),
#       `sd%` = round(sd / sd[1] * 100),
#       `length%` = round(length / length[1] * 100),
#     ) %>% 
#     ungroup() %>% 
#     # 3 relocate
#     select(
#       case, tau,
#       rerand, estimator, 
#       bias, `b/s%`, sd, `sd%`, rmse, 
#       cp, length, `length%`, shat,
#       power
#     )
# })
# toc()

tic()
sum_tab2 <- map_dfr(tau_seq, function(tau_alt) {
  Y1_alt <- Y1 - mean(Y1) + mean(Y0) + tau_alt
  # 2 analysis ----------------------------------------------------------------
  inf_tab <- map2_dfr(Z_design_rep, c("no", "yes"), function(Z_rep, rerand) {
    p_a <- ifelse(rerand == "no", 1, 0.001)
    mclapply(Z_rep, function(Z) {
      Y = Y1_alt * Z + Y0 * (1 - Z)
      rbind(
        stra_unadj(Y, Z, B, XK, p_a = p_a),
        stra_adj(Y, Z, B, X, "lasso", "1se")
      )
    }, mc.cores = n_cores) %>% 
      map_dfr(~.) %>% 
      mutate(case = job_id, rerand = rerand, .before = everything())
  }) %>% 
    mutate(across(c(rerand, estimator), as_factor))
  # 3 summary ----------------------------------------------------------------
  inf_tab %>% 
    mutate(tau = tau_alt) %>% 
    group_by(case, tau, rerand, estimator) %>% 
    summarise(
      bias = mean(tauhat - tau),
      sd = sd(tauhat),
      rmse = sqrt(mean((tauhat - tau)^2)),
      cp = mean(ci_l <= tau & tau <= ci_u),
      length = mean(ci_u - ci_l),
      shat = mean(shat),
      power = 1 - mean(ci_l <= 0 & 0 <= ci_u),
      .groups = "drop"
    ) %>% 
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
      case, tau,
      rerand, estimator, 
      bias, `b/s%`, sd, `sd%`, rmse, 
      cp, length, `length%`, shat,
      power
    )
})
toc()

# load("/Users/ke/Library/Mobile Documents/com~apple~CloudDocs/phd/3srrlasso/code/output/power_sum_c1.RData")
# theme_set(theme_bw())
# sum_tab2 %>%
#   arrange(estimator) %>%
#   mutate(
#     Method = as_factor(case_when(
#       rerand == "no" & estimator == "unadj"  ~ "unadj",
#       rerand == "yes" & estimator == "unadj"  ~ "rerand+unadj",
#       rerand == "no" & estimator == "lasso.1se"  ~ "lasso",
#       rerand == "yes" & estimator == "lasso.1se"  ~ "rerand+lasso"
#     ))
#   ) %>%
#   ggplot(aes(x = tau, y = power, color = Method)) +
#   geom_line() +
#   geom_point() +
#   xlim(0, 0.6) +
#   ylim(0,0.2)


save(sum_tab2, file = str_glue("output/power_sum_c{job_id}.RData"))
