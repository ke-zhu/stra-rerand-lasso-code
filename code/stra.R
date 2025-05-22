
# design: stra_rerand
# analysis: stra_unadj / stra_adj

# stra_rerand -------------------------------------------------------------

stra_compute_dif <- function(B, Z, R) {
  map(sort(unique(B)), function(m) {
    mean(B == m) * (mean(R[B == m & Z == 1]) - mean(R[B == m & Z == 0]))
  }) %>% reduce(`+`)
}

stra_compute_cov_inv <- function(B, n1, X) {
  map(sort(unique(B)), function(m) {
    pim <- mean(B == m)
    em <- n1[m] / sum(B == m)
    pim * cov(X[B == m, ,drop = F]) / (em * (1 - em)) / nrow(X)
  }) %>% reduce(`+`) %>% solve()
}

stra_compute_M <- function(B, Z, X, cov_inv = NULL) {
  if (is.null(cov_inv)) {
    n1 <- map_dbl(sort(unique(B)), function(m) {sum(Z[B == m])})
    cov_inv <- stra_compute_cov_inv(B, n1, X)
  }
  dif <- apply(X, 2, function(x) {stra_compute_dif(B, Z, x)})
  as.vector(t(dif) %*% cov_inv %*% dif)
}

stra_randomize <- function(B, n1) {
  z_pos <- map(sort(unique(B)), function(m) {
    sample(which(B == m), n1[m], replace = FALSE)
  }) %>% unlist
  Z <- rep(0, length(B))
  Z[z_pos] <- 1
  Z
}

stra_rerand_one <- function(B, n1, X, a, max_iter, cov_inv = NULL) {
  Z <- stra_randomize(B, n1)
  if (a != Inf) {
    M <- stra_compute_M(B, Z, X, cov_inv)
    if (max_iter == Inf) {
      while (M > a) {
        Z <- stra_randomize(B, n1)
        M <- stra_compute_M(B, Z, X, cov_inv)
      }
    } else {
      if (M > a) {
        Zbest <- Z
        Mbest <- M
        for (iter in 2:max_iter) {
          Z <- stra_randomize(B, n1)
          M <- stra_compute_M(B, Z, X, cov_inv)
          if (M < Mbest) {
            Zbest <- Z
            Mbest <- M
          }
          if (Mbest <= a) break
        }
        # output the best rather than the last
        Z <- Zbest
      }
    }
  }
  Z
}

stra_rerand <- function(B, n1, X = NULL, p_a = 0.001, n_rep = 1,
                        parallel = FALSE, n_cores = 4, max_iter = Inf) {
  if (p_a < 1) {
    a <- qchisq(p = p_a, df = ncol(X))
    cov_inv <- stra_compute_cov_inv(B, n1, X)
  } else {
    a <- Inf
    cov_inv <- NULL
  }
  if (parallel) {
    mclapply(1:n_rep, function(x) {stra_rerand_one(B, n1, X, a, max_iter, cov_inv)}, 
             mc.cores = n_cores)
  } else {
    lapply(1:n_rep, function(x) {stra_rerand_one(B, n1, X, a, max_iter, cov_inv)})
  }
}


# stra_unadj -------------------------------------------------------------

stra_estimate_var <- function(R, Z, B) {
  n <- length(R)
  stra_type <- map_chr(sort(unique(B)), function(m) {
    n1m <- sum(Z == 1 & B == m)
    n0m <- sum(Z == 0 & B == m)
    ifelse(n1m > 1 & n0m > 1, "coarse", "fine")
  })
  if (any(stra_type == "coarse")) {
    stra_id_c <- sort(unique(B))[stra_type == "coarse"]
    sighat2_c <- map(stra_id_c, function(m) {
      sighat2m <- map(c(0, 1), function(z) {
        var(R[B == m & Z == z]) / sum(B == m & Z == z)
      }) %>% reduce(`+`)
      pim <- mean(B == m)
      pim^2 * sighat2m
    }) %>% reduce(`+`)
  } else {
    sighat2_c <- 0
  }
  if (any(stra_type == "fine")) {
    stra_id_f <- sort(unique(B))[stra_type == "fine"]
    id_f <- which(B %in% stra_id_f)
    n_f <- length(id_f)
    tauhat_f <- stra_compute_dif(B[id_f], Z[id_f], R[id_f])
    sum1 <- map(stra_id_f, function(m) {
      nm <- sum(B == m)
      wm <- nm^2 / (n_f - 2 * nm)
      tauhat_m <- mean(R[B == m & Z == 1]) - mean(R[B == m & Z == 0])
      wm * (tauhat_m - tauhat_f)^2
    }) %>% reduce(`+`)
    sum2 <- map(stra_id_f, function(m) {
      nm <- sum(B == m)
      wm <- nm^2 / (n_f - 2 * nm)
      wm
    }) %>% reduce(`+`)
    sighat2_f <- (n_f / n)^2 / (n_f + sum2) * sum1 
  } else {
    sighat2_f <- 0
  }
  sighat2_c + sighat2_f
}

stra_unadj <- function(Y, Z, B, X = NULL, p_a = 1, alpha = 0.05, n_sim = 10^4) {
  tauhat <- stra_compute_dif(B, Z, Y)
  sighat <- sqrt(stra_estimate_var(Y, Z, B))
  if (p_a < 1) {
    # R2
    n1 <- map_dbl(sort(unique(B)), function(m) {sum(Z[B == m])})
    Sigma_XX_inv <- stra_compute_cov_inv(B, n1, X)
    Sigma_hat_tX <- map(sort(unique(B)), function(m) {
      pim <- mean(B == m)
      SmXY <- map(c(0, 1), function(z) {
        emz <- ifelse(z == 1, mean(Z[B == m]), 1 - mean(Z[B == m]))
        if (sum(B == m & Z == z) > 1) {
          SmXYz <- cov(X[B == m & Z == z, ,drop = F], Y[B == m & Z == z])
        } else {
          Xmc <- scale(X[B == m, ,drop = F], center = T, scale = F)
          nm <- sum(B == m)
          SmXYz <- nm / (nm - 1) * t(Xmc[Z[B == m] == z, ,drop = F]) %*% 
            Y[B == m & Z == z]
        }
        SmXYz / emz
      }) %>% reduce(`+`)
      pim * SmXY / nrow(X)
    }) %>% reduce(`+`)
    R2hat <- (t(Sigma_hat_tX) %*% Sigma_XX_inv %*% Sigma_hat_tX / sighat^2) %>% 
      min(1) %>% max(0)
    # simulated distribution
    k <- ncol(X)
    a <- qchisq(p_a, k)
    chisq <- rchisq(n_sim / p_a, k)
    trunc_chisq <- chisq[chisq <= a]
    n_sim_emp <- length(trunc_chisq)
    S <- sample(c(1, -1), n_sim_emp, replace = T)
    beta <- rbeta(n_sim_emp, 1 / 2, (k - 1) / 2)
    L <- trunc_chisq * S * sqrt(beta)
    e0 <- rnorm(n_sim_emp)
    dhat <- sqrt(1 - R2hat) * e0 + sqrt(R2hat) * L
    qhat <- (quantile(dhat, 1 - alpha / 2) - quantile(dhat, alpha / 2)) / 2
  } else {
    qhat <- qnorm(1 - alpha / 2)
  }
  tibble(
    estimator = "unadj",
    tauhat,
    ci_l = tauhat - qhat * sighat, 
    ci_u = tauhat + qhat * sighat,
    shat = 0
  )
}

# stra_adj ---------------------------------------------------------------

stra_adj <- function(Y, Z, B, X, base = "lasso", 
                     sel = "1se", s_res = 1, alpha = 0.05,
                     return_Sz = FALSE, force_adj = NULL) {
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
  } else if (base == "relaxed lasso") {
    fit_res <- map(c(0, 1), ~ {
      fit <- cv.glmnet(Xt[Z == .x,], Yt[Z == .x], intercept = T,
                       relax = TRUE, gamma = 0)
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
  } else if (base == "rescaled enet") {
    fit_res <- map(c(0, 1), ~ {
      fit <- cva.glmnet(Xt[Z == .x,], Yt[Z == .x], intercept = T)
      fit_mse <- map_dbl(fit$modlist, function(.m) {
        .m$cvm[which(.m$lambda == .m$lambda.1se)]
      })
      id_best <- which.min(fit_mse)
      fit_best <- fit$modlist[[id_best]]
      alpha_best <- fit$alpha[id_best]
      lambda_best <- fit_best$lambda.1se
      bhatz_naive <- fit_best %>% 
        coef("lambda.1se") %>% 
        as.matrix() %>%
        as_tibble(rownames = "var") %>% 
        filter(str_detect(var, "X")) %>% 
        pull(2)
      bhatz <- bhatz_naive * (1 + lambda_best * (1 - alpha_best))
      Sz <- which(bhatz != 0)
      lst(bhatz, Sz)
    })
    bhat <- fit_res[[1]]$bhatz + fit_res[[2]]$bhatz
  } else if (base == "lasso (force adj)") {
    fit_res <- map(c(0, 1), ~ {
      pf <- rep(1, ncol(Xt[Z == .x,]))
      pf[force_adj] <- 0
      fit <- cv.glmnet(Xt[Z == .x,], Yt[Z == .x], intercept = T,
                       penalty.factor = pf)
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
  if (base == "ridge" | base == "naive enet" | base == "enet") {
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
      shat
    )
  }
}
