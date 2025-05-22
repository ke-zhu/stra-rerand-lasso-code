feature_engineering <- function(dt_pre) {
  trans <- function(.x) {
    xt <- sign(.x) * log(abs(.x) + 1)
    structure(pmin(pmax(xt, mean(xt) - 3 * sd(xt)), mean(xt) + 3 * sd(xt)),
              label = attr(.x, "label")
    )
  }
  
  # 1 main effect ---------------------------------------------------------
  # dt_pre -> dt_main

  dt_main <- dt_pre %>%
    # remove nearly constant (01 predictor) -----------
    select(where(
      ~ !(attr(.x, "label") == "predictor" & all(.x %in% c(0, 1)) &
        abs(mean(.x) - 0.5) > 0.45)
    )) %>%
    # remove constant (continuous predictor) -----------
    select(where(
      ~ !(attr(.x, "label") == "predictor" & !all(.x %in% c(0, 1)) &
        var(.x) == 0)
    )) %>%
    # transform the (continuous predictor) with (kurt > 20) ----------------
    mutate(across(where(
      ~ attr(.x, "label") == "predictor" & !all(.x %in% c(0, 1)) &
        kurtosis(.x) > 20
    ), trans)) %>%
    # remove (continuous predictor) with (kurt > 20) -----------
    select(where(
      ~ !(attr(.x, "label") == "predictor" & !all(.x %in% c(0, 1)) &
        kurtosis(.x) > 20)
    )) %>%
    # normalize the (continuous predictor) ----------------
    mutate(across(where(
      ~ attr(.x, "label") == "predictor" & !all(.x %in% c(0, 1))
    ), ~ (.x - mean(.x)) / sd(.x)))

  # 2 quadratic terms & interaction ---------------------------------------
  # dt_main -> dt_full

  inter_name <- dt_main %>%
    select(where(~ attr(.x, "label") == "predictor")) %>%
    names() %>%
    combn(2)

  dt_inter <- map2_dfc(inter_name[1, ], inter_name[2, ], ~ {
    dt_main %>%
      mutate(
        "{.x}:{.y}" := structure(.data[[.x]] * .data[[.y]], label = "interaction")
      ) %>%
      select(last_col())
  })

  dt_full <- dt_main %>%
    # add quadratic terms ---------------------------------------
    mutate(across(where(
      ~ attr(.x, "label") == "predictor" & !all(.x %in% c(0, 1))
    ), ~ structure(.x^2, label = "quadratic"),
    .names = "{.col}^2"
    )) %>%
    # add interaction ---------------------------------------
    bind_cols(dt_inter) %>%
    # remove nearly constant (01 quad & inter) -----------
    select(where(
      ~ !(attr(.x, "label") %in% c("quadratic", "interaction") &
        all(.x %in% c(0, 1)) & abs(mean(.x) - 0.5) > 0.45)
    )) %>%
    # remove constant (continuous quad & inter) -----------
    select(where(
      ~ !(attr(.x, "label") %in% c("quadratic", "interaction") &
        !all(.x %in% c(0, 1)) & var(.x) == 0)
    )) %>%
    # transform the (continuous quad & inter) with (kurt > 20) ----------------
    mutate(across(where(
      ~ attr(.x, "label") %in% c("quadratic", "interaction") &
        !all(.x %in% c(0, 1)) & kurtosis(.x) > 20
    ), trans)) %>%
    # remove (continuous quad & inter) with (kurt > 20) -----------
    select(where(
      ~ !(attr(.x, "label") %in% c("quadratic", "interaction") &
        !all(.x %in% c(0, 1)) & kurtosis(.x) > 20)
    )) %>%
    # normalize the (continuous) (quad & inter) ------------------
    mutate(across(where(
      ~ attr(.x, "label") %in% c("quadratic", "interaction") &
        !all(.x %in% c(0, 1))
    ), ~ (.x - mean(.x)) / sd(.x)))


  # 3 remove |r| > 0.95 terms & zero variance -----------------------
  # dt_full -> dt_final

  var_highcor <- dt_full %>%
    select(where(
      ~ attr(.x, "label") %in% c("predictor", "quadratic", "interaction")
    )) %>%
    as.matrix() %>%
    fastCor(upperTri = T) %>%
    as.data.frame() %>%
    rownames_to_column(var = "var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "r") %>%
    drop_na() %>%
    filter(abs(r) > 0.95) %>%
    pull(var1) %>%
    unique()

  dt_final <- dt_full %>%
    select(!all_of(var_highcor))


  # 4 check & summary -------------------------------------------------------

  # check all (01 var)'s prop between [0.05, 0.95] -------
  tmp0 <- dt_final %>%
    summarise(across(where(
      ~ attr(.x, "label") %in% c("predictor", "quadratic", "interaction") & 
        all(.x %in% c(0, 1))
    ), mean))
  if (ncol(tmp0) == 0) {
    int_var_check <- TRUE
  } else {
    tmp <- dt_final %>%
      summarise(across(where(
        ~ attr(.x, "label") %in% c("predictor", "quadratic", "interaction") & 
          all(.x %in% c(0, 1))
      ), mean)) %>%
      pivot_longer(everything()) %>%
      arrange(value)
    int_var_check <- min(tmp$value) >= 0.05 & max(tmp$value) <= 0.95
  }
  
  # check all (continuous var)'s kurtosis <= 20 ------
  tmp <- dt_final %>%
    summarise(across(where(
      ~ attr(.x, "label") %in% c("predictor", "quadratic", "interaction") &
        !all(.x %in% c(0, 1))
    ), kurtosis)) %>%
    pivot_longer(everything()) %>%
    arrange(desc(value))
  cont_var_check <- tmp$value[1] <= 20

  # check all |cor| <= 0.95 ------
  tmp <- dt_final %>%
    select(where(
      ~ attr(.x, "label") %in% c("predictor", "quadratic", "interaction")
    )) %>%
    as.matrix() %>%
    fastCor(upperTri = T) %>%
    as.data.frame() %>%
    rownames_to_column(var = "var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "r") %>%
    drop_na() %>%
    mutate(r = abs(round(r, 2))) %>%
    arrange(desc(r))
  cor_check <- tmp$r[1] <= 0.95

  # dimension -------------------------------------------------------
  p_cont <- dt_final %>%
    select(where(
      ~ attr(.x, "label") == "predictor" & !all(.x %in% c(0, 1))
    )) %>% ncol()

  p_int <- dt_final %>%
    select(where(
      ~ attr(.x, "label") == "predictor" & all(.x %in% c(0, 1))
    )) %>% ncol()

  p_quad <- dt_final %>%
    select(where(
      ~ attr(.x, "label") == "quadratic"
    )) %>% ncol()

  p_inter <- dt_final %>%
    select(where(
      ~ attr(.x, "label") == "interaction"
    )) %>% ncol()

  p_all <- p_cont + p_int + p_quad + p_inter

  check_tab <- lst(
    int_var_check, cont_var_check, cor_check,
    p_cont, p_int, p_quad, p_inter, p_all
  )
  print(check_tab)
  return(dt_final)
}


impute_Y <- function(dt_final, X_fit, z, type, sel) {
  Y <- dt_final %>% 
    select(where(~ attr(.x, "label") == "outcome")) %>% 
    pull()
  Z <- dt_final %>% 
    select(where(~ attr(.x, "label") == "treatment")) %>% 
    pull()
  idz <- which(Z == z)
  # impute
  if (type == "sharp_null") {
    Yz <- Y
  } else if (type == "matching") {
    if (z == 0) {Zmatch <- Z} else {Zmatch <- 1 - Z}
    res_match <- Match(
      Tr = Zmatch, X = X_fit, M = 1, replace = TRUE,
      Weight = 2 # the Mahalanobis distance metric
    )
    Yz <- Y
    Yz[res_match$index.treated] <- Y[res_match$index.control]
  } else if (type == "lasso") {
    fit <- cv.glmnet(X_fit[idz,], Y[idz])
    Yz <- as.vector(predict(fit, s = str_glue("lambda.{sel}"), newx = X_fit))
    Yz[idz] <- Y[idz]
  }
  # summary
  fit <- cv.glmnet(X_fit, Yz)
  var_fx <- predict(fit, s = str_glue("lambda.{sel}"), newx = X_fit) %>% var() %>% as.vector
  lst(
    type = str_glue("{type}, {z}"),
    var = var(Yz),
    kurt = kurtosis(Yz),
    R2_ols = lm(Yz ~ X_fit) %>% summary %>% pluck("r.squared"),
    R2_lasso = var_fx / var(Yz),
    s = sum(coef(fit, s = str_glue("lambda.{sel}"))[-1] != 0)
  ) %>% print
  return(Yz)
}





