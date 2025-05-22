library(tidyverse)
library(glmnet)
library(moments) # kurtosis
library(HiClimR) # fastCor
library(haven) # read data
set.seed(2022)
source("utility.R")
source("stra.R")



# 1 wrangle ---------------------------------------------------------------

dt_raw <- read_dta("data/ok/OKgradesUpdate_Feb5_2010 anonymized.dta")

dt_pre <- dt_raw %>%
  filter(s_first_year == 0) %>% 
  mutate(
    s_motherbornchina = ifelse(s_mothercountofbirth == "china", 1, 0),
    s_fatherbornchina = ifelse(s_fathercountofbirth == "china", 1, 0)
  ) %>% 
  select(
    # stratum id
    s_group_quart,
    # treatment
    `T`,
    # outcome
    #avggradefall2008, 
    gpafall2008, 
    # stratum
    #s_male, 
    #s_quartile,
    # cont X
    s_hsgrade3,
    s_gpapreviousyear,
    s_age,
    s_expectedloans,
    s_expectedgrants,
    s_hourswork,
    s_attempted_credits_fall,
    # int X
    s_mtongue_english,
    s_test1correct,
    s_test2correct,
    s_mothergraddegree,
    s_mothercolldegree,
    s_motherhsdegree,
    s_motherbornchina,
    s_fathergraddegree,
    s_fathercolldegree,
    s_fatherhsdegree,
    s_fatherbornchina,
    s_liveathome,
    s_expecttowork,
    s_highfundsconcern,
    s_fouryrstograd,
    s_firstchoiceutsc
    ) %>%
  drop_na() %>% 
  # set label & type
  mutate(across(everything(), ~ {attributes(.x) <- NULL; .x})) %>% 
  mutate(
    B = structure(as.integer(as.factor(s_group_quart)), "label" = "stratum id"),
    Z = structure(as.integer(`T`), "label" = "treatment"),
    #Y = structure(avggradefall2008, "label" = "outcome"),
    Y = structure(gpafall2008, "label" = "outcome"),
    # across(
    #   c(s_male, s_quartile), 
    #   ~ structure(as.integer(.x), "label" = "stratum")
    # ),
    .keep = "unused",
    .before = everything()
  ) %>% 
  mutate(
    across(where(
      ~ is.null(attributes(.x))
    ), ~ structure(.x, "label" = "predictor")
    )
  ) %>% 
  rename_with(~ str_replace(.x, "^s_", ""))


# 2 feature engineering ---------------------------------------------------
# dt_pre -> dt_final

dt_final <- feature_engineering(dt_pre)


# 3 imputation -------------------------------------------------

# X, B, n1 --------------------------------------------------------------------

# XK <- dt_final %>% 
#   select(where(~ attr(.x, "label") == "predictor")) %>% 
#   as.matrix()
XK <- dt_final %>% 
  select(age, hsgrade3, test1correct, test2correct) %>% 
  as.matrix()
dim(XK)

X <- dt_final %>% 
  select(where(
    ~ attr(.x, "label") %in% c("predictor", "quadratic", "interaction")
  )) %>% 
  as.matrix()
dim(X)

B <- dt_final$B

n1 <- dt_final %>% 
  group_by(B) %>% 
  summarise(n = n(), n1 = sum(Z), prop = round(n1 / n, 2)) %>%
  pull(n1)

set.seed(2022)
(fit1 <- stra_adj(dt_final$Y, dt_final$Z, B, X, "lasso", "min", return_Sz = T))
(fit2 <- stra_adj(dt_final$Y, dt_final$Z, B, X, "lasso", "1se", return_Sz = T))
rbind(
  stra_unadj(dt_final$Y, dt_final$Z, B, p_a = 1),
  stra_adj(dt_final$Y, dt_final$Z, B, XK, "ols"),
  fit1$inf_res,
  fit2$inf_res
) %>% mutate(length = ci_u - ci_l)


# Y --------------------------------------------------------------------

X_imp <- dt_final %>% 
  group_by(B) %>% 
  mutate(across(
    where(~ attr(.x, "label") %in% c("predictor", "quadratic", "interaction")),  
    ~ {.x - mean(.x)}
  ), .keep = "none") %>%
  ungroup() %>% 
  select(-B) %>% 
  as.matrix()

Y1 <- impute_Y(dt_final, X_imp, 1, "lasso", "min")
Y0 <- impute_Y(dt_final, X_imp, 0, "lasso", "min")


save(
  dt_final, B, n1, XK, X, Y1, Y0,
  file = str_glue("data/ok.RData")
)

