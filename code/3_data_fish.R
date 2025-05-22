library(tidyverse)
library(glmnet)
library(moments) # kurtosis
library(HiClimR) # fastCor
library(haven) # read data
library(brglm) # propensity score
library(optmatch) # matching
set.seed(2022)
source("utility.R")
source("stra.R")

# 1 wrangle ---------------------------------------------------------------
# dt_pre

# X -----------------------------------------------------------------------

DEMO_H <- read_xpt("data/fish/DEMO_H.XPT") %>% 
  mutate(
    id = SEQN,
    gender_Male = ifelse(RIAGENDR == 1, 1, 0),
    age = RIDAGEYR,
    ratio_family_income_to_poverty = INDFMPIR,
    race_MexicanAmerican = ifelse(RIDRETH3 == 1, 1, 0),
    race_OtherHispanic = ifelse(RIDRETH3 == 2, 1, 0),
    race_NonHispanicWhite = ifelse(RIDRETH3 == 3, 1, 0),
    race_NonHispanicBlack = ifelse(RIDRETH3 == 4, 1, 0),
    race_NonHispanicAsian = ifelse(RIDRETH3 == 6, 1, 0),
    education = case_when(
      DMDEDUC2 == 1 | DMDEDUC3 < 9 | DMDEDUC3 == 55 | DMDEDUC3 == 66 ~ 1,
      DMDEDUC2 == 2 | DMDEDUC3 <= 12 ~ 2,
      DMDEDUC2 == 3 | DMDEDUC3 <= 14  ~ 3,
      DMDEDUC2 == 4 | DMDEDUC3 == 15  ~ 4,
      DMDEDUC2 == 5 ~ 5,
      TRUE ~as.numeric(NA)
    ),
    born_US = case_when(
      DMDBORN4 == 1 ~ 1,
      DMDBORN4 == 2 ~ 0,
      TRUE ~ as.numeric(NA)
    ),
    citizen_US = case_when(
      DMDCITZN == 1 ~ 1,
      DMDCITZN == 2 ~ 0,
      TRUE ~ as.numeric(NA)
    ),
    marital_Single = case_when(
      DMDMARTL %in% c(2, 3, 4, 5) ~ 1,
      DMDMARTL %in% c(1, 6) ~ 0,
      TRUE ~ as.numeric(NA)
    ),
    language_English = ifelse(SIALANG == 1, 1, 0),
    people_in_household = DMDHHSIZ,
    children_5_years_or_younger_in_HH = DMDHHSZA,
    children_6_17_years_in_HH = DMDHHSZB,
    adults_60_years_or_older_in_HH = DMDHHSZE,
    annual_household_income = case_when(
      !(INDHHIN2 %in% c(77, 99)) ~ INDHHIN2,
      INDHHIN2 == 12 ~ 5,
      INDHHIN2 == 13 ~ 4,
      TRUE ~ as.numeric(NA)
    ),
    .keep = "none"
  )

ALQ_H <- read_xpt("data/fish/ALQ_H.XPT") %>% 
  mutate(
    id = SEQN,
    avg_alcohol_per_day_past_12_mos = case_when(
      !(ALQ130 %in% c(777, 999)) ~ ALQ130,
      TRUE ~ as.numeric(NA)
    ),
    .keep = "none"
  )


DUQ_H <- read_xpt("data/fish/DUQ_H.XPT") %>% 
  mutate(
    id = SEQN,
    marijuana_hashish_ever = case_when(
      DUQ200 == 1 ~ 1,
      DUQ200 == 2 ~ 0,
      TRUE ~ as.numeric(NA)
    ),
    how_often_marijuana = case_when(
      marijuana_hashish_ever == 0 ~ 0,
      !(DUQ217 %in% c(7, 9)) ~ DUQ217,
      TRUE ~ as.numeric(NA)
    ),
    cocaine_heroin_methamphetamine_ever = case_when(
      DUQ240 == 1 ~ 1,
      DUQ240 == 2 ~ 0,
      TRUE ~ as.numeric(NA)
    ),
    inject_illegal_drug_ever = case_when(
      DUQ370 == 1 ~ 1,
      DUQ370 == 2 ~ 0,
      TRUE ~ as.numeric(NA)
    ),
    .keep = "none"
  )

SMQ_H <- read_xpt("data/fish/SMQ_H.XPT") %>% 
  mutate(
    id = SEQN,
    smoke_at_least_100_ever = case_when(
      SMQ020 == 1 ~ 1,
      SMQ020 == 2 ~ 0,
      TRUE ~ as.numeric(NA)
    ),
    .keep = "none"
  )


SMQFAM_H <- read_xpt("data/fish/SMQFAM_H.XPT") %>% 
  mutate(
    id = SEQN,
    people_household_smoke = case_when(
      !(SMD460 %in% c(777, 999)) ~ SMD460,
      TRUE ~ as.numeric(NA)
    ),
    .keep = "none"
  )

DLQ_H <- read_xpt("data/fish/DLQ_H.XPT") %>% 
  mutate(
    id = SEQN,
    difficulty_hearing = case_when(
      DLQ010 == 1 ~ 1,
      DLQ010 == 2 ~ 0,
      TRUE ~ as.numeric(NA)
    ),
    difficulty_seeing = case_when(
      DLQ020 == 1 ~ 1,
      DLQ020 == 2 ~ 0,
      TRUE ~ as.numeric(NA)
    ),
    difficulty_concentrating = case_when(
      DLQ040 == 1 ~ 1,
      DLQ040 == 2 ~ 0,
      TRUE ~ as.numeric(NA)
    ),
    difficulty_walking = case_when(
      DLQ050 == 1 ~ 1,
      DLQ050 == 2 ~ 0,
      TRUE ~ as.numeric(NA)
    ),
    difficulty_dressing_bathing = case_when(
      DLQ060 == 1 ~ 1,
      DLQ060 == 2 ~ 0,
      TRUE ~ as.numeric(NA)
    ),
    difficulty_doing_errands_alone = case_when(
      DLQ080 == 1 ~ 1,
      DLQ080 == 2 ~ 0,
      TRUE ~ as.numeric(NA)
    ),
    .keep = "none"
  )

HOQ_H <- read_xpt("data/fish/HOQ_H.XPT") %>% 
  mutate(
    id = SEQN,
    rooms_in_home = case_when(
      !(HOD050 %in% c(777, 999)) ~ HOD050,
      TRUE ~ as.numeric(NA)
    ),
    home_owned = case_when(
      HOQ065 == 1 ~ 1,
      HOQ065 == 2 | HOQ065 == 3 ~ 0,
      TRUE ~ as.numeric(NA)
    ),
    .keep = "none"
  )

WHQ_H <- read_xpt("data/fish/WHQ_H.XPT") %>% 
  mutate(
    id = SEQN,
    weight_1_yr_ago = case_when(
      !(WHD050 %in% c(7777, 9999)) ~ WHD050,
      TRUE ~ as.numeric(NA)
    ),
    height = case_when(
      !(WHD010 %in% c(7777, 9999)) ~ WHD010,
      TRUE ~ as.numeric(NA)
    ),
    .keep = "none"
  )


# Z -----------------------------------------------------------------------

DR1TOT_H <- read_xpt("data/fish/DR1TOT_H.XPT") %>% 
  mutate(
    id = SEQN,
    clam = DRD350AQ,
    crab = DRD350BQ,
    crayfish = DRD350CQ,
    lobster = DRD350DQ,
    mussel = DRD350EQ,
    oyster = DRD350FQ,
    scallop = DRD350GQ,
    shrimp = DRD350HQ,
    other_shellfish = DRD350IQ,
    unknown_shellfish = DRD350JQ,
    breaded_fish = DRD370AQ,
    tuna = DRD370BQ,
    bass = DRD370CQ,
    catfish = DRD370DQ,
    cod = DRD370EQ,
    flatfish = DRD370FQ,
    haddock = DRD370GQ,
    mackerel = DRD370HQ,
    perch = DRD370IQ,
    pike = DRD370JQ,
    pollock = DRD370KQ,
    progy = DRD370LQ,
    salmon = DRD370MQ,
    sardine = DRD370NQ,
    sea_bass = DRD370OQ,
    shark = DRD370PQ,
    swordfish = DRD370QQ,
    trout = DRD370RQ,
    walleye = DRD370SQ,
    other_fish = DRD370TQ,
    unknown_fish = DRD370UQ,
    .keep = "none"
  ) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
  mutate(
    fish_total = rowSums(across(!id)),
    fish_level_high = ifelse(fish_total > 12, 1, 0)
  ) %>% 
  select(id, fish_level_high)
  

# Y -----------------------------------------------------------------------


PBCD_H <- read_xpt("data/fish/PBCD_H.XPT") %>% 
  mutate(id = SEQN, blood_cadmium = LBXBCD, .keep = "none")



# dt_pre0 ------------------------------------------------------------------

dt_pre0 <- DEMO_H %>% 
  inner_join(HOQ_H, by = "id") %>% 
  inner_join(WHQ_H, by = "id") %>%
  inner_join(DLQ_H, by = "id") %>% 
  inner_join(ALQ_H, by = "id") %>% 
  inner_join(DUQ_H, by = "id") %>% 
  inner_join(SMQ_H, by = "id") %>% 
  inner_join(SMQFAM_H, by = "id") %>% 
  mutate(across(!id, ~ structure(., label = "predictor"))) %>% 
  inner_join(DR1TOT_H, by = "id") %>% 
  mutate(
    Z = structure(fish_level_high, label = "treatment"), 
    .keep = "unused"
  ) %>% 
  inner_join(PBCD_H, by = "id") %>% 
  mutate(
    Y = structure(log(blood_cadmium), label = "outcome"), 
    .keep = "unused"
  ) %>% 
  select(-id) %>% 
  drop_na()


# B (full matching on propensity score) ------------------------------------------------------------

fit_ps <- brglm(Z ~ gender_Male + age + ratio_family_income_to_poverty +
                  race_MexicanAmerican + race_OtherHispanic + 
                  race_NonHispanicWhite + race_NonHispanicBlack + 
                  race_NonHispanicAsian + education + smoke_at_least_100_ever, 
                family = "binomial", data = dt_pre0)

dt_ps <- dt_pre0 %>% mutate(prop_score = predict(fit_ps, type = "response"))
    
B_match <- match_on(Z ~ prop_score, data = dt_ps) %>%
  fullmatch(data = dt_ps, max.controls = 8) %>%
  as.numeric()

dt_pre <- dt_pre0 %>%
  mutate(B =  structure(B_match, label = "stratum id")) %>% 
  relocate(B, Z, Y)


dt_pre %>%
  group_by(B) %>%
  summarise(n1 = sum(Z == 1), n0 = sum(Z == 0))

dt_pre %>%
  count(Z)


# 2 feature engineering ---------------------------------------------------
# dt_pre -> dt_final

dt_final <- feature_engineering(dt_pre)


# 3 imputation -------------------------------------------------

# X, B, n1 --------------------------------------------------------------------

# XK <- dt_final %>% 
#   select(where(~ attr(.x, "label") == "predictor")) %>% 
#   as.matrix()
XK <- dt_final %>% 
  select(age, smoke_at_least_100_ever) %>% 
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

Y1 <- impute_Y(dt_final, X_imp, 1, "lasso", "1se")
Y0 <- impute_Y(dt_final, X_imp, 0, "lasso", "1se")

save(
  dt_final, B, n1, XK, X, Y1, Y0,
  file = str_glue("data/fish.RData")
)





