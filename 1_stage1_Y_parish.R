library(tidyverse)
library(ltmle)
library(SuperLearner)
`%nin%` <- negate(`%in%`)
source("0_helper_functions.R")

# See SAP here:
# https://arxiv.org/abs/2208.09508


d0 <- readRDS("data/_cleaned_ind_data.rds") %>% 
  mutate(Afake = ifelse(as.numeric(community_num) %% 2 == 0, 1, 0),
         agecat0 = age_0 < 13,               # 12, +1 yr year for TB baseline
         agecat1 = age_0 >= 13 & age_0 < 19, # +1
         agecat2 = age_0 >= 19 & age_0 < 26, # +1
         agecat3 = age_0 >= 26,              # +1
         ht12 = house_type == 1 | house_type == 2,
         Delta_0 = as.numeric(tstpos != -9 & !is.na(tstpos)),
         Delta_1 = as.numeric(tstpos_yr2 != -9 & !is.na(tstpos_yr2)),
         S_1_Delta_0 = as.numeric(Delta_0 & S),
         Y_0 = as.numeric(tstpos & S_1_Delta_0),
         Y_1 = as.numeric(tstpos_yr2 & Delta_1))

ps <- d0 %>% group_by(house_type, cp) %>%
  summarise(ps = sum(S) / n()) %>%
  group_by(house_type) %>% summarise(min = min(ps),
                                     q1 = quantile(ps, 0.25),
                                     median = median(ps),
                                     mean = mean(ps),
                                     q3 = quantile(ps, 0.75),
                                     max = max(ps))
d0s <- d0 %>% filter(S==1)

W <- c("agecat0",
       "agecat2",
       "agecat3", "ht12",
       "stable_0")
A <- "S_1_Delta_0"
Y <- "Y_0"
SLL <- c("SL.mean", "SL.glm", "SL.earth")

cps <- unique(d0$cp)

output_cols <- c("cp",
                 "min_g_denom", "min_g_q10_denom",
                 "min_g_num", "min_g_q10_num",
                 "population", "number_sampled_and_tested_BL",
                 "number_tested_FU",
                 "Y0_1_est", "Y0_1_est_hi", "Y0_1_est_lo",
                 "P_Y0_0", "P_Y0_0_and_Y1_1", "P_Y0_0_and_Y1_1_lo", "P_Y0_0_and_Y1_1_hi",
                 "Y_final_lo","Y_final_hi",
                 "Y_final")
output <- data.frame(matrix(NA, ncol = length(output_cols), nrow = length(cps)))
colnames(output) <- output_cols


do_Yc_all <- function(dat, seed){
  set.seed(seed)
  for(i in 1:length(cps)){
    print(paste0("############## community ", cps[i], ", ", i, " of ",length(cps)))
    output[i,] <- get_Yc(parish = cps[i], W = W, A = A, Y = Y, SLL = SLL,
                         column_names = output_cols,
                         d0 = dat %>% filter(cp == cps[i]))
  }
  return(output)
}

############################################### PRIMARY
# overall
output <- do_Yc_all(dat = d0, seed = 1)
saveRDS(get_As(d0 = d0, output = output), file = "data/_Yc_primary.rds")

# men
output <- do_Yc_all(dat = d0 %>% filter(sex_0 == 1), seed = 1)
saveRDS(get_As(d0 = d0, output = output), file = "data/_Yc_men.rds")

# women
output <- do_Yc_all(dat = d0 %>% filter(sex_0 == 0), seed = 1)
saveRDS(get_As(d0 = d0, output = output), file = "data/_Yc_women.rds")

# age <12
output <- do_Yc_all(dat = d0 %>% filter(age_0 < 11), seed = 1)
saveRDS(get_As(d0 = d0, output = output), file = "data/_Yc_agecat0.rds")

# age 12+
output <- do_Yc_all(dat = d0 %>% filter(age_0 >= 11), seed = 1)
saveRDS(get_As(d0 = d0, output = output), file = "data/_Yc_agecat123.rds")

# by HIV HH status positive
output <- do_Yc_all(dat = d0 %>% filter(ht12), seed = 1)
saveRDS(get_As(d0 = d0, output = output), file = "data/_Yc_hh_plus.rds")

# by HIV HH status negative
output <- do_Yc_all(dat = d0 %>% filter(!ht12), seed = 1)
saveRDS(get_As(d0 = d0, output = output), file = "data/_Yc_hh_neg.rds")

# HH HIV+ <12
output <- do_Yc_all(dat = d0 %>% filter(ht12, agecat0), seed = 1)
saveRDS(get_As(d0 = d0, output = output), file = "data/_Yc_hh_plus_agecat0.rds")

# HH HIV+ 12+
output <- do_Yc_all(dat = d0 %>% filter(ht12, !agecat0), seed = 1)
saveRDS(get_As(d0 = d0, output = output), file = "data/_Yc_hh_plus_agecat123.rds")

# HH HIV- <12
output <- do_Yc_all(dat = d0 %>% filter(!ht12, agecat0), seed = 1)
saveRDS(get_As(d0 = d0, output = output), file = "data/_Yc_hh_neg_agecat0.rds")

# HH HIV- 12+
output <- do_Yc_all(dat = d0 %>% filter(!ht12, !agecat0), seed = 1)
saveRDS(get_As(d0 = d0, output = output), file = "data/_Yc_hh_neg_agecat123.rds")

# HIV+ person (not HH)
# d0 %>% filter(S_1_Delta_0 == 1, hiv_0 == 1) %>%
#   group_by(cp) %>% summarise(n = n(), y1 = length(Y_1), y11 = sum(Y_1))
output <- do_Yc_all(dat = d0 %>% filter(hiv_0 == 1), seed = 1)
saveRDS(get_As(d0 = d0, output = output), file = "data/_Yc_hiv_pos_individual.rds")

# HIV- person
output <- do_Yc_all(dat = d0 %>% filter(hiv_0 != 1), seed = 1)
saveRDS(get_As(d0 = d0, output = output), file = "data/_Yc_hiv_neg_individual.rds")








############################################################################## Unadj Stage 1
do_Yc_all_unadjusted <- function(dat){
  unadj <- dat %>% filter(S == 1, Delta_0 == 1, Delta_1 == 1, Y_0 == 0) %>%
    group_by(cp) %>% summarise(Y_0_0 = n(), Y_1_1 = sum(Y_1),
                               Y_final = Y_1_1 / Y_0_0,
                               Y_final_lo = prop.test(x = Y_1_1, n = Y_0_0)$conf.int[1],
                               Y_final_hi = prop.test(x = Y_1_1, n = Y_0_0)$conf.int[2])
  return(unadj)
}

# overall
unadj <- do_Yc_all_unadjusted(dat = d0)
saveRDS(get_As(d0 = d0, output = unadj), file = "data/_Yc_unadjusted.rds")

# men
unadj <- do_Yc_all_unadjusted(dat = d0 %>% filter(sex_0 == 1))
saveRDS(get_As(d0 = d0, output = unadj), file = "data/_Yc_unadjusted_men.rds")

# women
unadj <- do_Yc_all_unadjusted(dat = d0 %>% filter(sex_0 == 0))
saveRDS(get_As(d0 = d0, output = unadj), file = "data/_Yc_unadjusted_women.rds")

# age <12
unadj <- do_Yc_all_unadjusted(dat = d0 %>% filter(agecat0))
saveRDS(get_As(d0 = d0, output = unadj), file = "data/_Yc_unadjusted_agecat0.rds")

# age 12+
unadj <- do_Yc_all_unadjusted(dat = d0 %>% filter(!agecat0))
saveRDS(get_As(d0 = d0, output = unadj), file = "data/_Yc_unadjusted_agecat123.rds")

# by HIV HH status positive
unadj <- do_Yc_all_unadjusted(dat = d0 %>% filter(ht12))
saveRDS(get_As(d0 = d0, output = unadj), file = "data/_Yc_unadjusted_hh_plus.rds")

# by HIV HH status negative
unadj <- do_Yc_all_unadjusted(dat = d0 %>% filter(!ht12))
saveRDS(get_As(d0 = d0, output = unadj), file = "data/_Yc_unadjusted_hh_neg.rds")

# HH HIV+ <12
unadj <- do_Yc_all_unadjusted(dat = d0 %>% filter(ht12, agecat0))
saveRDS(get_As(d0 = d0, output = unadj), file = "data/_Yc_unadjusted_hh_plus_agecat0.rds")

# HH HIV+ 12+
unadj <- do_Yc_all_unadjusted(dat = d0 %>% filter(ht12, !agecat0))
saveRDS(get_As(d0 = d0, output = unadj), file = "data/_Yc_unadjusted_hh_plus_agecat123.rds")

# HH HIV- <12
unadj <- do_Yc_all_unadjusted(dat = d0 %>% filter(!ht12, agecat0))
saveRDS(get_As(d0 = d0, output = unadj), file = "data/_Yc_unadjusted_hh_neg_agecat0.rds")

# HH HIV- 12+
unadj <- do_Yc_all_unadjusted(dat = d0 %>% filter(!ht12, !agecat0))
saveRDS(get_As(d0 = d0, output = unadj), file = "data/_Yc_unadjusted_hh_neg_agecat123.rds")

# HIV+ person
unadj <- do_Yc_all_unadjusted(dat = d0 %>% filter(hiv_0 == 1))
saveRDS(get_As(d0 = d0, output = unadj), file = "data/_Yc_unadjusted_hiv_pos_individual.rds")

# HIV- person
unadj <- do_Yc_all_unadjusted(dat = d0 %>% filter(hiv_0 != 1))
saveRDS(get_As(d0 = d0, output = unadj), file = "data/_Yc_unadjusted_hiv_neg_individual.rds")


# unadjusted, at 9 community level
############################################################################## Unadj Stage 1
do_Yc_all_unadjusted_comm <- function(dat){
  unadj <- dat %>% filter(S == 1, Delta_0 == 1, Delta_1 == 1, Y_0 == 0) %>%
    group_by(community_num) %>% summarise(Y_0_0 = n(), Y_1_1 = sum(Y_1),
                               Y_final = Y_1_1 / Y_0_0,
                               Y_final_lo = prop.test(x = Y_1_1, n = Y_0_0)$conf.int[1],
                               Y_final_hi = prop.test(x = Y_1_1, n = Y_0_0)$conf.int[2])
  return(unadj)
}
unadj <- do_Yc_all_unadjusted_comm(dat = d0)
temp <- get_As_comm(d0 = d0, output = unadj)

saveRDS(temp, file = "data/_Yc_unadjusted_comm.rds")





















