library(tidyverse)
library(ltmle)
library(SuperLearner)
`%nin%` <- negate(`%in%`)
source("0_helper_functions.R")

risk_job_categories <- c("2","3","4","5","6","7","8","9","10","11","13","14","15","16","17","21","23")

d0 <- readRDS("data/_cleaned_ind_data.rds") %>% 
  mutate(Afake = ifelse(as.numeric(community_num) %% 2 == 0, 1, 0),
         risky_job_0 = ifelse(occupation_0 %in% risk_job_categories, 1, 0),
         risky_job_3 = ifelse(occupation_3 %in% risk_job_categories, 1, 0),
         risky_job_0_NA = as.numeric(is.na(occupation_0)),
         risky_job_3_NA = as.numeric(is.na(occupation_3)),
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

W <- c("agecat0",
       "agecat2",
       "agecat3", "ht12",
       "stable_0")
A <- "S_1_Delta_0"
Y <- "Y_0"
SLL <- c("SL.mean", "SL.glm", "SL.earth")

cps <- unique(d0$community_num)

output_cols <- c("community_num",
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

for(i in 1:length(cps)){
  print(paste0("############## community ", cps[i], ", ", i, " of ",length(cps)))
  output[i,] <- get_Yc_comm(comm = cps[i], W = W, A = A, Y = Y, SLL = SLL,
                       column_names = output_cols,
                       d0 = d0 %>% filter(community_num == cps[i]))
}
saveRDS(get_As_comm(d0 = d0, output = output), file = "data/_Yc_comm_APS.rds")


get_As_comm <- function(d0, output){
  as <- d0 %>% group_by(Afake, community_num) %>% summarise(Afake = as.logical(first(Afake)),
                                                 community_num = first(community_num),
                                                 community_name = first(community_name),
                                                 Areal = first(intervention))
  out <- left_join(output, as)
  return(out)
}

get_Yc_comm <- function(comm, d0, W, A, Y, SLL,
                   age_strat = "all", column_names){
  output <- data.frame(matrix(NA, ncol = length(column_names), nrow = 1))
  colnames(output) <- column_names
  
  if(age_strat == "youth"){
    d0 <- d0 %>% filter(agecat0 == T)
  } else if (age_strat == "adult"){
    d0 <- d0 %>% filter(agecat0 == F)
  } else if(age_strat == "adol_young"){
    d0 <- d0 %>% filter(agecat0 == T | agecat1 == T)
  } else if(age_strat == "adol_old"){
    d0 <- d0 %>% filter(agecat2 == T | agecat3 == T)
  }
  v0 <- d0 %>% filter(community_num %in% comm) %>% select(all_of(c(W, A, Y)))
  v1 <- v0 %>% drop_na
  
  output[1,"community_num"] <- comm
  output[1,"population"] <- nrow(v1)
  output[1,"number_sampled_and_tested_BL"] <- sum(v1$S_1_Delta_0)
  
  print( "#### Getting Y denom")
  print(paste("n =", nrow(v1)))

  mod <- suppressWarnings(suppressMessages(ltmle(data = v1, SL.library = SLL,
                                                 Anodes = A, Ynodes = Y, estimate.time = F,
                                                 abar = 1, stratify = T)))
  output[1,"min_g_denom"] <- min(mod$cum.g.unbounded)
  output[1,"min_g_q10_denom"] <- quantile(mod$cum.g.unbounded, probs = .1, names = F)
  output[1,"Y0_1_est"] <- mod[["estimates"]][["tmle"]]
  output[1,"Y0_1_est_lo"] <- summary(mod)$treatment$CI[1]
  output[1,"Y0_1_est_hi"] <- summary(mod)$treatment$CI[2]
  ################## Denominator: Prob(Y_0* = 0) = 1 - Prob(Y_0* = 1)
  output[1,"P_Y0_0"] <- 1 - mod[["estimates"]][["tmle"]]
  
  # Numerator: Prob(Y_1*=1, Y_0*=0):
  # Modified from https://github.com/LauraBalzer/Far-From-MCAR/blob/master/Epi_Supp_Functions.R
  # do.dynamic
  
  print( "#### Getting Y num & ratio")
  
  data.temp0 <- d0 %>% filter(community_num %in% comm) %>% 
    select(all_of(c(W, "S_1_Delta_0", "Y_0", "risky_job_3",
                    "mobile_3", "Delta_1", "Y_1")))
  data.temp <- data.temp0 %>% drop_na
  
  output[1,"number_tested_FU"] <- sum(data.temp$Delta_1)
  
  print(paste("n =", nrow(data.temp)))
  abar <- matrix(nrow=nrow(data.temp), ncol=2)
  
  abar[,1] <- 1
  abar[,2] <- data.temp$Y_0 == 0 & data.temp$S_1_Delta_0
  
  LN <- c('Y_0')
  if(age_strat == "youth"){
    LN <- c('Y_0')
  }
  data.temp <- data.temp %>%
    select(all_of(c(W, 
                    "S_1_Delta_0", LN,
                    "Delta_1", "Y_1")))
  
  est.temp <- suppressWarnings(suppressMessages(ltmle(data = data.temp,
                                                      Anodes = c('S_1_Delta_0', 'Delta_1'),
                                                      Lnodes = LN, Ynodes = 'Y_1',
                                                      abar = abar, SL.library = SLL, estimate.time = F,
                                                      stratify = T, variance.method = 'ic')))

  output[1,"min_g_num"] <- min(est.temp$cum.g.unbounded)
  output[1,"min_g_q10_num"] <- quantile(est.temp$cum.g.unbounded, probs = .1, names = F)
  output[1,"P_Y0_0_and_Y1_1"] <- est.temp[["estimates"]][["tmle"]]
  output[1,"P_Y0_0_and_Y1_1_lo"] <- summary(est.temp)$treatment$CI[1]
  output[1,"P_Y0_0_and_Y1_1_hi"] <- summary(est.temp)$treatment$CI[2]
  output[1,"Y_final"] <- output[1,"P_Y0_0_and_Y1_1"] / output[1,"P_Y0_0"]
  var_ratio <- get.var.delta(mu1 = output[1,"P_Y0_0_and_Y1_1"],
                             mu0 = output[1,"P_Y0_0"],
                             IC1 = est.temp$IC$tmle,
                             IC0 = mod$IC$tmle)
  output[1,"Y_final_lo"] <- var_ratio$est$CI.lo
  output[1,"Y_final_hi"] <- var_ratio$est$CI.hi
  return(output)
}





























