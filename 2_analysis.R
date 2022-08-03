library(ltmle)
library(tidyverse)
library(knitr)
library(glmnet)
library(coefplot)
library(glmnet)
library(SuperLearner)
source("0_helper_functions.R")

#################################
blind <- F
#if(blind){
#  A <- "Afake"
#} else {
A <- "Areal"
#}
#################################

W <- c("hiv0", "alc0")#"tb0")#, "mob0", "hr0", "pphh", "nvill") # TB data not reliable.
Y <- "Y_final"
SLL <- c("SL.mean", "SL.glm", "SL.earth")#, "SL.glmnet.9folds")
#SLL <- c("SL.mean", "SL.glm", "SL.glmnet.9folds")

#############################################################################################
#############################################################################################
######################### Primary ###########################
#############################################################################################
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_primary.rds") %>% mutate(Areal = as.numeric(Areal),
                                                            Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, cv_folds = 9, seed = 1)
# Sensitive to the high outlier? Somewhat.
do_tmle(dat = dat %>% filter(cp %nin% 282), W = W, A = A, Y = Y, SLL = SLL, seed = 1)
# Sensitive to SL library? no... do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = c("SL.mean", "SL.glm"))
#do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, gcomp = T) #do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, iptw = T)

#############################################################################################
###### Stratifications

# Sex strat - men
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_men.rds") %>% mutate(Areal = as.numeric(Areal),
                                                        Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)
# Sex strat - women
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_women.rds") %>% mutate(Areal = as.numeric(Areal),
                                                          Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)
# Age stratified - Age 11 and under
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_agecat0.rds") %>% mutate(Areal = as.numeric(Areal),
                                                            Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)
# Age stratified - age 12 and over
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_agecat123.rds") %>% mutate(Areal = as.numeric(Areal),
                                                              Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)
# HIV HH strat - positive
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_hh_plus.rds") %>% mutate(Areal = as.numeric(Areal),
                                                            Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)
# Sensitive to the high outlier?
do_tmle(dat = dat %>% filter(cp %nin% 272), W = W, A = A, Y = Y, SLL = SLL, seed = 1)
do_tmle(dat = dat %>% filter(cp %nin% 282), W = W, A = A, Y = Y, SLL = SLL, seed = 1)


# HIV HH strat - negative
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_hh_neg.rds") %>% mutate(Areal = as.numeric(Areal),
                                                           Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)
# Sensitive to the high outlier? A little.
do_tmle(dat = dat %>% filter(cp %nin% 282), W = W, A = A, Y = Y, SLL = SLL, seed = 1)


# HIV people strat - positive
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_hiv_pos_individual.rds") %>% mutate(Areal = as.numeric(Areal),
                                                            Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)
# HIV people strat - negative
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_hiv_neg_individual.rds") %>% mutate(Areal = as.numeric(Areal),
                                                           Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)




#############################################################################################
#############################################################################################
######################### unadjusted Stage2 ###########################
#############################################################################################
# Adjusted Yc, but unadjusted second stage
dat <- readRDS("data/_Yc_primary.rds") %>% mutate(Areal = as.numeric(Areal),
                                                  Afake = as.numeric(Afake),
                                                  U = 1)
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)

#############################################################################################
###### Stratifications
# Sex strat - men
dat <- readRDS("data/_Yc_men.rds") %>% mutate(Areal = as.numeric(Areal),
                                                  Afake = as.numeric(Afake),
                                                  U = 1)
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)
# Sex strat - women
dat <- readRDS("data/_Yc_women.rds") %>% mutate(Areal = as.numeric(Areal),
                                              Afake = as.numeric(Afake),
                                              U = 1)
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)
# Age stratified - Age 11 and under
dat <- readRDS("data/_Yc_agecat0.rds") %>% mutate(Areal = as.numeric(Areal),
                                              Afake = as.numeric(Afake),
                                              U = 1)
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)
# Age stratified - age 12 and over
dat <- readRDS("data/_Yc_agecat123.rds") %>% mutate(Areal = as.numeric(Areal),
                                              Afake = as.numeric(Afake),
                                              U = 1)
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)
# HIV HH strat - positive
dat <- readRDS("data/_Yc_hh_plus.rds") %>% mutate(Areal = as.numeric(Areal),
                                              Afake = as.numeric(Afake),
                                              U = 1)
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)
# HIV HH strat - negative
dat <- readRDS("data/_Yc_hh_neg.rds") %>% mutate(Areal = as.numeric(Areal),
                                              Afake = as.numeric(Afake),
                                              U = 1)
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)
# HIV people strat - positive
dat <- readRDS("data/_Yc_hiv_pos_individual.rds") %>% mutate(Areal = as.numeric(Areal),
                                              Afake = as.numeric(Afake),
                                              U = 1)
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)
# HIV people strat - negative
dat <- readRDS("data/_Yc_hiv_neg_individual.rds") %>% mutate(Areal = as.numeric(Areal),
                                              Afake = as.numeric(Afake),
                                              U = 1)
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)

                                                                       



#############################################################################################
#############################################################################################
######################### unadjusted Stage1 ###########################
#############################################################################################
# Unadjusted Yc, but adjusted second stage
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_unadjusted.rds") %>% mutate(Areal = as.numeric(Areal),
                                                               Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)

# men
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_unadjusted_men.rds") %>% mutate(Areal = as.numeric(Areal),
                                                               Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)

# women
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_unadjusted_women.rds") %>% mutate(Areal = as.numeric(Areal),
                                                                   Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)

# age <12
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_unadjusted_agecat0.rds") %>% mutate(Areal = as.numeric(Areal),
                                                                     Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)

# age 12+
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_unadjusted_agecat123.rds") %>% mutate(Areal = as.numeric(Areal),
                                                                       Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)

# HIV HH strat - positive
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_unadjusted_hh_plus.rds") %>% mutate(Areal = as.numeric(Areal),
                                                                         Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)

# HIV HH strat - neg
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_unadjusted_hh_neg.rds") %>% mutate(Areal = as.numeric(Areal),
                                                                       Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)

# HIV people strat - positive
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_unadjusted_hiv_pos_individual.rds") %>% mutate(Areal = as.numeric(Areal),
                                                                      Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)

# HIV people strat - neg
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_unadjusted_hiv_neg_individual.rds") %>% mutate(Areal = as.numeric(Areal),
                                                                                  Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)



















######################### FULLY unadjusted ###########################
# FULLY unadjusted (Unadjusted Yc, unadjusted second stage)
dat <- readRDS("data/_Yc_unadjusted.rds") %>% mutate(Areal = as.numeric(Areal),
                                                     Afake = as.numeric(Afake),
                                                     U = 1)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)

# Sex strat - men
dat <- readRDS("data/_Yc_unadjusted_men.rds") %>% mutate(Areal = as.numeric(Areal),
                                              Afake = as.numeric(Afake),
                                              U = 1)
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)
# Sex strat - women
dat <- readRDS("data/_Yc_unadjusted_women.rds") %>% mutate(Areal = as.numeric(Areal),
                                                Afake = as.numeric(Afake),
                                                U = 1)
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)
# Age stratified - Age 11 and under
dat <- readRDS("data/_Yc_unadjusted_agecat0.rds") %>% mutate(Areal = as.numeric(Areal),
                                                  Afake = as.numeric(Afake),
                                                  U = 1)
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)
# Age stratified - age 12 and over
dat <- readRDS("data/_Yc_unadjusted_agecat123.rds") %>% mutate(Areal = as.numeric(Areal),
                                                    Afake = as.numeric(Afake),
                                                    U = 1)
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)
# HIV HH strat - positive
dat <- readRDS("data/_Yc_unadjusted_hh_plus.rds") %>% mutate(Areal = as.numeric(Areal),
                                                  Afake = as.numeric(Afake),
                                                  U = 1)
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)
# HIV HH strat - negative
dat <- readRDS("data/_Yc_unadjusted_hh_neg.rds") %>% mutate(Areal = as.numeric(Areal),
                                                 Afake = as.numeric(Afake),
                                                 U = 1)
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)
# HIV people strat - positive
dat <- readRDS("data/_Yc_unadjusted_hiv_pos_individual.rds") %>% mutate(Areal = as.numeric(Areal),
                                                             Afake = as.numeric(Afake),
                                                             U = 1)
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)
# HIV people strat - negative
dat <- readRDS("data/_Yc_unadjusted_hiv_neg_individual.rds") %>% mutate(Areal = as.numeric(Areal),
                                                             Afake = as.numeric(Afake),
                                                             U = 1)
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)






######################### FULLY unadjusted ###########################
# AT COMMUNITY LEVEL THIS TIME
# FULLY unadjusted (Unadjusted Yc, unadjusted second stage)
dat <- readRDS("data/_Yc_unadjusted_comm.rds") %>% mutate(Areal = as.numeric(Areal),
                                                     Afake = as.numeric(Afake),
                                                     U = 1)
do_tmle(dat = dat, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)





















#############################################################################################
# HIV HH strat - positive & young
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_hh_plus_agecat0.rds") %>% mutate(Areal = as.numeric(Areal),
                                                            Afake = as.numeric(Afake)))
plot_outcome(dat = dat, blind = F)
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)
do_tmle(dat = dat %>% filter(cp != "282"), W = W, A = A, Y = Y, SLL = SLL, seed = 1)
# HIV HH strat - positive & old
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_hh_plus_agecat123.rds") %>% mutate(Areal = as.numeric(Areal),
                                                                    Afake = as.numeric(Afake)))
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)
do_tmle(dat = dat %>% filter(cp != "282"), W = W, A = A, Y = Y, SLL = SLL, seed = 1)
# HIV HH strat - negative & young
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_hh_neg_agecat0.rds") %>% mutate(Areal = as.numeric(Areal),
                                                           Afake = as.numeric(Afake)))
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)
# HIV HH strat - negative & old
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_hh_neg_agecat123.rds") %>% mutate(Areal = as.numeric(Areal),
                                                                   Afake = as.numeric(Afake)))
do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, seed = 1)




# Ignoring all clustering... NOT RECOMMENDED!!!
# No adjustment for differential measurement/missingness, but
# adjusting for individual-level predictors of TB.
# Only among people tested twice.
dat <- left_join(readRDS("data/_W_uncl.rds"),
                 readRDS("data/_Yc_uncl.rds") %>% mutate(searchid = as.character(searchid),
                                                         Areal = as.numeric(Areal),
                                                         Afake = as.numeric(Afake))) %>% 
  drop_na
do_tmle(dat = dat, W = W, A = A, Y = "Y_1", SLL = SLL, seed = 1)










# Clustering at community and using adaptive pre-specification
# to select a single adjustment variable:
dat <- left_join(readRDS("data/_comm_APS_W.rds"),
                 readRDS("data/_Yc_comm_APS.rds") %>% mutate(A = as.numeric(Areal),
                                                             Y = Y_final,
                                                             U = 1,
                                                             alpha = 1,
                                                             nIndv = 1)) %>% 
  rename(id = community_num) %>% 
  select(all_of(c(W, "A", "Y", "id", "U", "alpha", "nIndv")))
############################ see "5_test_stage2.R" for more ###########
source("0_Adapt_Functions_cvv.R")
source("0_Stage2_Functions_scaled.R")
Stage2(weighting = "clust",
       goal = "aRR",
       data.input = dat,
       outcome = "Y", do.data.adapt = T,
       clust.adj = c("hiv0", "alc0", "U"),
       break.match = T)







