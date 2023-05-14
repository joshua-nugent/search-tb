library(tidyverse)
`%nin%` <- negate(`%in%`)
source("0_helper_functions.R")

load("data/outputs-withIntOnly-with kids and move and dead.RData")

tb_data <- read_csv("data/TSTdatabase_numeric_hhid split.csv") %>%
  mutate(cp = substr(hhid, start = 1, stop = 3))

fd <- outputs %>% mutate(cp = substr(hhid, start = 1, stop = 3)) %>% 
  filter(cp %in% tb_data$cp, dead_0 != 1, move_0 != 1, age_0 >= 4, data_flag != 1) %>%
  select(intervention, cp, community_num, hhid, searchid,
         chc_0, tr_0, hiv_0, age_0, #age_0 is one year before baseline HH TB measurements
         mobile_0, mobile_3, adult_0,
         occupation_0, occupation_3,
         nightsHome_0, nightsHome_3,
         stable_0, moAway_0,
         tb_0, sex_0, community_name) %>%
  mutate(S = ifelse(hhid %in% tb_data$hhid, 1, 0))


hhids <- unique(fd$hhid)
htypes <- readRDS("data/house_types.rds")
dat <- left_join(fd, htypes) # merge in house_type

final_dat <- left_join(dat,
                       tb_data %>% # Now merge in TB test data
                         mutate(hhid = as.character(hhid), searchid = paste0(hhid, "-", id)) %>% 
                         select(hhid, searchid, cp, tstpos, tstpos_yr2, tstpos_yr3))

# saveRDS(final_dat, "data/_cleaned_ind_data.rds")





















