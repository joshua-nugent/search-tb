library(tidyverse)
`%nin%` <- negate(`%in%`)

# Get baseline variables from SEARCH phase 1 at village level
load("data/outputs-withIntOnly.RData")

d <- read_csv("data/TSTdatabase_numeric_hhid split.csv") %>%
  mutate(cp = substr(hhid, start = 1, stop = 3),
         Afake = ifelse(community_num %% 2 == 0, 1, 0)) %>%
  mutate(Y = tstpos_yr2)

# strip down to just selected parishes
s1 <- outputs %>% mutate(community_num = substr(hhid, start = 1, stop = 2),
                         cp = substr(hhid, start = 1, stop = 3)) %>% 
  filter(cp %in% d$cp) %>% 
  mutate(hiv_0 = ifelse(is.na(hiv_0), 0, hiv_0),
         alcohol_0 = ifelse(is.na(alcohol_0), 0, alcohol_0)) %>% 
  filter(adult_0, resident_0, !dead_0, !move_0)

risk_job_categories <- c("2","3","4","5","6","7","8","9","10","11","13","14","15","16","17","21","23")

s1sum <- s1 %>% group_by(community_num) %>%
  summarise(tb0 = mean(tb_0 == T),
            hiv0 = mean(hiv_0 == T),
            mob0 = mean(mobile_0 == T),
            alc0 = mean(alcohol_0 == T),
            hr0 = mean(occupation_0 %in% risk_job_categories),
            pphh = n()/length(unique(hhid)),
            nvill = length(unique(substr(hhid, start = 1, stop = 5))))
# saveRDS(s1sum, "data/_comm_APS_W.rds")

















