library(tidyverse)
`%nin%` <- negate(`%in%`)

merge_muyembe_1 <- c("261", "262", "263", "264", "265")
merge_muyembe_2 <- c("266", "267", "268", "269")

# Get baseline variables from SEARCH phase 1 at village level
load("data/outputs-withIntOnly.RData")

d <- read_csv("data/TSTdatabase_numeric_hhid split.csv") %>%
  mutate(cp = substr(hhid, start = 1, stop = 3),
         cp = ifelse(cp %in% merge_muyembe_1, "261",
                     ifelse(cp %in% merge_muyembe_2, "266", cp)),
         Afake = ifelse(community_num %% 2 == 0, 1, 0)) %>%
  mutate(Y = tstpos_yr2)


# strip down to just selected parishes
s1 <- outputs %>% mutate(cp = substr(hhid, start = 1, stop = 3),
                         cp = ifelse(cp %in% merge_muyembe_1, "261",
                                     ifelse(cp %in% merge_muyembe_2, "266", cp))) %>% 
  filter(cp %in% d$cp) %>% 
  mutate(hiv_0 = ifelse(is.na(hiv_0), 0, hiv_0),
         alcohol_0 = ifelse(is.na(alcohol_0), 0, alcohol_0)) %>% 
  filter(adult_0, resident_0, !dead_0, !move_0)


risk_job_categories <- c("2","3","4","5","6","7","8","9","10","11","13","14","15","16","17","21","23")

s1sum <- s1 %>% group_by(cp) %>%
  summarise(tb0 = mean(tb_0 == T),
            hiv0 = mean(hiv_0 == T),
            mob0 = mean(mobile_0 == T),
            alc0 = mean(alcohol_0 == T),
            hr0 = mean(occupation_0 %in% risk_job_categories),
            pphh = n()/length(unique(hhid)),
            nvill = length(unique(substr(hhid, start = 1, stop = 5))))

# saveRDS(s1sum, "data/_parish_W_primary.rds")

############# UNclustered

s1sum_uncl <- s1 %>%
  mutate(riskjob = occupation_0 %in% risk_job_categories) %>% 
  dplyr::select(tb_0, hiv_0, mobile_0, alcohol_0, riskjob, searchid) %>% 
  rename(tb0 = tb_0, hiv0 = hiv_0, alc0 = alcohol_0)

# saveRDS(s1sum_uncl, "data/_W_uncl.rds")

















