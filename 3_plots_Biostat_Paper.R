library(tidyverse)
library(knitr)
library(ltmle)
library(glmnet)
library(coefplot)
library(forcats)
source("0_helper_functions.R")

#################################
blind <- F
A <- "Areal"
W <- c("hiv0", "alc0")
Y <- "Y_final"
SLL <- c("SL.mean", "SL.glm", "SL.earth")
######################### Primary ###########################
dat <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_primary.rds") %>% mutate(Areal = as.numeric(Areal)))

(results <- do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, return_unformatted = T, seed = 1))


# Unadjusted
datU <- readRDS("data/_Yc_unadjusted.rds") %>% mutate(Areal = as.numeric(Areal),
                                                     Afake = as.numeric(Afake),
                                                     U = 1)
(resultsU <- do_tmle(dat = datU, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1))

# Adjusted Yc, but unadjusted second stage
dat1 <- readRDS("data/_Yc_primary.rds") %>% mutate(Areal = as.numeric(Areal),
                                                  Afake = as.numeric(Afake),
                                                  U = 1)
results1 <- do_tmle(dat = dat1, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)

# Unadjusted Yc, but adjusted second stage
dat2 <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_unadjusted.rds") %>% mutate(Areal = as.numeric(Areal),
                                                               Afake = as.numeric(Afake)))
results2 <- do_tmle(dat = dat2, W = W, A = A, Y = Y, SLL = SLL, seed = 1)


# AT COMMUNITY LEVEL THIS TIME
# FULLY unadjusted (Unadjusted Yc, unadjusted second stage)
dat3 <- readRDS("data/_Yc_unadjusted_comm.rds") %>% mutate(Areal = as.numeric(Areal),
                                                          Afake = as.numeric(Afake),
                                                          U = 1)
results3 <- do_tmle(dat = dat3, W = "U", A = A, Y = Y, SLL = "SL.mean", seed = 1)


# Clustering at community and using adaptive pre-specification
# to select a single adjustment variable:
dat_aps <- left_join(readRDS("data/_comm_APS_W.rds"),
                 readRDS("data/_Yc_comm_APS.rds") %>% mutate(A = as.numeric(Areal),
                                                             Y = Y_final,
                                                             U = 1,
                                                             alpha = 1,
                                                             nIndv = 1)) %>% 
  rename(id = community_num) %>% 
  select(all_of(c(W, "A", "Y", "id", "U", "alpha", "nIndv")))



source("0_Adapt_Functions_cvv.R")
source("0_Stage2_Functions_scaled.R")
set.seed(1)
s2_aps <- Stage2(weighting = "clust",
       goal = "aRR",
       data.input = dat_aps,
       outcome = "Y", do.data.adapt = T,
       clust.adj = c("hiv0", "alc0", "U"),
       break.match = T)
s2_aps
(results_aps <- cbind.data.frame(psi = "RR",
                                pt.est = s2_aps$Effect.est,
                                CI.lo = s2_aps$Effect.CI.lo,
                                CI.hi = s2_aps$Effect.CI.hi,
                                se = sqrt(s2_aps$var.hat),
                                pval = s2_aps$pval))






############## VERTICAL PLOT
(plotdat2 <- rbind.data.frame(results %>% mutate(age = " Stage 1 & 2\nadjusted\n (parish)"),
                             resultsU %>% mutate(age = "Unadjusted\n\n (parish)"),
                             results1 %>% mutate(age = "Stage 1\nadjustment only\n(parish)"),
                             results2 %>% mutate(age = "Stage 2\nadjustment only\n(parish)"),
                             #results3 %>% mutate(age = "Unadjusted\n\n(community level)"),
                             results_aps %>% mutate(age = " Stage 1 & 2\nadjusted\n(community)"),
                             #results4 %>% mutate(age = " Stage 1 & 2\nadjusted\n(community level)"),
                             make.row.names = F) %>% 
    rename(arm = psi) %>%
  filter(arm == "RR") %>%
  select(-se) %>% 
    mutate(pval = ifelse(pval > .001, paste0("= ", round(pval, digits = 2)),"< .001")) %>% 
  mutate(labz = paste0("aRR: ", sprintf('%.2f',round(pt.est, digits = 2)),"\n(95% CI: ",
                         sprintf('%.2f',round(CI.lo, digits = 2)), " - ",
                         sprintf('%.2f',round(CI.hi, digits = 2)),")\np ", pval)) %>% 
  mutate(typ = 1:5))
  
  
(g <- ggplot(plotdat2) +
    geom_hline(aes(yintercept = 1), linetype = 1, color = "grey80") +
    #geom_point(aes(x = age, y = pt.est, color = (typ==1))) + 
    geom_point(aes(x = age, y = pt.est), color = "grey20") + 
    #geom_errorbar(aes(x = age, ymin = CI.lo, ymax = CI.hi, color = (typ==1)), linetype = 1, width = .2) +
    geom_errorbar(aes(x = age, ymin = CI.lo, ymax = CI.hi, linetype = (typ!=1)), color = "grey20",# linetype = 1,
                  width = .2) +
    geom_text(mapping = aes(x = age, label = labz, y = pt.est),
              size = 2.5, hjust = 0, nudge_x = .05, vjust = 1) +
    labs(title = "Estimated effects on the one-year incidence of TB in SEARCH-TB by method",
         y = "Adjusted Relative Risk (aRR)") + 
    theme_light() + 
    scale_x_discrete(expand=c(0, 1))+
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_text(colour = "black", size = 8),
          plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 7),
          axis.text = element_text(colour = "black", size = 8),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          panel.grid.minor.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()))
  
ggsave("_plots/results_plot_BW.png", plot = g,
       width = 7, height = 4.5,
       units = "in", dpi = 300)




