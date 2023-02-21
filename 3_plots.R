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
                 readRDS("data/_Yc_primary.rds") %>% mutate(Areal = as.numeric(Areal))) %>% mutate(A = Areal)

(results <- do_tmle(dat = dat, W = W, A = A, Y = Y, SLL = SLL, cv_folds = 9, seed = 1, return_unformatted = T))

# Age stratified - Age 11 and under
datY <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_agecat0.rds") %>% mutate(Areal = as.numeric(Areal)))
(resultsY <- do_tmle(dat = datY, W = W, A = A, Y = Y, SLL = SLL, cv_folds = 9, seed = 1, return_unformatted = T))

# Age stratified - age 12 and over
datA <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_agecat123.rds") %>% mutate(Areal = as.numeric(Areal)))
(resultsA <- do_tmle(dat = datA, W = W, A = A, Y = Y, SLL = SLL, cv_folds = 9, seed = 1, return_unformatted = T))

# Sex strat - men
datM <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_men.rds") %>% mutate(Areal = as.numeric(Areal),
                                                        Afake = as.numeric(Afake)))

(resultsM <- do_tmle(dat = datM, W = W, A = A, Y = Y, SLL = SLL, cv_folds = 9, seed = 1, return_unformatted = T))
# Sex strat - women
datW <- left_join(readRDS("data/_parish_W_primary.rds"),
                 readRDS("data/_Yc_women.rds") %>% mutate(Areal = as.numeric(Areal),
                                                          Afake = as.numeric(Afake)))

(resultsW <- do_tmle(dat = datW, W = W, A = A, Y = Y, SLL = SLL, cv_folds = 9, seed = 1, return_unformatted = T))



sz <- .75
############## by age: kids vs all
(plotdat_2a <- rbind.data.frame(results %>% mutate(age = "Overall"),
                             resultsA %>% mutate(age = "Ages 12+"),
                             resultsY %>% mutate(age = "Ages 5 - 11"), make.row.names = F) %>% 
    rename(arm = psi))

(g_2a <- ggplot(plotdat_2a %>% filter(arm == "RR")) +
    geom_point(aes(x = fct_rev(age), y = pt.est, color = age)) +
    geom_errorbar(aes(x = fct_rev(age), ymin = CI.lo, ymax = CI.hi, color = age),
                  width = .2, size = sz) +
    geom_text(data = plotdat_2a %>% filter(arm == "RR") %>% mutate(pval = paste0(#"one-sided p-value",
                                                                              ifelse(pval > .001, paste0("= ", round(pval, digits = 3)),
                                                                              "< .001"))),
              mapping = aes(x = fct_rev(age), label = paste0(#age, "\n\n
                "aRR: ",
                                                             sprintf('%.2f',round(pt.est, digits = 3)),"\n(95% CI: ",
                                                             sprintf('%.2f',round(CI.lo, digits = 3)), " - ",
                                                             sprintf('%.2f',round(CI.hi, digits = 3)),")\np ", pval),
                            y = pt.est), size = 3.3, hjust = 0,  vjust = 1, nudge_x = .05) +
    geom_hline(aes(yintercept = 1), linetype = 2) +
    labs(#title = "Relative risk of incident TB infection in the intervention vs. control arm of the SEARCH trial, overall and stratified by age",
         #subtitle = "Incident TB infection defined as tuberculin skin test conversion from negative at baseline to positive one year later.\nOne-sided p-values to test the null hypothesis that the SEARCH intervention did not reduce TB.",
         y = "Adjusted risk ratio (aRR)") + #ylim(c(.08, .35)) +
    theme_light() + 
    theme(legend.position = "none",
          axis.line = element_line(linewidth = .5, colour = "gray60"),
          #panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(colour = "black"),
          plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 7),
          panel.grid.minor.x = element_blank(),
          axis.ticks.x = element_blank()))

ggsave("results_plot_2a.png", plot = g_2a,
       width = 8, height = 5,
       units = "in", dpi = 300)


############## by age: kids vs all
(plotdat_2b <- rbind.data.frame(results %>% mutate(sex = " Overall"),
                             resultsM %>% mutate(sex = "Men"),
                             resultsW %>% mutate(sex = "Women"), make.row.names = F) %>% 
    rename(arm = psi))

(g_2b <- ggplot(plotdat_2b %>% filter(arm == "RR")) +
    geom_point(aes(x = sex,#fct_rev(sex),
                   y = pt.est, color = sex)) +#,
    geom_errorbar(aes(x = sex,#fct_rev(sex),
                      ymin = CI.lo, ymax = CI.hi, color = sex),
                  #position = position_dodge(width = -.2),
                  width = .2, size = sz) +
    geom_text(data = plotdat_2b %>% filter(arm == "RR") %>% mutate(pval = paste0(#"one-sided p-value",
      ifelse(pval > .001, paste0("= ", round(pval, digits = 3)),
             "< .001"))),
      mapping = aes(x = sex,#fct_rev(sex),
                    label = paste0(#age, "\n\n
        "aRR: ",
        sprintf('%.2f',round(pt.est, digits = 3)),"\n(95% CI: ",
        sprintf('%.2f',round(CI.lo, digits = 3)), " - ",
        sprintf('%.2f',round(CI.hi, digits = 3)),")\np ", pval),
        y = pt.est), size = 3.3, hjust = 0, vjust = 1,
      nudge_x = .05) +
    geom_hline(aes(yintercept = 1), linetype = 2) +
    labs(#title = "Relative risk of incident TB infection in the intervention vs. control arm of the SEARCH trial, overall and stratified by age",
      #subtitle = "Incident TB infection defined as tuberculin skin test conversion from negative at baseline to positive one year later.\nOne-sided p-values to test the null hypothesis that the SEARCH intervention did not reduce TB.",
      y = "Adjusted risk ratio (aRR)") + 
    #ylim(c(.08, .35)) +
    theme_light() + 
    theme(legend.position = "none",
          axis.line = element_line(linewidth = .5, colour = "gray60"),
          panel.grid.minor.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(colour = "black"),
          plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 7),
          panel.grid.minor.x = element_blank(),
          axis.ticks.x = element_blank()))



ggsave("results_plot_2b.png", plot = g_2b,
       width = 8, height = 5,
       units = "in", dpi = 300)



