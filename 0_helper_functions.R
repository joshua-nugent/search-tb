`%nin%` <- negate(`%in%`)

get_As <- function(d0, output){
  as <- d0 %>% group_by(Afake, cp) %>% summarise(Afake = as.logical(first(Afake)),
                                                 cp = first(cp),
                                                 community_num = first(community_num),
                                                 community_name = first(community_name),
                                                 Areal = first(intervention))
  out <- left_join(output, as)
  return(out)
}

get_As_comm <- function(d0, output){
  as <- d0 %>% group_by(Afake, community_num) %>% summarise(Afake = as.logical(first(Afake)),
                                                 community_num = first(community_num),
                                                 community_name = first(community_name),
                                                 Areal = first(intervention))
  out <- left_join(output, as)
  return(out)
}


get_Yc <- function(parish, d0, W, A, Y, SLL, ignore_S = F,
                   age_strat = "all", column_names){
  output <- data.frame(matrix(NA, ncol = length(column_names), nrow = 1))
  colnames(output) <- column_names
  
  v0 <- d0 %>% filter(cp %in% parish) %>% select(all_of(c(W, A, Y)))
  v1 <- v0 %>% drop_na
  
  output[1,"cp"] <- parish
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

  if(ignore_S){
    data.temp0 <- d0 %>% filter(cp %in% parish) %>% 
      select(all_of(c(W, "Delta_0", "Y_0",
                      "Delta_1", "Y_1")))
  } else {
    data.temp0 <- d0 %>% filter(cp %in% parish) %>% 
      select(all_of(c(W, "S_1_Delta_0", "Y_0",
                      "Delta_1", "Y_1")))
  }
  data.temp <- data.temp0 %>% drop_na
  
  output[1,"number_tested_FU"] <- sum(data.temp$Delta_1)
  
  print(paste("n =", nrow(data.temp)))
  abar <- matrix(nrow=nrow(data.temp), ncol=2)
  
  if(ignore_S){
    abar[,1] <- 1
    abar[,2] <- data.temp$Y_0 == 0 & data.temp$Delta_0
    data.temp <- data.temp %>%
      select(all_of(c(W, 
                      "Delta_0", "Y_0",
                      "Delta_1", "Y_1")))
    
    est.temp <- 
      ltmle(data = data.temp,
            Anodes = c('Delta_0', 'Delta_1'),
            Lnodes = 'Y_0', Ynodes = 'Y_1',
            abar = abar, SL.library = SLL, estimate.time = F,
            stratify = T, variance.method = 'ic')#))
  } else {
    abar[,1] <- 1
    abar[,2] <- data.temp$Y_0 == 0 & data.temp$S_1_Delta_0
    data.temp <- data.temp %>%
      select(all_of(c(W, 
                      "S_1_Delta_0", "Y_0",
                      "Delta_1", "Y_1")))
    
    est.temp <- #suppressWarnings(suppressMessages(
      ltmle(data = data.temp,
            Anodes = c('S_1_Delta_0', 'Delta_1'),
            Lnodes = 'Y_0', Ynodes = 'Y_1',
            abar = abar, SL.library = SLL, estimate.time = F,
            stratify = T, variance.method = 'ic')#))
    
  }
  print(est.temp$fit)
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





# get.var.delta - function to get inference via the delta method
# 	assumes inputed estimators are asymptotically linear
#		i.e. written in first order as an empircal mean of an influence curve (IC)
#	input:  point estimates (mu1, mu0), corresponding influence curves (IC1, IC0)
#		significance level
#	output: point estimate, var, wald-type CI 

get.var.delta <- function(mu1, mu0=NULL, IC1, IC0=NULL, alpha=0.05){
  
  mu1<- unlist(mu1)
  
  if(is.null(mu0)){ 
    # if single TMLE 
    psi <- mu1
    IC <- IC1
    log <- F
    
  } else { 
    # if ratio of TMLEs (i.e. target = psi/psi0)
    mu0 <- unlist(mu0)
    # get inference via the delta method on log scale
    psi <- log(mu1/mu0)
    IC <- 1/mu1*(IC1) - 1/mu0*IC0
    log <- T
  }
  
  # variance of asy lin est is var(IC)/n
  var<- var(IC)/length(IC)
  # testing and CI	
  cutoff <- qnorm(alpha/2, lower.tail=F)
  se <- sqrt(var)
  CI.lo <- psi - cutoff*se
  CI.hi <- psi + cutoff*se
  
  if(log){
    est <- data.frame(pt=exp(psi), CI.lo=exp(CI.lo), CI.hi=exp(CI.hi) ) 
  }else{
    est <- data.frame(pt=psi, CI.lo=CI.lo, CI.hi=CI.hi) 
  }
  
  list(est=est, IC=IC)
}


get_house_type <- function(hh_id){ # fd is "full data" set
  temp <- fd %>% filter(as.character(hhid) == as.character(hh_id),
                        adult_0 == T) %>% 
    mutate(chc_pos = !is.na(hiv_0) & chc_0 == TRUE & hiv_0 == T, # HType 1
           trk_pos = !is.na(hiv_0) & tr_0 == TRUE & hiv_0 == T,  # HType 2
           chc_test = !is.na(hiv_0) & chc_0 == TRUE,
           trk_test = !is.na(hiv_0) & tr_0 == TRUE)
  if(sum(temp$chc_pos) > 0){return("1")}
  if(sum(temp$trk_pos) > 0){return("2")}
  if(sum(temp$hiv_0, na.rm = T) == 0 & sum(temp$chc_test, na.rm = T) > 0){
    #print(temp)
    return("3")}
  if(sum(temp$hiv_0, na.rm = T) == 0 & sum(temp$chc_test, na.rm = T) == 0 & sum(temp$trk_test) > 0){
    return("4")}
  return("5")
}


do_tmle <- function(dat, W, A, Y, SLL, gcomp = F, iptw = F,
                    return_unformatted = T, one_sided = T, cv_folds = 9,
                    seed = NA){
  if(is.na(seed)){
    seed <- floor(runif(n = 1, min = 1, max = 9999999))
  }
  set.seed(seed)
  print("print(summary(dat[[min_g_q10_denom]]))")
  print(summary(dat[["min_g_q10_denom"]]))
  print(paste0("DENOM n min = ", min(dat[["number_sampled_and_tested_BL"]])))
  print(summary(dat[["Y_0_0"]]))
  print("print(summary(dat[[min_g_q10_num]]))")
  print(summary(dat[["min_g_q10_num"]]))
  print(paste0("NUM n min = ", min(dat[["number_tested_FU"]])))
  print(summary(dat[["Y_1_1"]]))
  mod1 <- suppressWarnings(suppressMessages(ltmle(data = dat %>% dplyr::select(all_of(c(W, A, Y))),
                Anodes = A, Ynodes = Y, abar = 1,
                variance.method = "ic", gcomp = gcomp, iptw.only = iptw,
                SL.library = SLL, SL.cvControl = list(V = cv_folds))))
  print(mod1$fit)
  print(summary(mod1$cum.g.unbounded))
  mod0 <- suppressWarnings(suppressMessages(ltmle(data = dat %>% dplyr::select(all_of(c(W, A, Y))),
                Anodes = A, Ynodes = Y, abar = 0,
                variance.method = "ic", gcomp = gcomp, iptw.only = iptw,
                SL.library = SLL, SL.cvControl = list(V = cv_folds))))
  the_end_result <- get.inference_31(txt = mod1, con = mod0, gcomp, iptw, one_sided = one_sided)
  if(return_unformatted){
    return(the_end_result)
  }
}



get.inference_31 <- function(txt, con, gcomp = F, iptw = F, alternate_null_negative = T, one_sided = T){
  if(gcomp){
    psi.1 <- txt$estimates["gcomp"]
    psi.0 <- con$estimates["gcomp"]
    IC.1 <- txt$IC$gcomp
    IC.0 <- con$IC$gcomp
  } else if(iptw){
    psi.1 <- txt$estimates["iptw"]
    psi.0 <- con$estimates["iptw"]
    IC.1 <- txt$IC$iptw
    IC.0 <- con$IC$iptw
  } else {
    psi.1 <- txt$estimates["tmle"]
    psi.0 <- con$estimates["tmle"]
    IC.1 <- txt$IC$tmle
    IC.0 <- con$IC$tmle
  }
  IC.diff <- IC.1 - IC.0
  # going after aRR, then get IC estimate on log scale
  #	i.e. Delta method for log(aRR) = log(R1) - log(R0)
  IC.ratio <- 1/psi.1*IC.1 - 1/psi.0*IC.0
  txt <- get.CI_31(psi.hat = psi.1, IC = IC.1, start = "intervention", do.relative = F, one_sided = one_sided)
  con <- get.CI_31(psi.hat = psi.0, IC = IC.0, start = "control", do.relative = F, one_sided = one_sided)
  RD <- get.CI_31(psi.hat = (psi.1 - psi.0), IC = IC.diff, one_sided = one_sided,
                  do.relative = F, alternate_null_negative, start = "RD")
  RR <- get.CI_31(psi.hat = log(psi.1 / psi.0), IC = IC.ratio, one_sided = one_sided,
                  do.relative = T, alternate_null_negative, start = "RR")
  
  out <- rbind.data.frame(txt, con,
               RR) %>% rownames_to_column("psi")
  return(out)
}

get.CI_31 <- function(psi.hat, IC, do.relative = T, start, one_sided = T, alternate_null_negative = T){
  J <- length(IC)
  var.IC <- var(IC) / J
  # cutoff based on t-dist for testing and CI	
  cutoff <- qt(0.05 / 2, df = (J - 2), lower.tail = F)
  # standard error (square root of the variance)
  se <- sqrt(var.IC)
  # test statistic (if goal=aRR then on the transformed scale)
  tstat <- psi.hat / se
  
  if(one_sided){
    pval <- pt(tstat, df = (J - 2), lower.tail = alternate_null_negative)
  } else {
    pval <- 2*pt(abs(tstat), df = (J - 2), lower.tail = F)
  }
  # 95% confidence interval
  CI.lo <- psi.hat - cutoff*se
  CI.hi <- psi.hat + cutoff*se
  if(do.relative){
    psi.hat <- exp(psi.hat)
    CI.lo <- exp(CI.lo)
    CI.hi <- exp(CI.hi)
  }
  out <- cbind.data.frame(#psi = start,
                          pt.est = psi.hat,
                          CI.lo = CI.lo, CI.hi = CI.hi,
                          se = se,
                          pval = pval)
  rownames(out) <- start     #colnames(out) <- paste(start, colnames(out), sep=".")
  return(out)
}



plot_outcome <- function(dat, blind){
  if(blind){
    dat <- dat %>% mutate(A = Afake,
                          grp_mean = ifelse(Afake == 0, mean(dat$Y_final[dat$Afake == 0]),
                                            mean(dat$Y_final[dat$Afake == 1])))
  } else {
    dat <- dat %>% mutate(A = Areal,
                          grp_mean = ifelse(Areal == 0, mean(dat$Y_final[dat$Areal == 0]),
                                            mean(dat$Y_final[dat$Areal == 1])))
  }
  ggplot(data = dat %>% mutate(A = ifelse(A == 1, "intervention", "control"))) +
    geom_errorbar(aes(x = cp, ymin = Y_final_lo, ymax = Y_final_hi,
                      color = community_num)) +
    geom_point(aes(x = cp, y = Y_final, color = community_num)) +
    labs(y = "New TB infection estimate", x = "parish code (color = community)") +
    theme_light() + 
    theme(legend.position = "none") +
    facet_grid(.~A, drop = T, scales = "free_x")
}



get_summ_stats <- function(parish, d0, age_strat = "all", column_names){
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
  v1 <- d0 %>% filter(cp %in% parish)
  
  output[1,"A"] <- v1$intervention[1]
  output[1,"cp"] <- parish
  output[1,"population"] <- nrow(v1)
  output[1,"number_sampled"] <- sum(v1$S)
  output[1,"number_sampled_and_tested_BL"] <- sum(v1$S_1_Delta_0)
  

  data.temp0 <- d0 %>% filter(cp %in% parish) %>% 
    select(all_of(c("S_1_Delta_0", "Y_0", "risky_job_3",
                    "mobile_3", "Delta_1", "Y_1")))
  data.temp <- data.temp0 %>% drop_na
  output[1,"number_tested_FU"] <- sum(data.temp$Delta_1)
  output[1,"number_tested_FU_neg_BL"] <- sum(data.temp$Delta_1 & data.temp$Y_0 == 0) 
  return(output)
}





