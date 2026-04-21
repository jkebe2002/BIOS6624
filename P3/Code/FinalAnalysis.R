library(tidyverse)
library(gt)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(gridExtra)


framdat_raw <- read_csv("~/Downloads/frmgham2.csv")
framdat <- read_csv("~/Downloads/frmgham2.csv")
library(gtsummary)


#Subjects with TIMESTRK 
#no subjects who already had a stroke
framdat <- framdat_raw %>%
  filter(TIMESTRK != 0 )

framdat <- framdat %>%
  filter(PERIOD==1 )

framdat <- framdat %>%
  mutate(time_obs = pmin(TIMESTRK, ceiling(10*365.25))) %>%
  select(-c(TIME,HDLC,LDLC,PERIOD,PREVSTRK)) %>%
  #censor at 10 years
  mutate(stroke_10yr = if_else(TIMESTRK > ceiling(10*365.25), 0, STROKE )
  )

framdat_m <- framdat %>%
  filter(SEX==1)
framdat_f <- framdat %>%
  filter(SEX==2)


# Filter to keep only complete rows
framdat_m <- framdat_m %>%
  filter(complete.cases(time_obs, STROKE, AGE, SYSBP, DIABETES, 
                        BPMEDS, PREVCHD, CURSMOKE, TOTCHOL, BMI))


# Filter to keep only complete rows
framdat_f <- framdat_f %>%
  filter(complete.cases(time_obs, STROKE, AGE, SYSBP, DIABETES, 
                        BPMEDS, PREVCHD, CURSMOKE, TOTCHOL, BMI))

# Some variables are of more interest 
# Age, diabetes, blood pressure (systolic).
# Several whether or not to include
# coronary heart disease, blood pressure meds (Y/N), smoking status, Total cholesterol, BMI.
# - baseline values for all.
library(survival)

################################
# Kaplan Meier Curves for Vars #
################################
library(survival)
library(survminer)
surv1<-survfit(Surv(time_obs, stroke_10yr)~AGE, data=framdat_m)
ggsurvplot(surv1, conf.int=FALSE,legend="none",
           xlim=c(0, max(framdat_m$time_obs)))
surv2<-survfit(Surv(time_obs, stroke_10yr)~SYSBP, data=framdat_m)
ggsurvplot(surv2, conf.int=FALSE,legend="none",
           xlim=c(0, max(framdat_m$time_obs)))
surv3<-survfit(Surv(time_obs, stroke_10yr)~DIABETES, data=framdat_m)
ggsurvplot(surv3, conf.int=FALSE,legend="none",
           xlim=c(0, max(framdat_m$time_obs)))
surv4<-survfit(Surv(time_obs, stroke_10yr)~CURSMOKE, data=framdat_m)
ggsurvplot(surv4, conf.int=FALSE,legend="none",
           xlim=c(0, max(framdat_m$time_obs)))
surv5<-survfit(Surv(time_obs, stroke_10yr)~AGE, data=framdat_f)
ggsurvplot(surv5, conf.int=FALSE,legend="none",
           xlim=c(0, max(framdat_m$time_obs)))
surv6<-survfit(Surv(time_obs, stroke_10yr)~SYSBP, data=framdat_f)
ggsurvplot(surv6, conf.int=FALSE,legend="none",
           xlim=c(0, max(framdat_m$time_obs)))
surv7<-survfit(Surv(time_obs, stroke_10yr)~DIABETES, data=framdat_f)
ggsurvplot(surv7, conf.int=FALSE,legend="none",
           xlim=c(0, max(framdat_m$time_obs)))
surv8<-survfit(Surv(time_obs, stroke_10yr)~CURSMOKE, data=framdat_f)
ggsurvplot(surv8, conf.int=FALSE,legend="none",
           xlim=c(0, max(framdat_m$time_obs)))




#model selection
full_model_m <- coxph(Surv(time_obs,stroke_10yr) ~ AGE+SYSBP+DIABETES+BPMEDS+PREVCHD+CURSMOKE+TOTCHOL+BMI, 
                      data = framdat_m)
null_model_m <- coxph(Surv(time_obs,stroke_10yr) ~ AGE+SYSBP+DIABETES, data = framdat_m)


# Selection based on lowest AIC
final_model_m <- step(
  full_model_m, 
  scope = list(lower = ~ AGE+SYSBP+DIABETES, upper = formula(full_model_m)), 
  direction = "backward", trace=0
)



full_model_f <- coxph(Surv(time_obs,stroke_10yr ) ~ AGE+SYSBP+DIABETES+BPMEDS+PREVCHD+CURSMOKE+TOTCHOL+BMI, 
                      data = framdat_f)
null_model_f <- coxph(Surv(time_obs,stroke_10yr ) ~ AGE+SYSBP+DIABETES, data = framdat_f)





final_model_f <- step(
  full_model_f, 
  scope = list(lower = ~ AGE+SYSBP+DIABETES, upper = formula(full_model_f)), 
  direction = "backward", trace = 0
)


#final_model_f_smoking <- coxph(Surv(time_obs,stroke_10yr ) ~ AGE+SYSBP+DIABETES+CURSMOKE, 
#                      data = framdat_f)






###################################################
# SURVIVAL CURVES FOR MEN AND WOMEN RISK PROFILES #
###################################################




newdata2<-data.frame(TRTMT=c(0,1), FUNCTCLS=0, DIGUSE=0, DIURET=0,
                     CHESTX_greater_55=0, EJF_PER_above_25=0)


#Males, SBP at median and 3rd quartile, diabetes
male_risk_profiles_SBP_diab_40 <- data.frame(AGE = 40, SYSBP = c(128.5,128.5, 160, 160), DIABETES = c(0,1,0,1), CURSMOKE = 0)

male_risk_profiles_SBP_diab_50 <- data.frame(AGE = 50, SYSBP = c(128.5,128.5, 160, 160), DIABETES = c(0,1,0,1), CURSMOKE = 0)

male_risk_profiles_SBP_diab_60 <- data.frame(AGE = 60, SYSBP = c(128.5,128.5, 160, 160), DIABETES = c(0,1,0,1), CURSMOKE = 0)

##############################
# SMOKING RISK PROFILES MALE #
##############################
male_risk_profiles_SBP_diab_40_sm <- data.frame(AGE = 40, SYSBP = c(128.5,128.5, 160, 160), DIABETES = c(0,1,0,1), CURSMOKE = 1)

male_risk_profiles_SBP_diab_50_sm <- data.frame(AGE = 50, SYSBP = c(128.5,128.5, 160, 160), DIABETES = c(0,1,0,1), CURSMOKE = 1)

male_risk_profiles_SBP_diab_60_sm <- data.frame(AGE = 60, SYSBP = c(128.5,128.5, 160, 160), DIABETES = c(0,1,0,1), CURSMOKE = 1)


#############################
# FEMALE RISK PROFILES      #
#############################
female_risk_profiles_SBP_diab_40 <- data.frame(AGE = 40, SYSBP = c(128.5,128.5, 160, 160), DIABETES = c(0,1,0,1), CURSMOKE = 0)

female_risk_profiles_SBP_diab_50 <- data.frame(AGE = 50, SYSBP = c(128.5,128.5, 160, 160), DIABETES = c(0,1,0,1), CURSMOKE = 0)

female_risk_profiles_SBP_diab_60 <- data.frame(AGE = 60, SYSBP = c(128.5,128.5, 160, 160), DIABETES = c(0,1,0,1), CURSMOKE = 0)

#smoking = 1
female_risk_profiles_SBP_diab_40_sm <- data.frame(AGE = 40, SYSBP = c(128.5,128.5, 160, 160), DIABETES = c(0,1,0,1), CURSMOKE = 1)

female_risk_profiles_SBP_diab_50_sm <- data.frame(AGE = 50, SYSBP = c(128.5,128.5, 160, 160), DIABETES = c(0,1,0,1), CURSMOKE = 1)

female_risk_profiles_SBP_diab_60_sm <- data.frame(AGE = 60, SYSBP = c(128.5,128.5, 160, 160), DIABETES = c(0,1,0,1), CURSMOKE = 1)



library(survminer)
library(survival)
library(cowplot)

xlim_max <- max(survfit(final_model_m)$time)
legend_labs <- c("Median SBP, No Diabetes", "Median SBP, Diabetes", 
                 "High SBP, No Diabetes", "High BP, Diabetes")

# ── Build ggsurvplots ──────────────────────────────────────────────────────────
plots <- list(
  ggsurvplot(survfit(final_model_m, newdata = male_risk_profiles_SBP_diab_40), 
             data = framdat_m, censor = F, conf.int = F, legend = "right",
             xlim = c(0, xlim_max), legend.labs = legend_labs,
             title = "Age 40",
             ggtheme = theme(plot.title = element_text(hjust = 0.5))),
  ggsurvplot(survfit(final_model_m, newdata = male_risk_profiles_SBP_diab_50), 
             data = framdat_m, censor = F, conf.int = F, legend = "right",
             xlim = c(0, xlim_max), legend.labs = legend_labs,
             title = "Age 50",
             ggtheme = theme(plot.title = element_text(hjust = 0.5))),
  ggsurvplot(survfit(final_model_m, newdata = male_risk_profiles_SBP_diab_60), 
             data = framdat_m, censor = F, conf.int = F, legend = "right",
             xlim = c(0, xlim_max), legend.labs = legend_labs,
             title = "Age 60",
             ggtheme = theme(plot.title = element_text(hjust = 0.5)))
)

plots_sm <- list(
  ggsurvplot(survfit(final_model_m, newdata = male_risk_profiles_SBP_diab_40_sm), 
             data = framdat_m, censor = F, conf.int = F, legend = "right",
             xlim = c(0, xlim_max), legend.labs = legend_labs,
             title = "Age 40",
             ggtheme = theme(plot.title = element_text(hjust = 0.5))),
  ggsurvplot(survfit(final_model_m, newdata = male_risk_profiles_SBP_diab_50_sm), 
             data = framdat_m, censor = F, conf.int = F, legend = "right",
             xlim = c(0, xlim_max), legend.labs = legend_labs,
             title = "Age 50",
             ggtheme = theme(plot.title = element_text(hjust = 0.5))),
  ggsurvplot(survfit(final_model_m, newdata = male_risk_profiles_SBP_diab_60_sm), 
             data = framdat_m, censor = F, conf.int = F, legend = "right",
             xlim = c(0, xlim_max), legend.labs = legend_labs,
             title = "Age 60",
             ggtheme = theme(plot.title = element_text(hjust = 0.5)))
)

plots_f <- list(
  ggsurvplot(survfit(final_model_f, newdata = female_risk_profiles_SBP_diab_40), 
             data = framdat_f, censor = F, conf.int = F, legend = "right",
             xlim = c(0, xlim_max), legend.labs = legend_labs,
             title = "Age 40",
             ggtheme = theme(plot.title = element_text(hjust = 0.5))),
  ggsurvplot(survfit(final_model_f, newdata = female_risk_profiles_SBP_diab_50), 
             data = framdat_f, censor = F, conf.int = F, legend = "right",
             xlim = c(0, xlim_max), legend.labs = legend_labs,
             title = "Age 50",
             ggtheme = theme(plot.title = element_text(hjust = 0.5))),
  ggsurvplot(survfit(final_model_f, newdata = female_risk_profiles_SBP_diab_60), 
             data = framdat_f, censor = F, conf.int = F, legend = "right",
             xlim = c(0, xlim_max), legend.labs = legend_labs,
             title = "Age 60",
             ggtheme = theme(plot.title = element_text(hjust = 0.5)))
)



# ── Extract plots and legend ───────────────────────────────────────────────────
p1 <- plots[[1]]$plot + theme(legend.position = "none")
p2 <- plots[[2]]$plot + theme(legend.position = "none")
p3 <- plots[[3]]$plot + theme(legend.position = "right")

p1_sm <- plots_sm[[1]]$plot + theme(legend.position = "none")
p2_sm <- plots_sm[[2]]$plot + theme(legend.position = "none")
p3_sm <- plots_sm[[3]]$plot + theme(legend.position = "right")

p1_f <- plots_f[[1]]$plot + theme(legend.position = "none")
p2_f <- plots_f[[2]]$plot + theme(legend.position = "none")
p3_f <- plots_f[[3]]$plot + theme(legend.position = "right")

clean <- function(p, legend = "none") {
  p + 
    theme_classic() +
    theme(
      legend.position  = legend,
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

p1    <- clean(plots[[1]]$plot)
p2    <- clean(plots[[2]]$plot)
p3    <- clean(plots[[3]]$plot, legend = "right")

p1_sm <- clean(plots_sm[[1]]$plot)
p2_sm <- clean(plots_sm[[2]]$plot)
p3_sm <- clean(plots_sm[[3]]$plot, legend = "right")

p1_f  <- clean(plots_f[[1]]$plot)
p2_f  <- clean(plots_f[[2]]$plot)
p3_f  <- clean(plots_f[[3]]$plot, legend = "right")

# extract legend from cleaned p3
shared_legend <- cowplot::get_plot_component(p3, "guide-box", return_all = TRUE)[[1]]

p3_no_legend    <- clean(plots[[3]]$plot)
p3_no_legend_sm <- clean(plots_sm[[3]]$plot)
p3_no_legend_f  <- clean(plots_f[[3]]$plot)



# Use get_plot_component() instead of get_legend() for newer cowplot
shared_legend <- cowplot::get_plot_component(p3, "guide-box", return_all = TRUE)[[1]]

p3_no_legend    <- p3    + theme(legend.position = "none")
p3_no_legend_sm <- p3_sm + theme(legend.position = "none")
p3_no_legend_f  <- p3_f  + theme(legend.position = "none")


# ── Row labels ─────────────────────────────────────────────────────────────────
label_m_ns <- ggdraw() + draw_label("Male\nNon-Sm", fontface = "bold", 
                                    hjust = 0.5, vjust = 0.5, size = 9)
label_m_sm <- ggdraw() + draw_label("Male\nSm",     fontface = "bold", 
                                    hjust = 0.5, vjust = 0.5, size = 9)
label_f    <- ggdraw() + draw_label("Female",            fontface = "bold", 
                                    hjust = 0.5, vjust = 0.5, size = 9)



# ── 3x3 grid: label | p1 | p2 | p3 | legend ───────────────────────────────────
row_m_ns <- plot_grid(label_m_ns, p1,    p2,    p3_no_legend,    nrow = 1, rel_widths = c(0.18, 1, 1, 1))
row_m_sm <- plot_grid(label_m_sm, p1_sm, p2_sm, p3_no_legend_sm, nrow = 1, rel_widths = c(0.18, 1, 1, 1))
row_f    <- plot_grid(label_f,    p1_f,  p2_f,  p3_no_legend_f,  nrow = 1, rel_widths = c(0.18, 1, 1, 1))

# Stack rows, add shared legend to the right
grid_4x3 <- plot_grid(row_m_ns, row_m_sm, row_f,nrow = 3)

plot_final_4x3 <- plot_grid(
  grid_4x3, shared_legend,
  nrow = 1,
  rel_widths = c(1, 0.18)
) + theme(plot.background = element_rect(fill = "white", color = NA))

#plot_final_4x3



# Table of risk profiles

# Survival curves M non smoker, M smoker, F 

# Table 1

# We can include one more table or graph...



########################
# CONFIDENCE INTERVALS #
########################

confint_coxph <- function(model, L, dat) {
  #returns exponentiated Wald CI's for estimate for different risk profiles
  # provides risk profile relative to nonsmoker, no diabetes, median blood pressure
  # L_diff <- L - t(c(median(dat$AGE), median(dat$SYSBP), 0, 0))
  # est  <- L_diff %*% model$coefficients
  # lower <- est - 1.96*sqrt(L_diff %*% model$var %*% t(L_diff))
  # upper <- est + 1.96*sqrt(L_diff %*% model$var %*% t(L_diff) )
  # return(exp(c(est, lower, upper)))
  
  if (length(L) == 4) {
    l_df <- data.frame(AGE = L[1], SYSBP = L[2], DIABETES = L[3], CURSMOKE = L[4])
  } else {
    l_df <- data.frame(AGE = L[1], SYSBP = L[2], DIABETES = L[3])
  }
  #l_df <- data.frame(AGE = L[1], SYSBP = L[2], DIABETES = L[3], CURSMOKE = L[4])
  survprobs <- survfit(model, newdata = l_df, data = dat)

  s <- summary(survprobs, times = 3652)
  
  return(c(
    sprintf("%.2f%%", 100*(1-s$surv[1])),
    sprintf("%.2f%%", 100*(1-s$upper[1])),
    sprintf("%.2f%%", 100*(1-s$lower[1]))
  ))
  # [1] "0.423" "0.281" "0.637"
  
}


confint_coxph_hazards <- function(model, L, dat) {
  #returns exponentiated Wald CI's for estimate for different risk profiles
  # provides risk profile relative to nonsmoker, no diabetes, median blood pressure
  if (length(L) == 4) {
    L_diff <- L - t(c(median(dat$AGE), median(dat$SYSBP), 0, 0))
    
  }
  if (length(L) == 3) {
    L_diff <- L - t(c(median(dat$AGE), median(dat$SYSBP), 0))
    
  }
  est  <- L_diff %*% model$coefficients
  lower <- est - 1.96*sqrt(L_diff %*% model$var %*% t(L_diff))
  upper <- est + 1.96*sqrt(L_diff %*% model$var %*% t(L_diff) )
  print(exp(c(est, lower, upper)))
  return(exp(c(est, lower, upper)))
  
  
}





# Male subgroups
male_matrix <- matrix(c(40,160, 0, 0,
                        40, quantile(framdat_m$SYSBP, .5), 1, 0,
                        40, quantile(framdat_m$SYSBP, .5), 0, 1,
                        40, 160, 1, 0,
                        40, 160, 1, 1,
                        50, 160, 0, 0,
                        50, quantile(framdat_m$SYSBP, .5), 1, 0,
                        50, quantile(framdat_m$SYSBP, .5), 0, 1,
                        50, 160, 1, 0,
                        50, 160, 1, 1,
                        60, 160, 0, 0,
                        60, 160, 1, 0,
                        60, quantile(framdat_m$SYSBP, .5), 0, 1,
                        60, 160, 1, 0,
                        60, 160, 1, 1,
                        median(framdat_m$AGE),quantile(framdat_m$SYSBP, .5),0 ,0
                        
                        ), byrow = TRUE, ncol = 4)

female_matrix <- matrix(c(40, 160, 0, 0,
                          40, quantile(framdat_f$SYSBP, .5), 1, 0,
                          40, quantile(framdat_f$SYSBP, .5), 0, 1,
                          40, 160, 1, 0,
                          40, 160, 1, 1,
                          50, 160, 0, 0,
                          50, quantile(framdat_f$SYSBP, .5), 1, 0,
                          50, quantile(framdat_f$SYSBP, .5), 0, 1,
                          50, 160, 1, 0,
                          50, 160, 1, 1,
                          60, 160, 0, 0,
                          60, quantile(framdat_f$SYSBP, .5), 1, 0,
                          60, quantile(framdat_f$SYSBP, .5), 0, 1,
                          60, 160, 1, 0,
                          60, 160, 1, 1,
                          median(framdat_f$AGE),quantile(framdat_f$SYSBP, .5),0 ,0
                          
), byrow = TRUE, ncol = 4)


female_matrix_no_sm <- matrix(c(40, 160, 0, 
                          40, quantile(framdat_f$SYSBP, .5), 1, 
                          40, quantile(framdat_f$SYSBP, .5), 0, 
                          40, 160, 1, 
                          40, 160, 1, 
                          50, 160, 0, 
                          50, quantile(framdat_f$SYSBP, .5), 1, 
                          50, quantile(framdat_f$SYSBP, .5), 0, 
                          50, 160, 1, 
                          50, 160, 1, 
                          60, 160, 0, 
                          60, quantile(framdat_f$SYSBP, .5), 1, 
                          60, quantile(framdat_f$SYSBP, .5), 0, 
                          60, 160, 1, 
                          60, 160, 1,
                          median(framdat_f$AGE),quantile(framdat_f$SYSBP, .5),0 
                          
), byrow = TRUE, ncol = 3)

df_cox_cis <- data.frame(`Males Prob. of Stroke @ 10yrs (95% CI)`=character(16), `Females Prob. (95% CI)`=character(16))

for (i in 1:16) {
  
  ci <- confint_coxph(L = male_matrix[i,], model = final_model_m, dat = framdat_m)
  df_cox_cis$`Males Prob. of Stroke @ 10yrs (95% CI)`[i] <- sprintf("%s (%s, %s)", ci[1], ci[2], ci[3])
  
  ci <- confint_coxph(L = female_matrix_no_sm[i,], model = final_model_f,dat = framdat_f)
  df_cox_cis$`Females Prob. (95% CI)`[i] <- sprintf("%s (%s, %s)", ci[1], ci[2], ci[3])
  if (i %in% c(3,5,8,10,13,15))
  {
    df_cox_cis$`Females Prob. (95% CI)`[i] <- "--"
    
  }
  
}
df_cox_cis$`Risk Group` <-
  c("High BP",
    "Diabetes",
    "Smoker",
    "High BP, Diabetes",
    "High BP, Diabetes, Smoker",
    "High BP",
    "Diabetes",
    "Smoker",
    "High BP, Diabetes",
    "High BP, Diabetes, Smoker",
    "High BP",
    "Diabetes",
    "Smoker",
    "High BP, Diabetes",
    "High BP, Diabetes, Smoker",
    "No Risk Factors")

col_order <- c("Risk Group", "Males Prob. of Stroke @ 10yrs (95% CI)", "Females Prob. (95% CI)")
my_data2 <- df_cox_cis[, col_order]
my_data2

gt_cox <- gt(my_data2)

gt_cox <-
  gt_cox |>
  tab_row_group(
    label = "Age = 40",
    rows = 1:5
  ) |>
  tab_row_group(
    label = "Age = 50",
    rows = 6:10
  ) |>
  tab_row_group(
    label = "Age = 60",
    rows = 11:15
    
  ) |>
  tab_row_group(
    label = "Median Age",
    rows = 16
  )  |> tab_caption(caption = "Probability of risk in a 10 year period for men and women in a sample from the Framingham Heart Study, using median age and systolic blood pressure, no diabetes, and for men, no smoking behavior, as the reference levels.")


gt_cox






#######################
# CHANGE OVER TIME    #
#######################
framdat_change <- framdat_raw %>%
  filter(TIMESTRK != 0 ) %>%
  filter(complete.cases(STROKE, AGE, SYSBP, DIABETES, 
                        BPMEDS, PREVCHD, CURSMOKE, TOTCHOL, BMI)) %>%
    filter(RANDID %in% framdat$RANDID) %>%
  group_by(RANDID) %>% filter(max(PERIOD)==3 & min(PERIOD)==1) %>% ungroup




changes <- data.frame(sex = framdat_change[framdat_change$PERIOD==1,]$SEX, 
                      delt_SYSBP = framdat_change[framdat_change$PERIOD==3,]$SYSBP - framdat_change[framdat_change$PERIOD==1,]$SYSBP,
                      delt_BMI = framdat_change[framdat_change$PERIOD==3,]$BMI - framdat_change[framdat_change$PERIOD==1,]$BMI,
                      delt_TOTCHOL = framdat_change[framdat_change$PERIOD==3,]$TOTCHOL - framdat_change[framdat_change$PERIOD==1,]$TOTCHOL,
                      delt_BPMEDS = framdat_change[framdat_change$PERIOD==3,]$BPMEDS - framdat_change[framdat_change$PERIOD==1,]$BPMEDS,
                      delt_PREVCHD = framdat_change[framdat_change$PERIOD==3,]$PREVCHD - framdat_change[framdat_change$PERIOD==1,]$PREVCHD,
                      delt_DIABETES = framdat_change[framdat_change$PERIOD==3,]$DIABETES - framdat_change[framdat_change$PERIOD==1,]$DIABETES,
                      delt_CURSMOKE = framdat_change[framdat_change$PERIOD==3,]$CURSMOKE - framdat_change[framdat_change$PERIOD==1,]$CURSMOKE,
                      delt_AGE = framdat_change[framdat_change$PERIOD==3,]$AGE - framdat_change[framdat_change$PERIOD==1,]$AGE)


############################
# CHANGES OVER STUDY PERIOD#
############################


table_change_m <- changes %>% filter(sex==1) |>
  tbl_summary(
    include = c(delt_SYSBP,delt_BMI,delt_TOTCHOL,delt_BPMEDS,delt_PREVCHD,delt_DIABETES,
                delt_CURSMOKE
    ), # Specify variables to include (optional)
    label = list(
      delt_SYSBP~"Change in SYSBP",
      delt_BMI~"Change in BMI",
      delt_TOTCHOL~"Change in total cholesterol",
      delt_BPMEDS~"Change in BP medication use",
      delt_PREVCHD~"CHD diagnosed after period 1",
      delt_DIABETES~"Diabetes diagnosed after period 1",
      delt_CURSMOKE~"Smoking Status Change"
    ),
    missing = "ifany",                      # Display missing values if any exist
    statistic = list(all_continuous() ~ "{mean} ({sd})")
  )

table_change_f <- changes %>% filter(sex==2) |>
  tbl_summary(
    include = c(delt_SYSBP,delt_BMI,delt_TOTCHOL,delt_BPMEDS,delt_PREVCHD,delt_DIABETES,
                delt_CURSMOKE
    ), # Specify variables to include (optional)
    label = list(
      delt_SYSBP~"Change in SYSBP",
      delt_BMI~"Change in BMI",
      delt_TOTCHOL~"Change in total cholesterol",
      delt_BPMEDS~"Change in BP medication use",
      delt_PREVCHD~"CHD diagnosed after period 1",
      delt_DIABETES~"Diabetes diagnosed after period 1",
      delt_CURSMOKE~"Smoking Status Change"
    ),
    missing = "ifany",                      # Display missing values if any exist
    statistic = list(all_continuous() ~ "{mean} ({sd})")
  )


tbl_merge(
  tbls = list(table_change_m, table_change_f),
  tab_spanner = c("**Male**", "**Female**")
)
###################
#SCHOENFELD PLOTS #
###################


ph_res_m<-cox.zph(final_model_m)

schoen_m <- ggcoxzph(ph_res_m, se=T)

ph_res_f<-cox.zph(final_model_f)

schoen_f <- ggcoxzph(ph_res_f, se=T)






