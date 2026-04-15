library(tidyverse)

framdat <- read_csv("~/Downloads/frmgham2.csv")
library(gtsummary)


#Subjects with TIMESTRK 
framdat <- framdat %>%
  filter(TIMESTRK != 0 )

framdat <- framdat %>%
  filter(PERIOD==1 )

framdat <- framdat %>%
  mutate(time_obs = pmin(TIMESTRK, ceiling(10*365.25))) %>%
  select(-c(TIME,HDLC,LDLC,PERIOD,PREVSTRK)) %>%
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




# Example using the 'lung' dataset
full_model_m <- coxph(Surv(time_obs,stroke_10yr) ~ AGE+SYSBP+DIABETES+BPMEDS+PREVCHD+CURSMOKE+TOTCHOL+BMI, 
                      data = framdat_m)
null_model_m <- coxph(Surv(time_obs,stroke_10yr) ~ AGE+SYSBP+DIABETES, data = framdat_m)


# Selection based on lowest AIC
final_model_m <- step(
  full_model_m, 
  scope = list(lower = ~ AGE+SYSBP+DIABETES, upper = formula(full_model_m)), 
  direction = "backward"
)



full_model_f <- coxph(Surv(time_obs,stroke_10yr ) ~ AGE+SYSBP+DIABETES+BPMEDS+PREVCHD+CURSMOKE+TOTCHOL+BMI, 
                      data = framdat_f)
null_model_f <- coxph(Surv(time_obs,stroke_10yr ) ~ AGE+SYSBP+DIABETES, data = framdat_f)





final_model_f <- step(
  full_model_f, 
  scope = list(lower = ~ AGE+SYSBP+DIABETES, upper = formula(full_model_f)), 
  direction = "backward"
)


final_model_f_smoking <- coxph(Surv(time_obs,stroke_10yr ) ~ AGE+SYSBP+DIABETES+CURSMOKE, 
                      data = framdat_f)






###################################################
# SURVIVAL CURVES FOR MEN AND WOMEN RISK PROFILES #
###################################################



# res.cox2<-coxph(Surv(DWHFDAYS,DWHF)~TRTMT+
#                   FUNCTCLS+ DIGUSE+DIURET+CHESTX_greater_55+
#                   EJF_PER_above_25,
#                 data=dig_data)
newdata2<-data.frame(TRTMT=c(0,1), FUNCTCLS=0, DIGUSE=0, DIURET=0,
                     CHESTX_greater_55=0, EJF_PER_above_25=0)


#Males, SBP at median and 3rd quartile, diabetes
male_risk_profiles_SBP_diab_40 <- data.frame(AGE = 40, SYSBP = c(128.5,128.5, 141, 141), DIABETES = c(0,1,0,1), CURSMOKE = 0)

male_risk_profiles_SBP_diab_50 <- data.frame(AGE = 50, SYSBP = c(128.5,128.5, 141, 141), DIABETES = c(0,1,0,1), CURSMOKE = 0)

male_risk_profiles_SBP_diab_60 <- data.frame(AGE = 60, SYSBP = c(128.5,128.5, 141, 141), DIABETES = c(0,1,0,1), CURSMOKE = 0)

##############################
# SMOKING RISK PROFILES MALE #
##############################
male_risk_profiles_SBP_diab_40_sm <- data.frame(AGE = 40, SYSBP = c(128.5,128.5, 141, 141), DIABETES = c(0,1,0,1), CURSMOKE = 1)

male_risk_profiles_SBP_diab_50_sm <- data.frame(AGE = 50, SYSBP = c(128.5,128.5, 141, 141), DIABETES = c(0,1,0,1), CURSMOKE = 1)

male_risk_profiles_SBP_diab_60_sm <- data.frame(AGE = 60, SYSBP = c(128.5,128.5, 141, 141), DIABETES = c(0,1,0,1), CURSMOKE = 1)




library(survminer)

xlim_max <- max(survfit(final_model_m)$time)
legend_labs <- c("Median SBP, No Diabetes", "Median SBP, Diabetes", 
                 "3rd Quartile SBP, No Diabetes", "3rd Quartile SBP, Diabetes")

# Store each plot in a list
plots <- list(
  ggsurvplot(survfit(final_model_m, newdata = male_risk_profiles_SBP_diab_40), 
             data = framdat_m, censor = F, conf.int = F, legend = "none",
             xlim = c(0, xlim_max), legend.labs = legend_labs,
             title = "Age 40",
             ggtheme = theme(plot.title = element_text(hjust = 0.5))),
  
  ggsurvplot(survfit(final_model_m, newdata = male_risk_profiles_SBP_diab_50), 
             data = framdat_m, censor = F, conf.int = F, legend = "none",
             xlim = c(0, xlim_max), legend.labs = legend_labs,
             title = "Age 50",
             ggtheme = theme(plot.title = element_text(hjust = 0.5))),
  
  ggsurvplot(survfit(final_model_m, newdata = male_risk_profiles_SBP_diab_60), 
             data = framdat_m, censor = F, conf.int = F, legend = "right",
             xlim = c(0, xlim_max), legend.labs = legend_labs,
             title = "Age 60",
             ggtheme = theme(plot.title = element_text(hjust = 0.5)))
)

arrange_ggsurvplots(plots, nrow = 1, ncol = 3)


plots_sm <- list(
  ggsurvplot(survfit(final_model_m, newdata = male_risk_profiles_SBP_diab_40_sm), 
             data = framdat_m, censor = F, conf.int = F, legend = "none",
             xlim = c(0, xlim_max), legend.labs = legend_labs,
             title = "Age 40",
             ggtheme = theme(plot.title = element_text(hjust = 0.5))),
  
  ggsurvplot(survfit(final_model_m, newdata = male_risk_profiles_SBP_diab_50_sm), 
             data = framdat_m, censor = F, conf.int = F, legend = "none",
             xlim = c(0, xlim_max), legend.labs = legend_labs,
             title = "Age 50",
             ggtheme = theme(plot.title = element_text(hjust = 0.5))),
  
  ggsurvplot(survfit(final_model_m, newdata = male_risk_profiles_SBP_diab_60_sm), 
             data = framdat_m, censor = F, conf.int = F, legend = "right",
             xlim = c(0, xlim_max), legend.labs = legend_labs,
             title = "Age 60",
             ggtheme = theme(plot.title = element_text(hjust = 0.5)))
)
library(cowplot)

# Extract ggplot objects from ggsurvplot
p1 <- plots[[1]]$plot + theme(legend.position = "none")
p2 <- plots[[2]]$plot + theme(legend.position = "none")
p3 <- plots[[3]]$plot + theme(legend.position = "right")

p1_sm <- plots_sm[[1]]$plot + theme(legend.position = "none")
p2_sm <- plots_sm[[2]]$plot + theme(legend.position = "none")
p3_sm <- plots_sm[[3]]$plot + theme(legend.position = "right")
# Extract legend separately
legend <- get_legend(p3)
p3_no_legend <- p3 + theme(legend.position = "none")

legend_sm <- get_legend(p3_sm)
p3_no_legend_sm <- p3_sm + theme(legend.position = "none")


# Arrange: 3 equal plots + legend as 4th element
plot_final_m_nonsmokers <- plot_grid(
  p1, p2, p3_no_legend, legend,
  nrow = 1,
  rel_widths = c(1, 1, 1, 0.5)  # adjust 0.5 to make legend wider/narrower
)
print(plot_final_m_nonsmokers)

plot_final_m_sm <- plot_grid(
  p1_sm, p2_sm, p3_no_legend_sm, legend_sm,
  nrow = 1,
  rel_widths = c(1, 1, 1, 0.5)  # adjust 0.5 to make legend wider/narrower
)
print(plot_final_m_sm)





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
  l_df <- data.frame(AGE = L[1], SYSBP = L[2], DIABETES = L[3], CURSMOKE = L[4])
  survprobs <- survfit(model, newdata = l_df, data = dat)

  s <- summary(survprobs, times = 3652)
  
  return(c(
    sprintf("%.3f%", 100*(1-s$surv[1])),
    sprintf("%.3f", 100*(1-s$upper[1])),
    sprintf("%.3f", 11-s$lower[1])
  ))
  # [1] "0.423" "0.281" "0.637"
  
}


# Male subgroups
male_matrix <- matrix(c(40, quantile(framdat_m$SYSBP, .75), 0, 0,
                        40, quantile(framdat_m$SYSBP, .5), 1, 0,
                        40, quantile(framdat_m$SYSBP, .5), 0, 1,
                        40, quantile(framdat_m$SYSBP, .75), 1, 0,
                        40, quantile(framdat_m$SYSBP, .75), 1, 1,
                        50, quantile(framdat_m$SYSBP, .75), 0, 0,
                        50, quantile(framdat_m$SYSBP, .5), 1, 0,
                        50, quantile(framdat_m$SYSBP, .5), 0, 1,
                        50, quantile(framdat_m$SYSBP, .75), 1, 0,
                        50, quantile(framdat_m$SYSBP, .75), 1, 1,
                        60, quantile(framdat_m$SYSBP, .75), 0, 0,
                        60, quantile(framdat_m$SYSBP, .5), 1, 0,
                        60, quantile(framdat_m$SYSBP, .5), 0, 1,
                        60, quantile(framdat_m$SYSBP, .75), 1, 0,
                        60, quantile(framdat_m$SYSBP, .75), 1, 1
                        ), byrow = TRUE, ncol = 4)

female_matrix <- matrix(c(40, quantile(framdat_f$SYSBP, .75), 0, 0,
                          40, quantile(framdat_f$SYSBP, .5), 1, 0,
                          40, quantile(framdat_f$SYSBP, .5), 0, 1,
                          40, quantile(framdat_f$SYSBP, .75), 1, 0,
                          40, quantile(framdat_f$SYSBP, .75), 1, 1,
                          50, quantile(framdat_f$SYSBP, .75), 0, 0,
                          50, quantile(framdat_f$SYSBP, .5), 1, 0,
                          50, quantile(framdat_f$SYSBP, .5), 0, 1,
                          50, quantile(framdat_f$SYSBP, .75), 1, 0,
                          50, quantile(framdat_f$SYSBP, .75), 1, 1,
                          60, quantile(framdat_f$SYSBP, .75), 0, 0,
                          60, quantile(framdat_f$SYSBP, .5), 1, 0,
                          60, quantile(framdat_f$SYSBP, .5), 0, 1,
                          60, quantile(framdat_f$SYSBP, .75), 1, 0,
                          60, quantile(framdat_f$SYSBP, .75), 1, 1
                          
), byrow = TRUE, ncol = 4)


df_cox_cis <- data.frame(`Males Prob. of Stroke @ 10yrs (95% CI)`=character(15), `Females Prob. (95% CI)`=character(15))

for (i in 1:15) {
  
  ci <- confint_coxph(L = male_matrix[i,], model = final_model_m, dat = framdat_m)
  df_cox_cis$`Males Prob. of Stroke @ 10yrs (95% CI)`[i] <- sprintf("%s (%s, %s)", ci[1], ci[2], ci[3])
  
  ci <- confint_coxph(L = female_matrix[i,], model = final_model_f_smoking,dat = framdat_f)
  df_cox_cis$`Females Prob. (95% CI)`[i] <- sprintf("%s (%s, %s)", ci[1], ci[2], ci[3])

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
    "High BP, Diabetes, Smoker")

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
  )

gt_cox





###################
#SCHOENFELD PLOTS #
###################


ph_res_m<-cox.zph(final_model_m)

schoen_m <- ggcoxzph(ph_res_m, se=T)

ph_res_f<-cox.zph(final_model_f_smoking)

schoen_f <- ggcoxzph(ph_res_f, se=T)

