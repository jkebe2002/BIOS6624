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

confint_coxph <- function(model, L) {
  #returns exponentiated Wald CI's for estimate for different risk profiles
  est  <- t(L) %*% model$coefficients
  lower <- est - 1.96*sqrt(t(L) %*% model$var %*% L)
  upper <- est + 1.96*sqrt(t(L) %*% model$var %*% L)  
  return(exp(c(est, lower, upper)))
}

print(confint_coxph(final_model_m, c(40, median(framdat_m$SYSBP), 0, 0)))

# Define baseline profile
baseline <- c(
  AGE      = median(framdat_m$AGE),
  SYSBP    = median(framdat_m$SYSBP),
  DIABETES = 0,
  CURSMOKE = 0
)

# Define risk profiles
risk_profiles <- expand.grid(
  AGE      = c(40, 50, 60),
  SYSBP    = c(median(framdat_m$SYSBP), quantile(framdat_m$SYSBP, 0.75)),
  DIABETES = c(0, 1),
  CURSMOKE = c(0, 1)
)

# Add labels
risk_profiles <- risk_profiles %>%
  mutate(
    SBP_label      = ifelse(SYSBP == median(framdat_m$SYSBP), "Median SBP", "75th Percentile SBP"),
    DIABETES_label = ifelse(DIABETES == 1, "Diabetes", "No Diabetes"),
    SMOKE_label    = ifelse(CURSMOKE == 1, "Smoker", "Non-Smoker"),
    profile_label  = paste(AGE, "yrs |", SBP_label, "|", DIABETES_label, "|", SMOKE_label)
  )

# Compute CI for each profile relative to baseline
results <- do.call(rbind, lapply(1:nrow(risk_profiles), function(i) {
  L <- as.numeric(risk_profiles[i, c("AGE", "SYSBP", "DIABETES", "CURSMOKE")]) - baseline
  ci <- confint_coxph(final_model_m, L)
  data.frame(
    profile = risk_profiles$profile_label[i],
    est     = ci[1],
    lower   = ci[2],
    upper   = ci[3]
  )
}))

print(results)

# Forest plot
ggplot(results, aes(x = est, y = profile)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "Hazard Ratio (vs. Median Age & SBP, Non-Smoker, No Diabetes)", 
       y = NULL, 
       title = "Risk Profiles - Male Cox Model") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

library(gt)

results_table <- results %>%
  mutate(
    `Hazard Ratio` = round(est, 3),
    `95% CI`       = paste0("(", round(lower, 3), ", ", round(upper, 3), ")")
  ) %>%
  select(Profile = profile, `Hazard Ratio`, `95% CI`)

results_table %>%
  gt() %>%
  tab_header(
    title    = "Risk Profiles - Male Cox Model",
    subtitle = "Hazard Ratios relative to median age & SBP, non-smoker, no diabetes"
  ) %>%
  cols_align(align = "left",   columns = Profile) %>%
  cols_align(align = "center", columns = c(`Hazard Ratio`, `95% CI`)) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  gtsave("risk_profiles_table.html")

browseURL("risk_profiles_table.html")
