framdat <- read_csv("~/Downloads/frmgham2.csv")
library(gtsummary)

#Subjects with TIMESTRK 
framdat <- framdat %>%
  filter(TIMESTRK != 0 )

framdat <- framdat %>%
  filter(PERIOD==1 )

framdat <- framdat %>%
  mutate(time_obs = pmin(TIMESTRK, ceiling(10*365.25))) %>%
  select(-c(TIME,HDLC,LDLC,PERIOD,PREVSTRK))
#framdat <- framdat %>%
#  filter(TIME < (365.25*10))

length(unique(framdat$RANDID))
#censored: stroke == 0 and last observation for subject.

library(naniar)



framdat[framdat== "."] <- NA

vis_miss(framdat)


framdat_m <- framdat %>%
  filter(SEX==1)
framdat_f <- framdat %>%
  filter(SEX==2)
# YOUR STATUS / EVENT VARIABLE IS STROKE

table1_m <- framdat_m |>
  tbl_summary(
    include = c(TOTCHOL, AGE, SYSBP, DIABP, CURSMOKE, CIGPDAY, BMI, DIABETES, BPMEDS, HEARTRTE,
                GLUCOSE, educ, PREVCHD, PREVAP, PREVMI,  PREVHYP, STROKE, HYPERTEN, time_obs 
    ), # Specify variables to include (optional)
    missing = "ifany"                      # Display missing values if any exist
  )

table1_f <- framdat_f |>
  tbl_summary(
    include = c(TOTCHOL, AGE, SYSBP, DIABP, CURSMOKE, CIGPDAY, BMI, DIABETES, BPMEDS, HEARTRTE,
                GLUCOSE, educ, PREVCHD, PREVAP, PREVMI,  PREVHYP, STROKE, HYPERTEN, time_obs 
    ), # Specify variables to include (optional)
    missing = "ifany"                      # Display missing values if any exist
  )

tbl_merge(
  tbls = list(table1_m, table1_f),
  tab_spanner = c("**Male**", "**Female**")
  
  
  
  #Correlation matrices
 #Age, diabetes, blood pressure (systolic).
#  Several whether or not to include
 # coronary heart disease, blood pressure meds (Y/N), smoking status
#, Total cholesterol, BMI.
#  - baseline values for all.
  
  cor_fram_m <- framdat_m %>%
    select(where(~ sd(., na.rm = TRUE) > 0)) %>%
    select(AGE, DIABETES, SYSBP, PREVCHD, BPMEDS, CURSMOKE,
           TOTCHOL, BMI) %>%
    cor(use = "complete.obs")  
  
  cor_fram_f <- framdat_f %>%
    select(where(~ sd(., na.rm = TRUE) > 0)) %>%
    select(AGE, DIABETES, SYSBP, PREVCHD, BPMEDS, CURSMOKE,
           TOTCHOL, BMI) %>%
    cor(use = "complete.obs")  
  
  
#Systolic BP, AGE correlated for women
#TOTC
  
  
  library(ggplot2)
  library(reshape2)
  
  # Melt the correlation matrix into long format
  cor_melted_m <- melt(cor_fram_m, varnames = c("Var1", "Var2"), value.name = "correlation")
  
  ggplot(cor_melted_m, aes(x = Var1, y = Var2, fill = correlation)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#d6604d",
      midpoint = 0, limits = c(-1, 1), name = "Correlation"
    ) +
    geom_text(aes(label = round(correlation, 2)), size = 2.5, color = "black") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    coord_fixed()
)

cor_melted_f <- melt(cor_fram_f, varnames = c("Var1", "Var2"), value.name = "correlation")

ggplot(cor_melted_f, aes(x = Var1, y = Var2, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#d6604d",
    midpoint = 0, limits = c(-1, 1), name = "Correlation"
  ) +
  geom_text(aes(label = round(correlation, 2)), size = 2.5, color = "black") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  coord_fixed()






