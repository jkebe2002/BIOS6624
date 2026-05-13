library(kableExtra)
# Load in the preliminary data:

memdat <- read.csv("~/Downloads/PrelimData.csv")

#find correlations for each variable pair
mem_corr <- cor(memdat)

#correlations for each outcome and predictor pair
corr_cvlt_il6  <-  mem_corr[1,3]
corr_cvlt_mcp1 <-  mem_corr[1,4]
corr_cort_il6  <-  mem_corr[2,3]
corr_cort_mcp1 <-  mem_corr[2,4]

corr_range_vec <- c(corr_cvlt_il6, corr_cvlt_mcp1, corr_cort_il6, corr_cort_mcp1)
library(powertools)

# rounded range of correlations
correlations <- seq(0.25, 0.70, by = 0.05)
alphas <- c(0.025, .05/6)
alpha_names <- c("Power, alpha = .05/2", "Power, alpha = .5/6")

# Build table
results_1samp <- data.frame(
  Correlation = paste0("Correlation = ", correlations),
  check.names = FALSE
)

for (i in seq_along(alphas)) {
  col_vals <- sapply(correlations, function(rho) {
    corr.1samp(N = 175, rhoA = rho, alpha = alphas[i], sides = 2)
  })
  results_1samp[[alpha_names[i]]] <- col_vals
}

#print(results, row.names = FALSE)


# MLR models using preliminary data to obtain R squared est.s
mod_ftest1 <- lm(CVLT_CNG3 ~ IL_6+MCP_1, data = memdat)
#summary(mod_ftest1)

mod_ftest2 <- lm(CORT_CNG3 ~ IL_6+MCP_1, data = memdat)
#summary(mod_ftest2)
#.48, .113



correlations <- seq(0.10, 0.50, by = 0.05)

results_rsq <- data.frame(
  Correlation = paste0("R Squared = ", correlations),
  "Power for Overall F, 4 predictors (alpha = 0.05/2)" = sapply(correlations, function(x) {
    round(mlrF.overall(
      N      = 175,
      p      = 4,
      Rsq    = x,
      fsq    = NULL,
      alpha  = 0.05/2,
      power  = NULL,
      random = FALSE,
      v      = FALSE
    ), 3)
  }),
  "Power for Overall F, 4 predictors (alpha = 0.05/6)" = sapply(correlations, function(x) {
    round(mlrF.overall(
      N      = 175,
      p      = 4,
      Rsq    = x,
      fsq    = NULL,
      alpha  = 0.05/6,
      power  = NULL,
      random = FALSE,
      v      = FALSE
    ), 3)
  }),
  check.names = FALSE
)




rho2_vals <- seq(0.4, 0.70, by = 0.05)

# Each ratio gives clean whole-number splits summing to 175:
# n.ratio = 1/6 -> n1 = 150, n2 = 25  (150 + 25 = 175)
# n.ratio = 2/5 -> n1 = 125, n2 = 50  (125 + 50 = 175)
# n.ratio = 3/4 -> n1 = 100, n2 = 75  (100 + 75 = 175)

n1_vals   <- c(150, 125, 100, 75, 50, 25)
n2_vals   <- c(25,  50,  75, 100, 125, 150)
ratios <- n2_vals/n1_vals
col_names <- paste0(
  
  " (n1=", n1_vals, ", n2=", n2_vals, ")"
)

#  "n.ratio = ", c("1/6", "2/5", "3/4", "4/3", "5/2", "6"),
results_2samp <- data.frame(
  "Correlation Difference" = paste0("", rho2_vals),
  check.names = FALSE
)

for (i in seq_along(ratios)) {
  results_2samp[[col_names[i]]] <- sapply(rho2_vals, function(x) {
    round(corr.2samp(
      n1      = n1_vals[i],
      n.ratio = ratios[i],
      rho1    = 0,
      rho2    = x,
      alpha   = 0.05/2,
      power   = NULL,
      sides   = 2,
      v       = FALSE
    ),3)
  })
}

#print(results, row.names = FALSE)

#kableExtra::kable(results)

results_display <- results_2samp
#bolding table values above .8
for (col in names(results_display)[-1]) {
  results_display[[col]] <- ifelse(
    results_display[[col]] >= 0.8,
    cell_spec(results_display[[col]], format = "latex", bold = TRUE),
    as.character(results_display[[col]])
  )
}




for (i in seq_along(ratios)) {
  results_2samp[[col_names[i]]] <- sapply(rho2_vals, function(x) {
    round(corr.2samp(
      n1      = n1_vals[i],
      n.ratio = ratios[i],
      rho1    = 0,
      rho2    = x,
      alpha   = 0.05/6,
      power   = NULL,
      sides   = 2,
      v       = FALSE
    ),3)
  })
}

#print(results, row.names = FALSE)

#kableExtra::kable(results)

results_display2 <- results_2samp
#bolding table values above .8
for (col in names(results_display2)[-1]) {
  results_display2[[col]] <- ifelse(
    results_display2[[col]] >= 0.8,
    cell_spec(results_display2[[col]], format = "latex", bold = TRUE),
    as.character(results_display2[[col]])
  )
}


# Perform selection using 

#F-tests and p-values (automated)
#library(StepReg)
#AIC
#BIC

#Choose 2 different choices for lambda
#alpha=1 LASSO
#alpha %in% (0,1) is elastic net
#lasso
## LASSO used extensively in high-dimensional problems...
#elastic net



# simulation chunk
#generate 6 datasets
#run function on each dataset
#add metrics to list


# summarize metrics


# Summarize your findings in terms of true positives 
# (what percentage of the time are variables X1 - X5 in the final model) 
# and false positive rates (confusion matrix is another name for this). 
# In other words, how often does your final model include any of the 15 variables 
# that are not associated with the outcome and how often does your final model include each of the 5 variables with a known relationship with Y. 
# For the variables that are retained in the model, what is the bias of the parameter estimates,
# what is the coverage of the 95% CI, and what is the type I and II error of the variables selected.  

# Type 1: 
# how often does the model
# For the variables retained, you can do testing on the final model 
# (in other words refit the final model like no model selection was done) 
# or propose an alternative approach to testing since we believe there is bias in just fitting the final model like there was no variable selection done.



# Questions for Carter:

#Type I and Type II error -  variable specific?
#Should the confusion matrix be variable specific
#What is the difference between doing this and the confusion matrix?
#coverage - variable specific
#Choosing alpha for elastic net... leave it at .5... can choose others 


#This is a report to a clinician.


#Two lasso lambda's: the estimate that minimizes the mse and the lambda 1 standard deviation higher than this one.

# is parameter significant in simulation
# refit lm with the variables selected from each procedure
# find bias of estimates, errors, coverage...
# percentage of coverage for each variable







# Perform selection using 

#F-tests and p-values (automated)
#library(StepReg)
#AIC
#BIC

#Choose 2 different choices for lambda
#alpha=1 LASSO
#alpha %in% (0,1) is elastic net
#lasso
## LASSO used extensively in high-dimensional problems...
#elastic net



# simulation chunk
#generate 6 datasets
#run function on each dataset
#add metrics to list


# summarize metrics


# Summarize your findings in terms of true positives 
# (what percentage of the time are variables X1 - X5 in the final model) 
# and false positive rates (confusion matrix is another name for this). 
# In other words, how often does your final model include any of the 15 variables 
# that are not associated with the outcome and how often does your final model include each of the 5 variables with a known relationship with Y. 
# For the variables that are retained in the model, what is the bias of the parameter estimates,
# what is the coverage of the 95% CI, and what is the type I and II error of the variables selected.  

# Type 1: 
# how often does the model
# For the variables retained, you can do testing on the final model 
# (in other words refit the final model like no model selection was done) 
# or propose an alternative approach to testing since we believe there is bias in just fitting the final model like there was no variable selection done.



# Questions for Carter:

#Type I and Type II error -  variable specific?
#Should the confusion matrix be variable specific
#What is the difference between doing this and the confusion matrix?
#coverage - variable specific
#Choosing alpha for elastic net... leave it at .5... can choose others 


#This is a report to a clinician.


#Two lasso lambda's: the estimate that minimizes the mse and the lambda 1 standard deviation higher than this one.

# is parameter significant in simulation
# refit lm with the variables selected from each procedure
# find bias of estimates, errors, coverage...
# percentage of coverage for each variable



