# Models of outcomes as function of hard drugs etc.
# Note that code from Camille Moore, PhD was implemented heavily in this section!

# Dependencies 
library(cmdstanr)
library(bayesplot)  # diagnostic plots of the MCMC chains
library(posterior)  # for summarizing posterior draws
library(bayestestR) # for calculating higest density posterior intervals
library(mcmcse)     # for calculating MCMCSE's
library(loo)        # for getting model fit statistics (WAIC and LOO-IC)
library(dplyr)
library(tibble)
library(ggplot2)
library(gridExtra)


#############
#CREATE HISTOGRAMS SHOWING DIFFERENCES

#Here, we create histgrams of distribution of each outcome at year 0 and year 2 stratified by drug use.
#############
create_histograms <- function{
  p_lvl0 <- ggplot(data=hivdat_clean, aes(x=lVload_0, group=hard_drugs_0, fill=hard_drugs_0)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_bw()+xlab("Log Viral Load (Baseline)")+
    labs(fill = "Hard Drug Use")
  
  p_lvl2 <- ggplot(data=hivdat_clean, aes(x=lVload_2, group=hard_drugs_0, fill=hard_drugs_0)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_bw()+xlab("Log Viral Load (Year 2)")+
    labs(fill = "Hard Drug Use")
  
  p_lvl2 <- ggplot(data=hivdat_clean, aes(x=lVload_2, group=hard_drugs_0, fill=hard_drugs_0)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_bw()+xlab("Log Viral Load (Year 2)")+
    labs(fill = "Hard Drug Use")
  
  
  p_leu3n0 <- ggplot(data=hivdat_clean, aes(x=LEU3N_0, group=hard_drugs_0, fill=hard_drugs_0)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_bw()+xlab("CD4+ T Cell Count (Baseline)")+
    labs(fill = "Hard Drug Use")
  
  
  p_leu3n2 <- ggplot(data=hivdat_clean, aes(x=LEU3N_2, group=hard_drugs_0, fill=hard_drugs_0)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_bw()+xlab("CD4+ T Cell Count ( Year 2)")+
    labs(fill = "Hard Drug Use")
  
  
  p_aggment0 <- ggplot(data=hivdat_clean, aes(x=AGG_MENT_0, group=hard_drugs_0, fill=hard_drugs_0)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_bw()+xlab("Mental Quality of Life Score (Baseline)")+
    labs(fill = "Hard Drug Use")
  
  
  p_aggment2 <- ggplot(data=hivdat_clean, aes(x=AGG_MENT_2, group=hard_drugs_0, fill=hard_drugs_0)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_bw()+xlab("Mental Quality of Life Score (Year 2)")+
    labs(fill = "Hard Drug Use")
  
  
  
  p_aggphys0 <- ggplot(data=hivdat_clean, aes(x=AGG_PHYS_0, group=hard_drugs_0, fill=hard_drugs_0)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_bw()+xlab("Physical Quality of Life Score (Baseline)")+
    labs(fill = "Hard Drug Use")
  
  p_aggphys2 <- ggplot(data=hivdat_clean, aes(x=AGG_PHYS_2, group=hard_drugs_0, fill=hard_drugs_0)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_bw()+xlab("Physical Quality of Life Score (Baseline)")+
    labs(fill = "Hard Drug Use")
  
  #Put this line in the main report to show histograms
  grid.arrange(p_lvl0, p_leu3n0, p_aggment0, p_aggphys0, p_lvl2, p_leu3n2, p_aggment2, p_aggphys2, nrow=2, ncol=4)
}



############
#STAN FILE##
############
stan_file <- write_stan_file("data {
  int<lower=0> N;                  // number of observations
  int<lower=0> P;                  // number of predictors including intercept
  matrix[N, P] X;                  // design matrix (first column = intercept)
  vector[N] y;                     // outcome

  vector[P] prior_mean;            // prior means for each beta
  vector<lower=0>[P] prior_sd;     // prior SDs for each beta

  real<lower=0> sigma_prior_sd;    // SD for half-normal prior on sigma
}

parameters {
  vector[P] beta;                  // regression coefficients
  real<lower=0> sigma;             // residual SD
}

model {
  // Vectorized priors for regression coefficients
  beta ~ normal(prior_mean, prior_sd);

  // Half-normal prior for sigma
  sigma ~ normal(0, sigma_prior_sd);

  // Likelihood
  y ~ normal(X * beta, sigma);
  

}

generated quantities {
  // log likelihood for each observation for calculating model fit stats
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | X[n] * beta, sigma);}
  }", dir="~/Documents/GitHub/BIOS6624/P1/Code", basename='linear_regression_half_normal'
)

#compile stan prog.
mod <- cmdstan_model('~/Documents/GitHub/BIOS6624/P1/Code/linear_regression_half_normal.stan')

#outcomes:

lvload <- hivdat_clean[!is.na(hivdat_clean$BMI),]$diff_lvload
leu3n <- hivdat_clean[!is.na(hivdat_clean$BMI),]$diff_LEU3N
agg_ment <- hivdat_clean[!is.na(hivdat_clean$BMI),]$diff_AGG_MENT
agg_phys <- hivdat_clean[!is.na(hivdat_clean$BMI),]$diff_AGG_PHYS

# Design matrix in our linear regression
X <- model.matrix(~ hard_drugs_0+educ_cat+race_cat+smoke_cat+BMI+age_0,
                  data = hivdat_clean,
                  na.action=na.pass) 

#Design matrix with no hard drugs to compute change in loo-ic
X_no_hard_drugs <- model.matrix(~educ_cat+race_cat+smoke_cat+BMI+age_0,
                                data = hivdat_clean,
                                na.action=na.pass) 

#Desing matrix with adherence variable to compute effect of adherence on hard_drugs estimate
X_adh <- model.matrix(~ hard_drugs_0+ADH+educ_cat+race_cat+smoke_cat+BMI+age_0,
                      data = hivdat_clean,
                      na.action=na.pass) 

# N is the number of observations in your data set 
# P is the number of columns in the design matrix
N <- nrow(X)
P <- ncol(X)

P_adh <- ncol(X_adh)
P_no_hard_drugs <- ncol(X_no_hard_drugs)



m <- c(1, rep(0, (P-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
m_adh <- c(1, rep(0, (P_adh-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
m_no_hard_drugs <- c(1, rep(0, (P_no_hard_drugs-1))) # mvnorm mean  (mean in the prior on the regression coefficients)

s <- rep(100,P) # SD in the prior on regression coefficients --> variance 100^2
s_adh <- rep(100,P_adh) # SD in the prior on regression coefficients --> variance 100^2
s_no_hard_drugs<- rep(100,P_no_hard_drugs) # SD in the prior on regression coefficients --> variance 100^2

sigma_sd <- 100

# create data lists to pass to STAN for each outcome
data_list_lvload <- list(
  N = N,
  P = P,
  X = X,
  y = lvload,
  prior_mean = m,
  prior_sd = s,
  sigma_prior_sd = sigma_sd
)

data_list_leu3n <- list(
  N = N,
  P = P,
  X = X,
  y = leu3n,
  prior_mean = m,
  prior_sd = s,
  sigma_prior_sd = sigma_sd
)

data_list_agg_ment <- list(
  N = N,
  P = P,
  X = X,
  y = agg_ment,
  prior_mean = m,
  prior_sd = s,
  sigma_prior_sd = sigma_sd
)

data_list_agg_phys <- list(
  N = N,
  P = P,
  X = X,
  y = agg_phys,
  prior_mean = m,
  prior_sd = s,
  sigma_prior_sd = sigma_sd
)

#Model fits for each outcome

#log viral load
fit_lvload <- mod$sample(
  data = data_list_lvload ,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 25000,
  seed = 2319
)

#t cell count
fit_leu3n <- mod$sample(
  data = data_list_leu3n,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 25000,
  seed = 2319
)

#mental health score
fit_agg_ment <- mod$sample(
  data = data_list_agg_ment,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 25000,
  seed = 2319
)

#phys health score
fit_agg_phys <- mod$sample(
  data = data_list_agg_phys,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 25000,
  seed = 2319
)



###################################
######ADHERENCE MODELS#############
###################################

# create data lists to pass to STAN for each outcome
data_list_lvload_adh <- list(
  N = N,
  P = P_adh,
  X = X_adh,
  y = lvload,
  prior_mean = m_adh,
  prior_sd = s_adh,
  sigma_prior_sd = sigma_sd
)

data_list_leu3n_adh <- list(
  N = N,
  P = P_adh,
  X = X_adh,
  y = leu3n,
  prior_mean = m_adh,
  prior_sd = s_adh,
  sigma_prior_sd = sigma_sd
)

data_list_agg_ment_adh <- list(
  N = N,
  P = P_adh,
  X = X_adh,
  y = agg_ment,
  prior_mean = m_adh,
  prior_sd = s_adh,
  sigma_prior_sd = sigma_sd
)

data_list_agg_phys_adh <- list(
  N = N,
  P = P_adh,
  X = X_adh,
  y = agg_phys,
  prior_mean = m_adh,
  prior_sd = s_adh,
  sigma_prior_sd = sigma_sd
)

#Model fits for each outcome

#log viral load
fit_lvload_adh <- mod$sample(
  data = data_list_lvload_adh,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 25000,
  seed = 2319
)

#t cell count
fit_leu3n_adh <- mod$sample(
  data = data_list_leu3n_adh,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 25000,
  seed = 2319
)

#mental health score
fit_agg_ment_adh <- mod$sample(
  data = data_list_agg_ment_adh,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 25000,
  seed = 2319
)

#phys health score
fit_agg_phys_adh <- mod$sample(
  data = data_list_agg_phys_adh,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 25000,
  seed = 2319
)




#######################################
######NO HARD DRUGS MODELS#############
#######################################

# create data lists to pass to STAN for each outcome
data_list_lvload_no_hard_drugs <- list(
  N = N,
  P = P_no_hard_drugs,
  X = X_no_hard_drugs,
  y = lvload,
  prior_mean = m_no_hard_drugs,
  prior_sd = s_no_hard_drugs,
  sigma_prior_sd = sigma_sd
)

data_list_leu3n_no_hard_drugs <- list(
  N = N,
  P = P_no_hard_drugs,
  X = X_no_hard_drugs,
  y = leu3n,
  prior_mean = m_no_hard_drugs,
  prior_sd = s_no_hard_drugs,
  sigma_prior_sd = sigma_sd
)

data_list_agg_ment_no_hard_drugs <- list(
  N = N,
  P = P_no_hard_drugs,
  X = X_no_hard_drugs,
  y = agg_ment,
  prior_mean = m_no_hard_drugs,
  prior_sd = s_no_hard_drugs,
  sigma_prior_sd = sigma_sd
)

data_list_agg_phys_no_hard_drugs <- list(
  N = N,
  P = P_no_hard_drugs,
  X = X_no_hard_drugs,
  y = agg_phys,
  prior_mean = m_no_hard_drugs,
  prior_sd = s_no_hard_drugs,
  sigma_prior_sd = sigma_sd
)

#Model fits for each outcome

#log viral load
fit_lvload_no_hard_drugs <- mod$sample(
  data = data_list_lvload_no_hard_drugs,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 25000,
  seed = 2319
)

#t cell count
fit_leu3n_no_hard_drugs <- mod$sample(
  data = data_list_leu3n_no_hard_drugs,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 25000,
  seed = 2319
)

#mental health score
fit_agg_ment_no_hard_drugs <- mod$sample(
  data = data_list_agg_ment_no_hard_drugs,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 25000,
  seed = 2319
)

#phys health score
fit_agg_phys_no_hard_drugs <- mod$sample(
  data = data_list_agg_phys_no_hard_drugs,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 25000,
  seed = 2319
)

# function to return LOO and WAIC:

bayes_fit_stats <- function(fit){
  # Extract the log likelihood info
  loglik_mat <- as_draws_matrix(fit$draws("log_lik"))  # iterations x N
  
  # Get LOO-CV and WAIC
  loo_res <- loo(loglik_mat)
  print(loo_res)
  
  waic_res <- waic(loglik_mat)
  print(waic_res)
  return(loo_res)
}


####################
#COMPUTE DELTA LOOIC
####################
loo1_nhd <- bayes_fit_stats(fit_lvload_no_hard_drugs)
loo2_nhd <- bayes_fit_stats(fit_leu3n_no_hard_drugs)
loo3_nhd <- bayes_fit_stats(fit_agg_ment_no_hard_drugs)
loo4_nhd <- bayes_fit_stats(fit_agg_phys_no_hard_drugs)

loo1 <- bayes_fit_stats(fit_lvload)
loo2 <-bayes_fit_stats(fit_leu3n)
loo3 <-bayes_fit_stats(fit_agg_ment)
loo4 <-bayes_fit_stats(fit_agg_phys)

loo1_delt <- loo1_nhd$looic-loo1$looic
loo2_delt <- loo2_nhd$looic-loo2$looic
loo3_delt <- loo3_nhd$looic-loo3$looic
loo4_delt <- loo4_nhd$looic-loo4$looic
#######################################
#function to print results summary table 

bayes_table <- function(fit) {
  # Extract posterior draws from cmdstanr fit
  draws <- fit$draws()  
  
  draws_mat <- as_draws_matrix(draws)
  params <- colnames(draws_mat)
  
  # Exclude non-parameter columns from summarization: lp__ and log_lik[n]
  params <- params[!grepl("lp__|log_lik", params)]
  
  
  # Build summary table
  summary_table <- lapply(params, function(p) {
    vals <- as.numeric(draws_mat[, p])
    
    # Monte Carlo standard error
    # Used to determine if we have run enough iterations
    # rule of thumb: MCSE should be less than 6.27% of the posterior standard deviation
    mcse_val <- mcmcse::mcse(vals)$se
    
    # Effective sample size
    ess_val <- ess_bulk(vals)
    
    # 95% HPDI
    hpd <- hdi(vals, ci = 0.95)
    
    tibble(
      Parameter = p,
      Estimate  = mean(vals),
      MCSE      = mcse_val,
      Std_Dev   = sd(vals),
      HPDI_2.5  = hpd$CI_low,
      HPDI_97.5 = hpd$CI_high,
      ESS       = ess_val
    )
  }) %>% bind_rows()
  
  print(summary_table)
  return(summary_table)
}

b1 <- bayes_table(fit_lvload)
b2 <-bayes_table(fit_leu3n)
b3 <- bayes_table(fit_agg_ment)
b4 <-bayes_table(fit_agg_phys)


b1_adh <- bayes_table(fit_lvload_adh)
b2_adh <- bayes_table(fit_leu3n_adh)
b3_adh <- bayes_table(fit_agg_ment_adh)
b4_adh <- bayes_table(fit_agg_phys_adh)

# check that mcmcse < 6% of posterior SD
(100*summary_table$MCSE/summary_table$Std_Dev) < 6







# fit plots and diagnostics function

bayes_diagnostics <- function(fit){
  draws <- as_draws_array(fit$draws())  # dimensions: iterations x chains x parameters
  draws_df <- as_draws_df(fit$draws())  # tidy format for bayesplot / ggplot
  params <- c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "beta[6]", "beta[7]", "sigma")
  
  # Trace plots
  # Trace plots show whether chains are mixing well and exploring the parameter space.
  # Chains should mix well (lines overlap)
  # No long trends (no “stickiness”)
  # Ideally all chains overlap around the same mean
  
  mcmc_trace(draws, pars = params)
  
  
  # R-hat and ESS diagnostics
  # R-hat ≈ 1.00 → converged
  # ESS should be reasonably large (e.g., >100–200 per parameter)
  fit$summary(variables=params)
  
  
  # Posterior Density Plots
  # Overlaid densities from multiple chains should match closely
  # If densities differ substantially between chains → poor mixing
  mcmc_dens_overlay(draws, pars = params)
  
  # Auto-correlation 
  # Autocorrelation should decay quickly
  # High autocorrelation → may need more iterations or thinning
  mcmc_acf(draws, pars = params)
  
  # Diagnostics from cmdstan
  fit$cmdstan_diagnose()
  
}


#frequentist models for each outcome:
freq_lvload   <- lm(diff_lvload~hard_drugs_0+educ_cat+race_cat+smoke_cat+BMI+age_0, data=hivdat_clean)
freq_leu3n    <- lm(diff_LEU3N~hard_drugs_0+educ_cat+race_cat+smoke_cat+BMI+age_0, data=hivdat_clean)
freq_agg_ment <- lm(diff_AGG_MENT~hard_drugs_0+educ_cat+race_cat+smoke_cat+BMI+age_0, data=hivdat_clean)
freq_agg_phys <- lm(diff_AGG_PHYS~hard_drugs_0+educ_cat+race_cat+smoke_cat+BMI+age_0, data=hivdat_clean)



#Bayesian Regression output table
clinic_vals<- c(.5,50,2,2)

postprob <- function(q, mean, sd) {
  return(1+pnorm(-q,mean=mean, sd=sd, lower.tail = FALSE)-pnorm(q,mean=mean,sd=sd,lower.tail = TRUE))
}
bayes_table_matrix <- matrix(c(b1$Estimate[2],
                               b1$HPDI_2.5[2],
                               b1$HPDI_97.5[2],
                               postprob(q=clinic_vals[1],mean=b1$Estimate[2],sd=b1$Std_Dev[2]),
                               loo1$looic,
                               loo1_delt,
                               b2$Estimate[2],
                               b2$HPDI_2.5[2],
                               b2$HPDI_97.5[2],
                               pnorm(clinic_vals[2],mean=b2$Estimate[2],sd=b2$Std_Dev[2],lower.tail = TRUE),
                               loo2$looic,
                               loo2_delt,
                               b3$Estimate[2],
                               b3$HPDI_2.5[2],
                               b3$HPDI_97.5[2],
                               pnorm(clinic_vals[3],mean=b3$Estimate[2],sd=b3$Std_Dev[2],lower.tail = TRUE),
                               loo3$looic,
                               loo3_delt,
                               b4$Estimate[2],
                               b4$HPDI_2.5[2],
                               b4$HPDI_97.5[2],
                               pnorm(clinic_vals[4],mean=b4$Estimate[1],sd=b4$Std_Dev[2],lower.tail = TRUE),
                               loo4$looic,
                               loo4_delt,),
                             byrow = TRUE, ncol = 6)
colnames(bayes_table_matrix) <- c("Hard Drugs Effect Estimate", "2.5% HPDI", "97.5% HPDI", "Post. Prob.", "Model LOO-IC", "Delta LOO-IC (Model w/o Hard Drug Use)")
rownames(bayes_table_matrix) <- c("log Viral Load", "CD4+ T Cell Count", "Mental Quality of Life Score", "Physical Quality of Life Score")
knitr::kable(bayes_table_matrix)
#Table 1:

library(table1)

label(hivdat_clean$lVload) <- "log(HIV copies per mL of blood)"
label(hivdat_clean$LEU3N) <- "CD4+ T cell Count"
label(hivdat_clean$AGG_MENT) <- "Mental Quality of Life Score"
label(hivdat_clean$AGG_PHYS) <- "Physical Quality of Life Score"
label(hivdat_clean$educ_cat) <- "Education Level"
label(hivdat_clean$race_cat) <- "Race"
label(hivdat_clean$smoke_cat) <- "Smoking Status"
label(hivdat_clean$BMI) <- "Body Mass Index (BMI)"
label(hivdat_clean$age) <- "Age"
hivdat_clean$hard_drugs <- as.factor(hivdat_clean$hard_drugs)
label(hivdat_clean$hard_drugs) <- "Hard Drug Usage (Yes=1)"


t1 <- table1(~lVload+LEU3N+AGG_MENT+AGG_PHYS+educ_cat+race_cat+smoke_cat+BMI+age|hard_drugs, hivdat_clean,
             caption="Table of Descriptive Statistics by Drug Use.")

t1<- t1flex(t1,tablefn = "qflextable")
t1 %>% flextable::fontsize(size = 7, part = "all") %>% 
  # Reduce font size
  flextable::autofit()


#t1 is your table 1 object.


##########################
#Freq. Vs. Bayes Table####
##########################
lvload_freq_est <- summary(freq_lvload)$coefficients
leu3n_freq_est <- summary(freq_leu3n)$coefficients
aggment_freq_est <- summary(freq_agg_ment)$coefficients
aggphys_freq_est <- summary(freq_agg_phys)$coefficients
bayes_vs_freq_matrix <- round(matrix(c(b1$Estimate[2],
                                       b1$Std_Dev[2],
                                       b1$HPDI_2.5[2],
                                       b1$HPDI_97.5[2],
                                       lvload_freq_est[2,"Estimate"],
                                       lvload_freq_est[2,"Std. Error"],
                                       lvload_freq_est[2,"Estimate"]-1.96*lvload_freq_est[2,"Std. Error"],
                                       lvload_freq_est[2,"Estimate"]+1.96*lvload_freq_est[2,"Std. Error"],
                                       b2$Estimate[2],
                                       b2$Std_Dev[2],
                                       b2$HPDI_2.5[2],
                                       b2$HPDI_97.5[2],
                                       leu3n_freq_est[2,"Estimate"],
                                       leu3n_freq_est[2,"Std. Error"],
                                       leu3n_freq_est[2,"Estimate"]-1.96*leu3n_freq_est[2,"Std. Error"],
                                       leu3n_freq_est[2,"Estimate"]+1.96*leu3n_freq_est[2,"Std. Error"],
                                       b3$Estimate[2],
                                       b3$Std_Dev[2],
                                       b3$HPDI_2.5[2],
                                       b3$HPDI_97.5[2],
                                       aggment_freq_est[2,"Estimate"],
                                       aggment_freq_est[2,"Std. Error"],
                                       aggment_freq_est[2,"Estimate"]-1.96*aggment_freq_est[2,"Std. Error"],
                                       aggment_freq_est[2,"Estimate"]+1.96*aggment_freq_est[2,"Std. Error"],
                                       b4$Estimate[2],
                                       b4$Std_Dev[2],
                                       b4$HPDI_2.5[2],
                                       b4$HPDI_97.5[2],
                                       aggphys_freq_est[2,"Estimate"],
                                       aggphys_freq_est[2,"Std. Error"],
                                       aggphys_freq_est[2,"Estimate"]-1.96*aggphys_freq_est[2,"Std. Error"],
                                       aggphys_freq_est[2,"Estimate"]+1.96*aggphys_freq_est[2,"Std. Error"]),
                                     byrow = TRUE, ncol = 8),3)


colnames(bayes_vs_freq_matrix) <- c("Hard Drugs Effect Posterior Mean (Bayes)","Bayes Std. Dev.", "2.5% HPDI", "97.5% HPDI", 
                                    "Hard Drugs Effect Estimate (Freq.)", "Freq. Std. Dev.", "2.5% Wald CI", "97.5% Wald CI")
rownames(bayes_vs_freq_matrix) <- c("log Viral Load", "CD4+ T Cell Count", "Mental Quality of Life Score", "Physical Quality of Life Score")


#This is your frequentist comparison table:
kableExtra::kable(bayes_vs_freq_matrix)



#####################
#ADHERENCE TABLE#####
#####################
b1$Estimate[2],
b1$HPDI_2.5[2],
b1$HPDI_97.5[2],
b1_adh$Estimate[2],
b1_adh$HPDI_2.5[2],
b1_adh$HPDI_97.5[2],
b1$Estimate[2]-b1_adh$Estimate[2],
abs(b1$Estimate[2]-b1_adh$Estimate[2])/b1$Estimate[2],

b2$Estimate[2],
b2$HPDI_2.5[2],
b2$HPDI_97.5[2],
b2_adh$Estimate[2],
b2_adh$HPDI_2.5[2],
b2_adh$HPDI_97.5[2],
b2$Estimate[2]-b2_adh$Estimate[2],
abs(b2$Estimate[2]-b2_adh$Estimate[2])/b2$Estimate[2],

b3$Estimate[2],
b3$HPDI_2.5[2],
b3$HPDI_97.5[2],
b3_adh$Estimate[2],
b3_adh$HPDI_2.5[2],
b3_adh$HPDI_97.5[2],
b3$Estimate[2]-b3_adh$Estimate[2],
abs(b1$Estimate[2]-b3_adh$Estimate[2])/b3$Estimate[2],

b4$Estimate[2],
b4$HPDI_2.5[2],
b4$HPDI_97.5[2],
b4_adh$Estimate[2],
b4_adh$HPDI_2.5[2],
b4_adh$HPDI_97.5[2],
b4$Estimate[2]-b4_adh$Estimate[2],
abs(b4$Estimate[2]-b4_adh$Estimate[2])/b4$Estimate[2],





adh_tab <- matrix(c(
  sprintf("%.3f (%.3f, %.3f)", b1$Estimate[2], b1$HPDI_2.5[2], b1$HPDI_97.5[2]),
  sprintf("%.3f (%.3f, %.3f)", b1_adh$Estimate[2], b1_adh$HPDI_2.5[2], b1_adh$HPDI_97.5[2]),
  sprintf("%.3f", b1$Estimate[2]-b1_adh$Estimate[2]),
  sprintf("%.1f%%", 100*abs(b1$Estimate[2]-b1_adh$Estimate[2])/b1$Estimate[2]),
  
  sprintf("%.3f (%.3f, %.3f)", b2$Estimate[2], b2$HPDI_2.5[2], b2$HPDI_97.5[2]),
  sprintf("%.3f (%.3f, %.3f)", b2_adh$Estimate[2], b2_adh$HPDI_2.5[2], b2_adh$HPDI_97.5[2]),
  sprintf("%.3f", b2$Estimate[2]-b2_adh$Estimate[2]),
  sprintf("%.1f%%", 100*abs(b2$Estimate[2]-b2_adh$Estimate[2])/b2$Estimate[2]),
  
  sprintf("%.3f (%.3f, %.3f)", b3$Estimate[2], b3$HPDI_2.5[2], b3$HPDI_97.5[2]),
  sprintf("%.3f (%.3f, %.3f)", b3_adh$Estimate[2], b3_adh$HPDI_2.5[2], b3_adh$HPDI_97.5[2]),
  sprintf("%.3f", b3$Estimate[2]-b3_adh$Estimate[2]),
  sprintf("%.1f%%", 100*abs(b3$Estimate[2]-b3_adh$Estimate[2])/b3$Estimate[2]),
  
  sprintf("%.3f (%.3f, %.3f)", b4$Estimate[2], b4$HPDI_2.5[2], b4$HPDI_97.5[2]),
  sprintf("%.3f (%.3f, %.3f)", b4_adh$Estimate[2], b4_adh$HPDI_2.5[2], b4_adh$HPDI_97.5[2]),
  sprintf("%.3f", b4$Estimate[2]-b4_adh$Estimate[2]),
  sprintf("%.1f%%", 100*abs(b4$Estimate[2]-b4_adh$Estimate[2])/b4$Estimate[2])
), byrow = TRUE, ncol = 4)

colnames(table_mat) <- c("Crude Estimate (95% HPDI)", "Adjusted Estimate (95% HPDI)", "Difference", "Percent Difference")
rownames(table_mat) <- c("Log Viral Load", "CD4+ T cell Count", "Mental Q.O.L. Score", "Physical Q.O.L. Score")

table_mat
###adherence into 2 categories
## Descriptively talk about adherence 
#t1