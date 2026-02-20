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

# N is the number of observations in your data set 
# P is the number of columns in the design matrix
N <- nrow(X)
P <- ncol(X)


m <- c(1, rep(0, (P-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
s <- rep(100,P) # SD in the prior on regression coefficients --> variance 100^2
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


# check that mcmcse < 6% of posterior SD
(100*summary_table$MCSE/summary_table$Std_Dev) < 6


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

loo1 <- bayes_fit_stats(fit_lvload)
loo2 <-bayes_fit_stats(fit_leu3n)
loo3 <-bayes_fit_stats(fit_agg_ment)
loo4 <-bayes_fit_stats(fit_agg_phys)


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
freq_lvload   <- lm(diff_lvload~hard_drugs+educ_cat+race_cat+smoke_cat+BMI+age, data=hivdat_clean)
freq_leu3n    <- lm(diff_LEU3N~hard_drugs+educ_cat+race_cat+smoke_cat+BMI+age, data=hivdat_clean)
freq_agg_ment <- lm(diff_AGG_MENT~hard_drugs+educ_cat+race_cat+smoke_cat+BMI+age, data=hivdat_clean)
freq_agg_phys <- lm(diff_AGG_PHYS~hard_drugs+educ_cat+race_cat+smoke_cat+BMI+age, data=hivdat_clean)



#Bayesian Regression output table
clinic_vals<- c(.5,50,2,2)
bayes_table_matrix <- matrix(c(b1$Estimate[1],
                               b1$HPDI_2.5[1],
                               b1$HPDI_97.5[1],
                               1-2*pnorm(clinic_vals[1],mean=b1$Estimate[1],sd=b1$Std_Dev[1],lower.tail = TRUE),
                             loo1$looic,
                             b2$Estimate[1],
                             b2$HPDI_2.5[1],
                             b2$HPDI_97.5[1],
                             1-2*pnorm(clinic_vals[2],mean=b2$Estimate[1],sd=b2$Std_Dev[1],lower.tail = TRUE),
                             loo2$looic,
                             b3$Estimate[1],
                             b3$HPDI_2.5[1],
                             b3$HPDI_97.5[1],
                             1-2*pnorm(clinic_vals[3],mean=b3$Estimate[1],sd=b3$Std_Dev[1],lower.tail = TRUE),
                             loo3$looic,
                             b4$Estimate[1],
                             b4$HPDI_2.5[1],
                             b4$HPDI_97.5[1],
                             1-2*pnorm(clinic_vals[4],mean=b4$Estimate[1],sd=b4$Std_Dev[1],lower.tail = TRUE),
                             loo4$looic),
                             byrow = TRUE, ncol = 5)
colnames(bayes_table_matrix) <- c("Hard Drugs Effect Estimate", "2.5% HPDI", "97.5% HPDI", "Post. Prob.", "Model LOO-IC")
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
#t1