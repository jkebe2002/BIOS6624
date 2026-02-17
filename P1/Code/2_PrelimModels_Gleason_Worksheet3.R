###########################################################
## Program Name: 2_PrelimModels_Gleason_Worksheet3.R     ##
## Purpose: Investigate effect of Gleason on PSA.        ##
## Created by: Camille Moore                             ##
###########################################################

## WE WILL USE STAN FOR BAYESIAN ANALYSIS WITH MCMC
## CmdStan is a stand alone program and cmdstanr is an R interface

###MUST INSTALL cmdstanr: https://mc-stan.org/cmdstanr/ for full instructions
###FOR WINDOWS USERS: Rtools must be installed (you may already have this)	
###https://cran.r-project.org/bin/windows/Rtools/
###FOR MAC USERS: Xcode Command Line Tools must be installed (you may have this)
### if not, from the terminal run: xcode-select –install 

# STAN recommends running this in a fresh R session or restarting your current session
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

# Run this function to install the CmdStan software that the R package relies on
library(cmdstanr)
#install_cmdstan()

# This website has examples and vignettes that might be useful:
# https://mc-stan.org/cmdstanr/articles/cmdstanr.html

# Dependencies 
library(cmdstanr)
library(bayesplot)  # diagnostic plots of the MCMC chains
library(posterior)  # for summarizing posterior draws
library(bayestestR) # for calculating higest density posterior intervals
library(mcmcse)     # for calculating MCMCSE's
library(loo)        # for getting model fit statistics (WAIC and LOO-IC)
library(dplyr)
library(tibble)

# Working directory and data
#setwd("/Users/mooreca/Documents/BIOS6624/PSAExample/")
psaclean2 <- read.csv("~/Downloads/prostate.csv")


###########################################################################
#################LOG PSA AND GLEASON SCORE ANALYSIS########################
###########################################################################
####Run a linear model just to orient us using standard methods
psaclean2$gleason<-as.factor(psaclean2$gleason)
psaclean2$lpsa <- log(psaclean2$PSA)
gleasonmodel <-lm(lpsa ~ gleason, data=psaclean2)

sink("Output/WS3_GleasonLinearModelOutput.txt")
print(gleasonmodel)
summary(gleasonmodel)
confint(gleasonmodel)
anova(gleasonmodel)
sink()

####Bayesian Analysis
###########################################################################
# STEP 1: Define the model
# This is a general linear regression set up for STAN
# It is using a half normal prior on sigma
# and normal priors on the regression coefficients
# once the stan file is written it can be reused
###########################################################################

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
  }", dir="Code/STAN", basename='linear_regression_half_normal'
)

###########################################################################
# STEP 2: Compile the Stan program
###########################################################################
mod <- cmdstan_model('Code/STAN/linear_regression_half_normal.stan')


###########################################################################
# STEP 3: Prepare your data to pass to Stan
###########################################################################

# Outcome data
y <- psaclean2$lpsa

# Design matrix in our linear regression
X <- model.matrix(~ grade7 + grade8, data = psaclean2) 

# These are variables needed in Stan
# N is the number of observations in your data set 
# P is the number of columns in the design matrix
N <- nrow(X)
P <- ncol(X)

# Hyperparameters for prior distributions. These change depending on your assumptions
# For this problem the prior for model error SD (sigma) is a half normal and priors for
# the coefficients are normal with a mean and SD, which can be written flexibly 
# with a mean vector (m below) and a SD vector (s below).  

m <- c(1, rep(0, (P-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
s <- rep(100,P) # SD in the prior on regression coefficients --> variance 100^2
sigma_sd <- 100

# create data list to pass to STAN
data_list <- list(
  N = N,
  P = P,
  X = X,
  y = y,
  prior_mean = m,
  prior_sd = s,
  sigma_prior_sd = sigma_sd
)


###########################################################################
# STEP 4: Fit Model
###########################################################################

# 25k iterations is overkill for this problem
# fewer would suffices, especially for preliminary modeling and
# explorartory analysis
fit <- mod$sample(
  data = data_list,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 25000,
  seed = 123
)


###########################################################################
# STEP 5: Summarize the posterior
###########################################################################

# For a quick posterior summary
# If you do not specify variables, this can take a log time 
# since it will summarize the log likelihoods for each observation that we saved
fit$summary(variables=c("beta[1]", "beta[2]", "sigma"))


# Fancier Table of Results of Interest from the Posterior
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

# check that mcmcse < 6% of posterior SD
(100*summary_table$MCSE/summary_table$Std_Dev) < 6

###########################################################################
# STEP 6: Model fit statistics (you only need 1 of these, not all 3)
###########################################################################
# Extract the log likelihood info
loglik_mat <- as_draws_matrix(fit$draws("log_lik"))  # iterations x N

# Get LOO-CV and WAIC
loo_res <- loo(loglik_mat)
print(loo_res)

waic_res <- waic(loglik_mat)
print(waic_res)

# DIC must be computed by hand
# Deviance per iteration: -2 * log-likelihood summed over observations
D_theta <- -2 * rowSums(loglik_mat)
mean_D <- mean(D_theta)

# Deviance at posterior mean
beta_mean <- colMeans(draws_mat[, grep("^beta\\[", colnames(draws_mat))])
sigma_mean <- mean(draws_mat[, "sigma"])
D_bar_theta <- -2 * sum(dnorm(y, mean = X %*% beta_mean, sd = sigma_mean, log = TRUE))

DIC <- 2 * mean_D - D_bar_theta

print(DIC)

###########################################################################
# STEP 6: MCMC Diagnostics
###########################################################################

draws <- as_draws_array(fit$draws())  # dimensions: iterations x chains x parameters
draws_df <- as_draws_df(fit$draws())  # tidy format for bayesplot / ggplot
params <- c("beta[1]", "beta[2]", "beta[3]", "sigma")

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

## Save output if desired ###
sink("Output/WS3_GleasonMCMCOutput.txt")
print(summary_table)
print(DIC)
print(waic_res)
print(loo_res)
print(fit$cmdstan_diagnose())
print(fit$summary(variables=params))
sink()

pdf('Output/WS3_MCMC_DiagnosticPlots.pdf')
mcmc_trace(draws, pars = params)
mcmc_dens_overlay(draws, pars = params)
mcmc_acf(draws, pars = params)
dev.off()


