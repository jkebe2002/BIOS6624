#Model Selection Code

#data is stored in a local drive on my computer

#use hivdat_0 from EDA file 

hivdat <- read.csv("~/Downloads/hiv_6624_final.csv")


hivdat <- hivdat%>%
  mutate(CESD = na_if(CESD, -1))    %>%
  mutate(LIV34 = na_if(LIV34, 9))   %>%
  mutate(FRP = na_if(FRP, 9))       %>%
  mutate(FP = na_if(FP, 9))         %>%
  mutate(HBP = na_if(HBP, 9))       %>%
  mutate(DYSLIP = na_if(DYSLIP, 9)) %>%
  mutate(DIAB = na_if(DIAB, 9))     %>%
  mutate(BMI = na_if(BMI, 999))     %>%
  mutate(income = na_if(income, 9)) %>%
  mutate(HEROPIATE = na_if(HEROPIATE, -9)) 


diff_LEU3N <- hivdat %>%
  filter(years %in% c(0, 2)) %>%
  select(newid, years, LEU3N) %>%
  tidyr::pivot_wider(names_from = years, values_from = LEU3N) %>%
  mutate(diff_LEU3N = `2` - `0`)

diff_VLOAD <- hivdat %>%
  filter(years %in% c(0, 2)) %>%
  select(newid, years, VLOAD) %>%
  tidyr::pivot_wider(names_from = years, values_from = VLOAD) %>%
  mutate(diff_VLOAD = `2` - `0`)

diff_AGG_MENT <- hivdat %>%
  filter(years %in% c(0, 2)) %>%
  select(newid, years, AGG_MENT) %>%
  tidyr::pivot_wider(names_from = years, values_from = AGG_MENT) %>%
  mutate(diff_AGG_MENT = `2` - `0`)

diff_AGG_PHYS <- hivdat %>%
  filter(years %in% c(0, 2)) %>%
  select(newid, years, AGG_PHYS) %>%
  tidyr::pivot_wider(names_from = years, values_from = AGG_PHYS) %>%
  mutate(diff_AGG_PHYS = `2` - `0`)

hivdat_0 <- hivdat[hivdat$years==0,]
hivdat_0$diff_LEU3N<-diff_LEU3N$diff_LEU3N
hivdat_0$diff_VLOAD<-diff_VLOAD$diff_VLOAD
hivdat_0$diff_AGG_MENT<-diff_AGG_MENT$diff_AGG_MENT
hivdat_0$diff_AGG_PHYS<-diff_AGG_PHYS$diff_AGG_PHYS
hivdat_0<- subset(hivdat_0, select=-c(DIAB,LIV34,FP,TCHOL,TRIG,LDL,DYSLIP))
hivdat_0 <- hivdat_0[(!is.na(hivdat_0$diff_LEU3N)&
                       !is.na(hivdat_0$diff_VLOAD)&!is.na(hivdat_0$diff_AGG_MENT)&!is.na(hivdat_0$diff_AGG_PHYS)),]


vars <- c("HASHV","HASHF","income","BMI","HBP","KID","FRP",
          "CESD","SMOKE","DKGRP","HEROPIATE","IDU",
          "ADH","RACE","EDUCBAS","age","ART","hard_drugs", "diff_LEU3N", "diff_VLOAD", "diff_AGG_MENT", "diff_AGG_PHYS")

hivdat_0 <- hivdat_0[complete.cases(hivdat_0[,vars]), ]
#convert variables to factors


#Stan file for lasso regression
stan_file <- write_stan_file("data {
  int<lower=0> N;                  // number of observations
  int<lower=0> P;                  // number of predictors including intercept
  matrix[N, P] X;                  // design matrix (first column = intercept)
  vector[N] y;                     // outcome
  real<lower=0> tau;               // Lasso param.

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
  beta ~ double_exponential(0, tau);

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
  }", dir="Code/", basename='lasso_stan'
)



library(cmdstanr)
library(bayesplot)

# 1. Prepare Data (Standardization is crucial)



X_hiv <- model.matrix(~HASHV+HASHF+income+BMI+HBP+KID+FRP+CESD+SMOKE+DKGRP+HEROPIATE+IDU+ADH+
                          RACE+EDUCBAS+age+ART+hard_drugs, data = hivdat_0)

N <- nrow(X_hiv)
P <- ncol(X_hiv)

#high correlation of 
#FRP, (FP)
#TCHOL , (LDL)
#(IDU), hard_drugs
#ART, (everART)
#paranthesed vars excluded from analusis


m <- c(1, rep(0, (P-1))) # mvnorm mean  (mean in the prior on the regression coefficients)
s <- rep(100,P) # SD in the prior on regression coefficients --> variance 100^2
sigma_sd <- 100



y1 <- hivdat_0$diff_LEU3N
y2 <- hivdat_0$diff_VLOAD
y3 <- hivdat_0$diff_AGG_MENT
y4 <- hivdat_0$diff_AGG_PHYS

data_list1 <- list(N = N, P = P, X = X_hiv, y = y1, tau = .00001,
                   prior_mean = m,
                   prior_sd = s,
                   sigma_prior_sd = sigma_sd)
data_list2 <- list(N = N, P = P, X = X_hiv, y = y2, tau = .00001,
                   prior_mean = m,
                   prior_sd = s,
                   sigma_prior_sd = sigma_sd)
data_list3 <- list(N = N, P = P, X = X_hiv, y = y3, tau = .00001,
                   prior_mean = m,
                   prior_sd = s,
                   sigma_prior_sd = sigma_sd)
data_list4 <- list(N = N, P = P, X = X_hiv, y = y4, tau = .00001,
                   prior_mean = m,
                   prior_sd = s,
                   sigma_prior_sd = sigma_sd)


# 2. Compile Model
mod <- cmdstan_model("Code/lasso_stan.stan")

# 3. Fit Model
fit1 <- mod$sample(
  data = data_list1,
  seed = 123,
  chains = 4,
  parallel_chains = 4
)

fit2 <- mod$sample(
  data = data_list2,
  seed = 123,
  chains = 4,
  parallel_chains = 4
)
fit3 <- mod$sample(
  data = data_list3,
  seed = 123,
  chains = 4,
  parallel_chains = 4
)
fit4 <- mod$sample(
  data = data_list4,
  seed = 123,
  chains = 4,
  parallel_chains = 4
)

# 4. Analyze Results
print(fit1, pars = c("beta", "sigma"))
print(fit2, pars = c("beta", "sigma"))
print(fit3, pars = c("beta", "sigma"))
print(fit4, pars = c("beta", "sigma"))
mcmc_areas(fit$draws("beta"))