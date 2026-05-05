library(hdrm)
library(olsrr)


#set that seed yo!
set.seed(23192319)

#dataframe for study results!
sim_results <- data.frame(var = NULL, bias = NULL, coverage = NULL,
                          type1 = NULL, type2 = NULL, trueP = NULL, falseP = NULL,
                          method = NULL, case = NULL)

#List of truly significant variables
true_var <- c("V01", "V02", "V03", "V04", "V05")
#List of true coefficients
true_coef <- seq(.5/3, 2.5/3, .5/3)

# here I write a function to refit an lm() and look at true positive, false positive, type I, II error, bias and coverage.

selection_diag <- function(model, true_var, true_coef, X, y, method, case) {
  # model.    : model object
  # true_var  : list of strings for the names of the variables with association with y
  # true_coef : the true coefficients for the relationship between true_var and y
  # X,y       : simulated data we are regressing
  # method    : which method are we testing selection diagnostics?
  # case      : 1a,1b,2a,2b (varying N, correlation structure...)
  
  trueP = 0
  type1 = 0
  type2 = 0
  falseP= 0
  
  # handle two cases. Extracting included coefficient names from either a glmnet model or an lm model.
  if (method %in% c("LASSO1","LASSO2","ENET1","ENET2")) {
    # Extract coefficients for a specific lambda
    coeffs <- coef(model)
    # Get only the names of nonzero coefficients (excluding intercept)
    select_var <- rownames(coeffs)[which(coeffs != 0)]
    select_var <- select_var[select_var != "(Intercept)"]
  }
  if(method %in% c("AIC","BIC","P-Value")) {
    model_tidy <- broom::tidy(model)
    select_var <- model_tidy$term[-1]
  }
  new_X <- subset(X, select = select_var)
  df <- as.data.frame(cbind(y, new_X))
  new_model <- lm(y ~ ., data = df)
  #refit model with selected variables
  #new_model <- lm(y ~ new_X)
  tidy_model <- tidy(new_model)
  sig_var   <- tidy_model[tidy_model$p.value < .05 & tidy_model$term != "(Intercept)",]$term 
  
  #instantiate biases and coverage
  bias <- rep(0, length(true_var))
  coverage <- rep(0, length(true_var))
  type2 <- rep(0, length(true_var))
  for (i in 1: length(true_var)) {
    #calculate bias and coverage for each true_var item
    if (true_var[i] %in% select_var)
    {
      bias[i] = coef(new_model)[true_var[i]] - true_coef[i]
      coverage[i] = as.numeric( (true_coef[i] >= confint(new_model)[true_var[i],][1]) & 
                                  (true_coef[i] <= confint(new_model)[true_var[i],][2]))
    }
    
    else {
      bias[i] = -true_coef[i] # in the case that the true variable is not in the model, we treat the new model's fitted coefficient as 0
      # coverage stays 0 because our confidence interval only contains 0.
    }
    if(!(true_var %in% sig_var)) {
      type2[i]=1
    }

  }
 
  if (all(true_var %in% select_var)) { 
    trueP = 1
  }
  if (any(!(select_var %in% true_var))) { 
    falseP = 1
  }

  if(!all(sig_var %in% true_var)) {
    type1 = length(setdiff(sig_var, true_var)
  }
  
  type1_var <- setdiff(select_var, true_var)
  
  type1_var_bias <- 0
  
  if (length(type1_var) != 0){
  
  for (i in 1:length(type1_var)) {
    type1_var_bias = type1_var_bias + tidy_model[tidy_model$term == type1_var[i],]$estimate
      
  }
  type1_var_bias = type1_var_bias / 15
  }
  
  
  results <- data.frame(var = true_var, bias = bias, coverage = coverage,
                        type1 = type1, type2 = type2, trueP = trueP, falseP = falseP,
                        method = method, case = case, VE = type1_var_bias)
  
  return(results)
  
}

selection_simulation <- function(n,p,p1,beta,corr,rho, case) {
  #This function takes input for parameters of data generation mechanism
  #It generates the data, performs modelling, and returns metrics of interest.
  sim_results <- data.frame(var = NULL, bias = NULL, coverage = NULL,
                            type1 = NULL, type2 = NULL, trueP = NULL, falseP = NULL,
                            method = NULL, case = NULL)
  simdat <- gen_data(
    n=n,
    p=p,
    p1 =p1,
    beta = beta,
    family = "gaussian",
    signal = "homogeneous",
    corr = corr,
    rho = rho
  )
  df <- as.data.frame(cbind(y = simdat$y, simdat$X))
  lm_mod <- lm(y ~ ., data = df)

  #perform variable selection
  
  #p-values
#takes lm object
  p_model_run <- ols_step_backward_p(
    model = lm_mod,
    p_val = 0.1,
    include = NULL,
    exclude = NULL,
    hierarchical = FALSE,
    progress = FALSE,
    details = FALSE
  )
  p_model <- p_model_run$model
  
  AIC_model <- step(object=lm_mod, method="backward", k=2, trace = 0)
  
  BIC_model <- step(object=lm_mod, method="backward", k=log(nrow(df)), trace = 0)
  
  cv_lasso <- cv.glmnet(x=simdat$X, y=simdat$y, alpha = 1)
  
  lasso_model1 <- glmnet(x=simdat$X, y=simdat$y, alpha = 1, lambda = cv_lasso$lambda.min)
  lasso_model2 <- glmnet(x=simdat$X, y=simdat$y, alpha = 1, lambda = cv_lasso$lambda.1se)
  
  
  cv_enet <- cv.glmnet(x=simdat$X, y=simdat$y, alpha = .5)
  
  # Set alpha to ????? for elastic net model...
  enet_model1 <- glmnet(x=simdat$X, y=simdat$y, alpha = .5, lambda = cv_enet$lambda.1se)
  enet_model2 <- glmnet(x=simdat$X, y=simdat$y, alpha = .5, lambda = cv_enet$lambda.min)
  

  sim_results <- rbind(sim_results, 
    selection_diag(
    model=p_model,
    true_var=true_var,
    true_coef = true_coef,
    X=simdat$X,
    y=simdat$y,
    method="P-Value",
    case = case
  ), 
  selection_diag(
    model=AIC_model,
    true_var=true_var,
    true_coef = true_coef,
    X=simdat$X,
    y=simdat$y,
    method="AIC",
    case = case
  ), 
  selection_diag(
    model=BIC_model,
    true_var=true_var,
    true_coef = true_coef,
    X=simdat$X,
    y=simdat$y,
    method="BIC",
    case = case
  ), 
  selection_diag(
    model=lasso_model1,
    true_var=true_var,
    true_coef = true_coef,
    X=simdat$X,
    y=simdat$y,
    method="LASSO1",
    case = case
  ), 
  selection_diag(
    model=lasso_model2,
    true_var=true_var,
    true_coef = true_coef,
    X=simdat$X,
    y=simdat$y,
    method="LASSO2",
    case = case
  ),
  selection_diag(
    model=enet_model1,
    true_var=true_var,
    true_coef = true_coef,
    X=simdat$X,
    y=simdat$y,
    method="ENET1",
    case = case
  ),
  selection_diag(
    model=enet_model2,
    true_var=true_var,
    true_coef = true_coef,
    X=simdat$X,
    y=simdat$y,
    method="ENET2",
    case = case
  ))
  return(sim_results)
}
  

system.time( { 
#Run the simulation
for (i in 1:10000) {
  print(i)
  
  # simdat1 <-gen_data(
  #   n=250,
  #   p=20,
  #   p1 = 5,
  #   beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),
  #   family = "gaussian",
  #   SNR = 1,
  #   signal = "homogeneous",
  #   rho = c(0)
  # )
  sim_results <- rbind(sim_results,
                      
   selection_simulation(n=250,p=20,p1=5,
                       beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),corr="exchangeable",
                       rho=0, case = "1a"),
  
  selection_simulation(n=250,p=20,p1=5,
                       beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),corr="exchangeable",
                       rho=.35, case = "1b"),
  
  selection_simulation(n=250,p=20,p1=5,
                       beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),corr="exchangeable",
                       rho=.7, case = "1c"),
  
  selection_simulation(n=500,p=20,p1=5,
                       beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),corr="exchangeable",
                       rho=0, case = "2a"),
  
  selection_simulation(n=500,p=20,p1=5,
                       beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),corr="exchangeable",
                       rho=0, case = "2b"),
  
  selection_simulation(n=500,p=20,p1=5,
                       beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),corr="exchangeable",
                       rho=0, case = "2c"))
}
  
  
})

#############
# SUMMARIES #
#############

sim_results01 <- data.frame(var = NULL, bias = NULL, coverage = NULL,
                            type1 = NULL, type2 = NULL, trueP = NULL, falseP = NULL,
                            method = NULL, case = NULL)
for (i in 1:100) {
  print(i)
  
  # simdat1 <-gen_data(
  #   n=250,
  #   p=20,
  #   p1 = 5,
  #   beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),
  #   family = "gaussian",
  #   SNR = 1,
  #   signal = "homogeneous",
  #   rho = c(0)
  # )
  

  
  
  sim_results01 <- rbind(sim_results01,
                       
                       selection_simulation(n=250,p=20,p1=5,
                                            beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),corr="exchangeable",
                                            rho=0, case = "1a"),
                       
                       selection_simulation(n=250,p=20,p1=5,
                                            beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),corr="exchangeable",
                                            rho=.35, case = "1b"),
                       
                       selection_simulation(n=250,p=20,p1=5,
                                            beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),corr="exchangeable",
                                            rho=.7, case = "1c"),
                       
                       selection_simulation(n=500,p=20,p1=5,
                                            beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),corr="exchangeable",
                                            rho=0, case = "2a"),
                       
                       selection_simulation(n=500,p=20,p1=5,
                                            beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),corr="exchangeable",
                                            rho=.35, case = "2b"),
                       
                       selection_simulation(n=500,p=20,p1=5,
                                            beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),corr="exchangeable",
                                            rho=.7, case = "2c"))
}







library(doParallel)
library(foreach)
library(doRNG)

# Set up cluster - leave 1-2 cores free for your OS
n_cores <- detectCores() - 1
library(doSNOW)
cl <- makeCluster(n_cores)
registerDoSNOW(cl)
pb <- txtProgressBar(max=10000, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
#registerDoParallel(cl)

beta_vec <- c(.5/3, 1/3, 1.5/3, 2/3, 2.5/3, rep(0, 15))
set.seed(23192319)
sim_results <- foreach(i = 1:10000, 
                       .combine = rbind, 
                       .export = c("selection_simulation", "selection_diag"),
                       .packages = c("hdrm", "olsrr","dplyr","glmnet"),
                       .options.snow = opts
                       ) %dopar% {
                         rbind(
                           selection_simulation(n=250, p=20, p1=5, beta=beta_vec, corr="exchangeable", rho=0,   case="1a"),
                           selection_simulation(n=250, p=20, p1=5, beta=beta_vec, corr="exchangeable", rho=.35, case="1b"),
                           selection_simulation(n=250, p=20, p1=5, beta=beta_vec, corr="exchangeable", rho=.7,  case="1c"),
                           selection_simulation(n=500, p=20, p1=5, beta=beta_vec, corr="exchangeable", rho=0,   case="2a"),
                           selection_simulation(n=500, p=20, p1=5, beta=beta_vec, corr="exchangeable", rho=.35,   case="2b"),
                           selection_simulation(n=500, p=20, p1=5, beta=beta_vec, corr="exchangeable", rho=.7,   case="2c")
                         )
                       }

stopCluster(cl)




old_pkgs <- list.dirs("/Library/Frameworks/R.framework/Versions/4.6/Resources/library_old", 
                      recursive = FALSE, full.names = FALSE)

# Remove base packages (already installed with R)
base_pkgs <- rownames(installed.packages(priority = "base"))
recommended_pkgs <- rownames(installed.packages(priority = "recommended"))
to_install <- setdiff(old_pkgs, c(base_pkgs, recommended_pkgs))

install.packages(to_install)
























library(hdrm)
library(doSNOW)
library(foreach)
library(doRNG)

# Setup
true_var  <- c("V01", "V02", "V03", "V04", "V05")
true_coef <- seq(.5/3, 2.5/3, .5/3)
beta_vec  <- c(.5/3, 1/3, 1.5/3, 2/3, 2.5/3, rep(0, 15))

n_cores <- detectCores() - 2
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

pb <- txtProgressBar(min=0, max=10000, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

system.time({
  set.seed(23192319)
  sim_results <- foreach(i = 1:10000,
                         .combine = rbind,
                         .packages = c("hdrm", "glmnet", "olsrr", "broom"),
                         .export  = c("selection_simulation", "selection_diag",
                                      "true_var", "true_coef"),
                         .options.snow = opts) %dorng% {
                           rbind(
                             selection_simulation(n=250, p=20, p1=5, beta=beta_vec, corr="exchangeable", rho=0,   case="1a"),
                             selection_simulation(n=250, p=20, p1=5, beta=beta_vec, corr="exchangeable", rho=.35, case="1b"),
                             selection_simulation(n=250, p=20, p1=5, beta=beta_vec, corr="exchangeable", rho=.7,  case="1c"),
                             selection_simulation(n=500, p=20, p1=5, beta=beta_vec, corr="exchangeable", rho=0,   case="2a"),
                             selection_simulation(n=500, p=20, p1=5, beta=beta_vec, corr="exchangeable", rho=.35,   case="2b"),
                             selection_simulation(n=500, p=20, p1=5, beta=beta_vec, corr="exchangeable", rho=.7,   case="2c")
                           )
                         }
})

close(pb)
stopCluster(cl)


system.time(
  selection_simulation(n=250, p=20, p1=5, beta=beta_vec, 
                       corr="exchangeable", rho=0, case="1a")
)


