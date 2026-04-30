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
  
  #instantiate biases and coverage
  bias <- rep(0, length(true_var))
  coverage <- rep(0, length(true_var))
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
  }
  tidy_model <- broom::tidy(new_model)
  sig_var   <- tidy_model[tidy_model$p.value < .05 & tidy_model$term != "(Intercept)",]$term 
  
  if (all(true_var %in% select_var)) { 
    trueP = 1
  }
  if (any(!(select_var %in% true_var))) { 
    falseP = 1
  }
  if (!all(true_var %in% sig_var)) {
    type2 = 1
  }
  if(!all(sig_var %in% true_var)) {
    type1 = 1
  }
  
  results <- data.frame(var = true_var, bias = bias, coverage = coverage,
                        type1 = type1, type2 = type2, trueP = trueP, falseP = falseP,
                        method = method, case = case)
  
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
  lm_mod <- lm(y ~ ., data = df1)

  #perform variable selection
  
  #p-values
#takes lm object
  p_model_run <- ols_step_backward_p(
    model = lm_mod,
    p_val = 0.15,
    include = NULL,
    exclude = NULL,
    hierarchical = FALSE,
    progress = FALSE,
    details = FALSE
  )
  p_model <- p_model_run$model
  
  AIC_model <- step(object=lm_mod1, method="backward", k=2)
  
  BIC_model <- step(object=lm_mod, method="backward", k=log(nrow(df)))
  
  cv_lasso <- cv.glmnet(x=simdat$X, y=simdat$y, alpha = 1)
  
  lasso_model1 <- glmnet(x=simdat$X, y=simdat$y, alpha = 1, lambda = cv_lasso$lambda.min)
  lasso_model2 <- glmnet(x=simdat$X, y=simdat$y, alpha = 1, lambda = cv_lasso$lambda.1se)
  
  
  cv_enet <- cv.glmnet(x=simdat$X, y=simdat$y, alpha = .5)
  
  # Set alpha to ????? for elastic net model...
  enet_model1 <- glmnet(x=simdat$X, y=simdat$y, alpha = .5, lambda = cv_enet$lambda.1se)
  enet_model2 <- glmnet(x=simdat$X, y=simdat$y, alpha = .5, lambda = cv_enet$lambda.min)
  

  sim_results <- rbind(sim_results, selection_diag(
    model=p_model,
    true_var=true_var,
    true_coef = true_coef,
    X=simdat$X,
    y=simdat$y,
    method="P-Value",
    case = case
  ))
  
  sim_results <- rbind(sim_results, selection_diag(
    model=AIC_model,
    true_var=true_var,
    true_coef = true_coef,
    X=simdat$X,
    y=simdat$y,
    method="AIC",
    case = case
  ))
  
  sim_results <- rbind(sim_results, selection_diag(
    model=BIC_model,
    true_var=true_var,
    true_coef = true_coef,
    X=simdat$X,
    y=simdat$y,
    method="BIC",
    case = case
  ))
  
  sim_results <- rbind(sim_results, selection_diag(
    model=lasso_model1,
    true_var=true_var,
    true_coef = true_coef,
    X=simdat$X,
    y=simdat$y,
    method="LASSO1",
    case = case
  ))
  
  sim_results <- rbind(sim_results, selection_diag(
    model=lasso_model2,
    true_var=true_var,
    true_coef = true_coef,
    X=simdat$X,
    y=simdat$y,
    method="LASSO2",
    case = case
  ))
  
  sim_results <- rbind(sim_results, selection_diag(
    model=enet_model1,
    true_var=true_var,
    true_coef = true_coef,
    X=simdat$X,
    y=simdat$y,
    method="ENET1",
    case = case
  ))
  
  sim_results <- rbind(sim_results, selection_diag(
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
  


#Run the simulation
for (i in 1:10) {
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
  rho=0, case = "1a"))
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









