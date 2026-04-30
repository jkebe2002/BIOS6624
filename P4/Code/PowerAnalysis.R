# Here I investigated bias and coverage simulation sample sizes and monte carlo standard errors. 
# I found that 10000 is smaller than the coverage neee

for (i in seq(100,10000, 100)) {
  for(j in 1:i)
    ests <- seq(0,i,1)
    sds <- seq(0,i,1)
  
  case_1a <- gen_data(
    n=250,
    p=20,
    p1 = 5,
    beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),
    family = "gaussian",
    SNR = 1,
    signal = "homogeneous",
    corr = c("exchangeable"),
    rho = .7
  )
  df <- as.data.frame(cbind(y = simdat$y, simdat$X))
  lm_mod <- lm(y ~ ., data = df)
  AIC_model <- step(object=lm_mod, method="backward", k=2, trace = 0)
  tidy_AIC <- broom::tidy(AIC_model)
  if ("V05" %in% tidy_AIC$term) {
    ests[j] <- tidy_AIC[tidy_AIC$term=="V05",]$estimate
    sds[j] <- tidy_AIC[tidy_AIC$term=="V05",]$std.error
  }
  
  print(i)
  print(sqrt( (1/(i*(i-1))) * sd(ests)) )
  print("percent:")
  print(100* sqrt( (1/(i*(i-1))) * sd(ests)) / mean(sds) )
  print(sqrt( (1/(i*(i-1))) * sd(ests)) < .1*mean(sds))
  

  
}




# Monte Carlo SE for estimating coverage


for (i in seq(100,10000, 100)) {
  
  print("Coverage MC SE")
  print(sqrt(95*5/i))
  
}


95*5 / (.1)^2

