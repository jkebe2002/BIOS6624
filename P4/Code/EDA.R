simdat1 <-gen_data(
  n=250,
  p=20,
  p1 = 5,
  beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),
  family = "gaussian",
  SNR = 1,
  signal = "homogeneous",
  rho = c(0)
)

simdat2 <-gen_data(
  n=250,
  p=20,
  p1 = 5,
  beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),
  family = "gaussian",
  SNR = 1,
  signal = "homogeneous",
  corr = c("exchangeable"),
  rho = .35
)

simdat3 <-gen_data(
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

simdat4 <-gen_data(
  n=500,
  p=20,
  p1 = 5,
  beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),
  family = "gaussian",
  SNR = 1,
  signal = "homogeneous",
  rho = c(0)
)

simdat5<-gen_data(
  n=500,
  p=20,
  p1 = 5,
  beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),
  family = "gaussian",
  SNR = 1,
  signal = "homogeneous",
  corr = c("exchangeable"),
  rho = .35
)

simdat6 <-gen_data(
  n=500,
  p=20,
  p1 = 5,
  beta = c(.5/3,1/3,1.5/3,2/3,2.5/3, rep(0, 15)),
  family = "gaussian",
  SNR = 1,
  signal = "homogeneous",
  corr = c("exchangeable"),
  rho = .7
)











