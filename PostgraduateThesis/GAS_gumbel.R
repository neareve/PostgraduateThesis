##################定义分布和密度函数形式###################
gumbel_cdf <- function(u, v, k1) {
  cdf <- exp(-((-log(u)) ^ k1 + (-log(v)) ^ k1) ^ (k1 ^ (-1)))
  return(cdf)
}

gumbel_pdf <- function(u, v, k1) {
  pdf <- gumbel_cdf(u, v, k1) * (u * v) ^ (-1) * (((-log(u)) * (-log(v))) ^ (k1 - 1))
  pdf <- pdf * (((-log(u)) ^ k1 + (-log(v)) ^ k1) ^ (k1 ^ (-1) - 2))
  pdf <- pdf * ((((-log(u)) ^ k1 + (-log(v)) ^ k1) ^ (k1 ^ (-1))) + k1 - 1)
  return(pdf)
}

gumbelCL <- function(theta, data) {
  ut <- -log(data[, 1])
  vt <- -log(data[, 2])
  CL <- log(gumbel_cdf(data[, 1], data[, 2], theta)) - log(data[, 1]) - log(data[, 2])
  CL <- CL + (theta - 1) * (log(ut) + log(vt)) - (2 - 1 / theta) * (log(ut ^ theta + vt ^ theta))
  CL <- CL + log((ut ^ theta + vt ^ theta) ^ (1 / theta) + theta - 1)
  CL <- sum(CL)
  if (is.na(CL)) {
    CL <- -1e8
  }
  return(CL = -CL)
}

gumbelCLa <- function(theta, data) {
  ut <- -log(data[, 1])
  vt <- -log(data[, 2])
  CL <- log(gumbel_cdf(data[, 1], data[, 2], theta)) - log(data[, 1]) - log(data[, 2])
  CL <- CL + (theta - 1) * (log(ut) + log(vt)) - (2 - 1 / theta) * (log(ut ^ theta + vt ^ theta))
  CL <- CL + log((ut ^ theta + vt ^ theta) ^ (1 / theta) + theta - 1)
  return(CLa = -CL)
}

#################定义格点样条插值函数####################
LLgrad_g <- function(theta, data) {
  l <- length(data[, 1])
  B <- matrix(NA, nrow = l, ncol = 2)
  B[, 1] <- gumbelCLa(theta, data)
  h <- 2.2204e-16 ^ (1 / 3) * max(abs(theta), 0.01) #相当于导数定义的取微小值
  thetad <- theta - h
  B[, 2] <- gumbelCLa(thetad, data)
  B2 <- (B[, 1] - B[, 2]) / h
  return(grad = B2)
}
###########################################################
# 生成hessian矩阵
reps <- 1e4
K <- matrix(seq(1.01, 20, length.out = 100), ncol = 1)
HESSgumbel <- matrix(NA, nrow = length(K), ncol = 1)
for (i in 1:length(K)) {
  type <- gumbelCopula(K[i], dim = 2)
  Usim <- rCopula(reps, type)
  scoreSIM <- LLgrad_g(K[i], Usim) #相当于求delta
  HESSgumbel[i] <- mean(scoreSIM ^ 2, na.rm = T) #相当于求delta平方
}
hessINFO4 <- cbind(K, HESSgumbel)
###################定义时变极大似然函数###############
gumbel_GAS_CL <- function(theta, data, hessINFO) {
  t <- nrow(data)
  w <- theta[1]
  a <- theta[2]
  b <- theta[3]
  h <- 0.0001 #计算得分时的步长
  ft <- rep(NA, t)
  rho <- rep(NA, t) #gumbel copula的初始值，可选用常系数模型的估计值
  rho[1] <- rho0
  ft[1] <- log(rho[1] - 1) #根据转换函数反函数得到
  for (i in 2:t) {
    It <- ml_spline(hessINFO[, 1], hessINFO[, 2], rho[i - 1])
    DELTAt <- (-gumbelCL(rho[i - 1] + h, t(data[i - 1,])) - -gumbelCL(rho[i - 1], t(data[i - 1,]))) / h #估计kappa的得分矩阵
    drhodf <- exp(ft[i - 1])
    Itildat <- It / (drhodf ^ 2) #根据kappa和f的函数与反函数求导关系求f的hessian矩阵，注意不同转换函数具有不同形式导数关系
    DELTAtildat <- DELTAt / drhodf #根据kappa和f的函数与反函数求导关系求f的得分矩阵
    ft[i] <- w + b * ft[i - 1] + a * DELTAtildat / sqrt(Itildat)
    rho[i] <- 1 + exp(ft[i])
  }
  return(rho)
}

###################给出极大似然函数表达式###################
gumbel_GAS_maxlik <- function(theta, data, hessINFO) {
  rho <- gumbel_GAS_CL(theta, data, hessINFO)
  CL <- gumbelCL(rho, data)
}
