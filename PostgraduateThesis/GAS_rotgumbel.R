# 生成hessian矩阵
reps <- 1e4
K <- matrix(seq(1.01, 20, length.out = 100), ncol = 1)
HESSgumbel <- matrix(NA, nrow = length(K), ncol = 1)
for (i in 1:length(K)) {
  type <- gumbelCopula(K[i], dim = 2)
  Usim <- rCopula(reps, type)
  Usim <- 1 - Usim
  scoreSIM <- LLgrad_g(K[i], Usim) #相当于求delta
  HESSgumbel[i] <- mean(scoreSIM ^ 2, na.rm = T) #相当于求delta平方
}
hessINFO14 <- cbind(K, HESSgumbel)
###################定义时变极大似然函数###############
rotgumbel_GAS_CL <- function(theta, data, hessINFO) {
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
rotgumbel_GAS_maxlik = function(theta, data, hessINFO) {
  data <- 1 - data
  rho <- rotgumbel_GAS_CL(theta, data, hessINFO)
  CL <- gumbelCL(rho, data)
}
