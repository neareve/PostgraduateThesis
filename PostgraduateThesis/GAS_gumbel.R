##################����ֲ����ܶȺ�����ʽ###################
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

#################������������ֵ����####################
LLgrad_g <- function(theta, data) {
  l <- length(data[, 1])
  B <- matrix(NA, nrow = l, ncol = 2)
  B[, 1] <- gumbelCLa(theta, data)
  h <- 2.2204e-16 ^ (1 / 3) * max(abs(theta), 0.01) #�൱�ڵ��������ȡ΢Сֵ
  thetad <- theta - h
  B[, 2] <- gumbelCLa(thetad, data)
  B2 <- (B[, 1] - B[, 2]) / h
  return(grad = B2)
}
###########################################################
# ����hessian����
reps <- 1e4
K <- matrix(seq(1.01, 20, length.out = 100), ncol = 1)
HESSgumbel <- matrix(NA, nrow = length(K), ncol = 1)
for (i in 1:length(K)) {
  type <- gumbelCopula(K[i], dim = 2)
  Usim <- rCopula(reps, type)
  scoreSIM <- LLgrad_g(K[i], Usim) #�൱����delta
  HESSgumbel[i] <- mean(scoreSIM ^ 2, na.rm = T) #�൱����deltaƽ��
}
hessINFO4 <- cbind(K, HESSgumbel)
###################����ʱ�伫����Ȼ����###############
gumbel_GAS_CL <- function(theta, data, hessINFO) {
  t <- nrow(data)
  w <- theta[1]
  a <- theta[2]
  b <- theta[3]
  h <- 0.0001 #����÷�ʱ�Ĳ���
  ft <- rep(NA, t)
  rho <- rep(NA, t) #gumbel copula�ĳ�ʼֵ����ѡ�ó�ϵ��ģ�͵Ĺ���ֵ
  rho[1] <- rho0
  ft[1] <- log(rho[1] - 1) #����ת�������������õ�
  for (i in 2:t) {
    It <- ml_spline(hessINFO[, 1], hessINFO[, 2], rho[i - 1])
    DELTAt <- (-gumbelCL(rho[i - 1] + h, t(data[i - 1,])) - -gumbelCL(rho[i - 1], t(data[i - 1,]))) / h #����kappa�ĵ÷־���
    drhodf <- exp(ft[i - 1])
    Itildat <- It / (drhodf ^ 2) #����kappa��f�ĺ����뷴�����󵼹�ϵ��f��hessian����ע�ⲻͬת���������в�ͬ��ʽ������ϵ
    DELTAtildat <- DELTAt / drhodf #����kappa��f�ĺ����뷴�����󵼹�ϵ��f�ĵ÷־���
    ft[i] <- w + b * ft[i - 1] + a * DELTAtildat / sqrt(Itildat)
    rho[i] <- 1 + exp(ft[i])
  }
  return(rho)
}

###################����������Ȼ��������ʽ###################
gumbel_GAS_maxlik <- function(theta, data, hessINFO) {
  rho <- gumbel_GAS_CL(theta, data, hessINFO)
  CL <- gumbelCL(rho, data)
}