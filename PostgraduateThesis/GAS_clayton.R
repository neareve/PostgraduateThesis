clayton_pdf <- function(u, v, theta) {
  pdf <- (u * v) ^ (-1 - theta) * (u ^ (-theta) + v ^ (-theta) - 1) ^ (-1 / theta - 2) * (1 + theta)
  return(pdf)
}

claytonCL <- function(theta, data) {
  u <- data[, 1]
  v <- data[, 2]
  out1 <- (u * v) ^ (-1 - theta)
  out2 <- (u ^ (-theta) + v ^ (-theta) - 1) ^ (-1 / theta - 2)
  out <- (1 + theta) * out1 * out2
  CL <- sum(log(out))
  if (is.na(CL)) {
    CL <- -1e8
  }
  return(CL = -CL)
}

claytonCLa <- function(theta, data) {
  u <- data[, 1]
  v <- data[, 2]
  out1 <- (u * v) ^ (-1 - theta)
  out2 <- (u ^ (-theta) + v ^ (-theta) - 1) ^ (-1 / theta - 2)
  out <- (1 + theta) * out1 * out2
  CLa <- log(out)
  return(CLa = -CLa)
}

# ������ֵ����
LLgrad_c <- function(theta, data) {
  l <- length(data[, 1])
  B <- matrix(NA, nrow = l, ncol = 2)
  B[, 1] <- claytonCLa(theta, data)
  h <- 2.2204e-16 ^ (1 / 3) * max(abs(theta), 0.01) #�൱�ڵ��������ȡ΢Сֵ
  thetad <- theta - h
  B[, 2] <- claytonCLa(thetad, data)
  B2 <- (B[, 1] - B[, 2]) / h
  return(grad = B2)
}

##########################################################
reps <- 1e4 # ����hessian����
# K<-matrix(seq(.01,20,length.out = 100),ncol = 1)
K <- matrix(c(seq(-.99, -.01, length.out = 20), seq(.01, 20, length.out = 80)), ncol = 1)
HESSclayton <- matrix(NA, nrow = length(K), ncol = 1)
for (i in 1:length(K)) {
  type <- claytonCopula(K[i], dim = 2)
  Usim <- rCopula(reps, type)
  scoreSIM <- LLgrad_c(K[i], Usim) #�൱����delta
  HESSclayton[i] <- mean(scoreSIM ^ 2, na.rm = T) #�൱����deltaƽ��
}
hessINFO3 <- cbind(K, HESSclayton)
# ʱ�亯��
clayton_GAS_CL <- function(theta, data, hessINFO) {
  RBAR <- 0.9999
  t <- nrow(data)
  w <- theta[1]
  a <- theta[2]
  b <- theta[3]
  h <- 0.00001 #����÷����ò���
  # ����rhotʱ������
  ft <- rep(0, t)
  rho <- rep(0, t)
  rho[1] <- rho0
  ft[1] <- log(1 + rho[1])
  for (i in 2:t) {
    It <- ml_spline(hessINFO[, 1], hessINFO[, 2], rho[i - 1])
    DELTAt <- (-claytonCL(rho[i - 1] + h, t(data[i - 1,])) - -claytonCL(rho[i - 1], t(data[i - 1,]))) / h
    drhodf <- -exp(ft[i - 1])
    Itildat <- It / (drhodf ^ 2)
    DELTAtildat <- DELTAt / drhodf
    ft[i] <- w + b * ft[i - 1] + a * DELTAtildat / sqrt(Itildat)
    rho[i] <- -1 + exp(ft[i])
  }
  return(rho)
}
# ������Ȼ����
clayton_GAS_maxlik <- function(theta, data, hessINFO) {
  rho <- clayton_GAS_CL(theta, data, hessINFO)
  CL <- claytonCL(rho, data)
}