# ����hessian����
reps <- 1e4
K <- matrix(seq(1.01, 20, length.out = 100), ncol = 1)
HESSgumbel <- matrix(NA, nrow = length(K), ncol = 1)
for (i in 1:length(K)) {
  type <- gumbelCopula(K[i], dim = 2)
  Usim <- rCopula(reps, type)
  Usim <- 1 - Usim
  scoreSIM <- LLgrad_g(K[i], Usim) #�൱����delta
  HESSgumbel[i] <- mean(scoreSIM ^ 2, na.rm = T) #�൱����deltaƽ��
}
hessINFO14 <- cbind(K, HESSgumbel)
###################����ʱ�伫����Ȼ����###############
rotgumbel_GAS_CL <- function(theta, data, hessINFO) {
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
rotgumbel_GAS_maxlik = function(theta, data, hessINFO) {
  data <- 1 - data
  rho <- rotgumbel_GAS_CL(theta, data, hessINFO)
  CL <- gumbelCL(rho, data)
}