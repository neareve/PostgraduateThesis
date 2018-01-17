library(copula)
library(VineCopula)
library(nloptr)
library(Rcpp)
sourceCpp("MlSpline.cpp")
source("GAS_clayton.R") #3
source("GAS_frank.R") #5
source("GAS_gumbel.R") #4
source("GAS_joe.R") #6
source("GAS_normal.R") #1
source("GAS_t.R") #2
source("ARMA_Gumbel.R")
source("ARMA_Normal.R")
source("ARMA_Joe.R")
source("ARMA_t.R")
source("ARMA_Clayton.R")
source("ARMA_Frank.R")

# simulation
######################## ncopula, 1###########################################
l <- 1000
par1 <- c(.2 + runif(n = l/2, min = -.2, max = .2), .6 + runif(n = l/2, -.2, .2))
plot(par1, type = "l", xlab = "", ylab = "", ylim = c(min(par1) - .1, max(par1) + .2))
data1 <- matrix(NA, ncol = 2, nrow = l)
for(i in 1:l){
  cop <- normalCopula(par1[i], dim = 2)
  if(i <= l/2){
    data1[i,] <- rCopula(1, cop)  
  }else{
    data1[i,] <- rCopula(1, cop)
  }
}

fit1 <- BiCopSelect(data1[,1], data1[,2], familyset = 1,indeptest = T)
abline(h = fit1$par, col = "blue")
rho0 <- fit1$par
BIC1_S <- fit1$BIC
AIC1_S <- fit1$AIC
w0<--log((1-rho0)/(1+rho0))
out <- nloptr(x0 = c(w0,.1,.1), eval_f = ncopula_GAS_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
              opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval=10000),
              data = data1, hessINFO = hessINFO1)
theta1 <- out$solution
rho1 <- ncopula_GAS_CL(theta1, data1, hessINFO1)
BIC1_GAS <- log(l)*3 - 2 * (-ncopula_GAS_maxlik(theta1, data1, hessINFO1))
AIC1_GAS <- 2 * 3 - 2 * (-ncopula_GAS_maxlik(theta1, data1, hessINFO1))
lines(rho1, col = "red")

out_ARMA <- nloptr(x0 =c(w0,.1,.1), eval_f = ncopula_ARMA_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
                opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval=10000),
                data = data1)
theta1_ARMA <- out_ARMA$solution
rho1_ARMA <- ncopula_ARMA_CL(theta1_ARMA, data1)
BIC1_ARMA <- log(l)*3 - 2 * (-ncopula_ARMA_maxlik(theta1_ARMA, data1))
AIC1_ARMA <- 2 * 3 - 2 * (-ncopula_ARMA_maxlik(theta1_ARMA, data1))
lines(rho1_ARMA, col = "green")

MSE_NStation <- sum((par1 - rho0)^2)/l
MSE_NGAS <- sum((par1 - rho1)^2)/l
MSE_NARMA <- sum((par1 - rho1_ARMA)^2)/l


################################## tcopula, 2#######################################
# hessian matrix
reps<-1e4
R<-matrix(c(seq(-.99,-.5,length=40),seq(-.49,.5,length=20),seq(.51,.99,length=40)),ncol=1)
HESSstudt<-matrix(NA,length(R),ncol=1)
for(k in 1:length(R)){
  type<-tCopula(R[k],df = nu0,dim = 2)
  Usim<-rCopula(reps,type)
  scoreSIM<-LLgrad_t(R[k],Usim)
  HESSstudt[k]<-mean(scoreSIM^2,na.rm = T)
} 
hessINFO2<-cbind(R,HESSstudt[,1]) #唯一确定t的随机矩阵位置


par2 <- c(.2 + runif(n = l/2, min = -.2, max = .2), .7 + runif(n = l/2, -.2, .2))
plot(par2, type = "l", xlab = "", ylab = "")
data2 <- matrix(NA, ncol = 2, nrow = l)
for(i in 1:l){
  cop <- tCopula(par2[i], dim = 2)
  if(i <= l/2){
    data2[i,] <- rCopula(1, cop)  
  }else{
    data2[i,] <- rCopula(1, cop)
  }
}

fit2 <- BiCopSelect(data2[,1], data2[,2], familyset = 2,indeptest = T)
abline(h = fit2$par, col = "blue")
rho0 <- fit2$par
nu0 <- fit2$par2
BIC2_S <- fit2$BIC
AIC2_S <- fit2$AIC
w0 <- log((1+rho0)/(1-rho0))

# opt
out <- nloptr(x0 =c(w0,.1,.1), eval_f =tcopula_GAS_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
              opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval=10000),
              data = data2, hessINFO = hessINFO2)
theta2 <- out$solution
rho2 <- tcopula_GAS_CL(theta2, data2, hessINFO2)
BIC2_GAS <- log(l)*3 - 2 * (-tcopula_GAS_maxlik(theta2, data2, hessINFO2))
AIC2_GAS <- 2 * 3 - 2 * (-tcopula_GAS_maxlik(theta2, data2, hessINFO2))
lines(rho2, col = "red")

out_ARMA <- nloptr(x0 =c(w0,.1,.1), eval_f = tcopula_ARMA_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
                opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval=10000),
                data = data2)
theta2_ARMA <- out_ARMA$solution
rho2_ARMA <- tcopula_ARMA_CL(theta2_ARMA, data2)
BIC2_ARMA <- log(l)*3 - 2 * (-tcopula_ARMA_maxlik(theta2_ARMA, data2))
AIC2_ARMA <- 2 * 3 - 2 * (-tcopula_ARMA_maxlik(theta2_ARMA, data2))
lines(rho2_ARMA, col = "green")

MSE_tStation <- sum((par2 - rho0)^2)/l
MSE_tGAS <- sum((par2 - rho2)^2)/l
MSE_tARMA <- sum((par2 - rho2_ARMA)^2)/l

###################################### clayton, 3#################################
par3 <- c(.5 + runif(n = l/2, min = -.2, max = .2), 1 + runif(n = l/2, -.2, .2))
plot(par3, type = "l", xlab = "", ylab = "")
data3 <- matrix(NA, ncol = 2, nrow = l)
for(i in 1:l){
  cop <- claytonCopula(par3[i], dim = 2)
  if(i <= l/2){
    data3[i,] <- rCopula(1, cop)  
  }else{
    data3[i,] <- rCopula(1, cop)
  }
}

fit3 <- BiCopSelect(data3[,1], data3[,2], familyset = 3, indeptest = T, rotations = F)
abline(h = fit3$par, col = "blue")
rho0 <- fit3$par
BIC3_S <- fit3$BIC
AIC3_S <- fit3$AIC
w0 <- log(1 + rho0)

out <- nloptr(x0 =c(w0,.1,.1), eval_f =clayton_GAS_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
              opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval=10000),
              data = data3, hessINFO = hessINFO3)
theta3 <- out$solution
rho3 <- clayton_GAS_CL(theta3, data3, hessINFO3)
BIC3_GAS <- log(l)*3 - 2 * (-clayton_GAS_maxlik(theta3, data3, hessINFO3))
AIC3_GAS <- 2 * 3 - 2 * (-clayton_GAS_maxlik(theta3, data3, hessINFO3))
lines(rho3, col = "red")

###
out_ARMA <- nloptr(x0 =c(w0,.1,.1), eval_f = clayton_ARMA_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
                opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval=10000),
                data = data3)
theta3_ARMA <- out_ARMA$solution
rho3_ARMA <- clayton_ARMA_CL(theta3_ARMA, data3)
BIC3_ARMA <- log(l)*3 - 2 * (-clayton_ARMA_maxlik(theta3_ARMA, data3))
AIC3_ARMA <- 2 * 3 - 2 * (-clayton_ARMA_maxlik(theta3_ARMA, data3))
lines(rho3_ARMA, col = "green")

MSE_cStation <- sum((par3 - rho0)^2)/l
MSE_cGAS <- sum((par3 - rho3)^2)/l
MSE_cARMA <- sum((par3 - rho3_ARMA)^2)/l

########################################### gumbel, 4################################################
par4 <- c(2 + runif(l/2), 3 + runif(l/2))
data4 <- matrix(NA, nrow = l, ncol = 2)
plot(par4, type = "l", xlab = "", ylab = "")
for(i in 1:l){
  cop <- gumbelCopula(par4[i],2)
  if(i <= l/2){
    data4[i,] <- rCopula(1, cop)
  }else{
    data4[i,] <- rCopula(1, cop)
  }
}
fit4 <- BiCopSelect(data4[,1], data4[,2], familyset = 4,indeptest = T)
abline(h = fit4$par, col = "blue")
rho0 <- fit4$par
BIC4_S <- fit4$BIC
AIC4_S <- fit4$AIC
w0<-log(rho0-1)
out <- nloptr(x0 =c(w0,.1,.1), eval_f = gumbel_GAS_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
              opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval=10000),
              data = data4, hessINFO = hessINFO4)
theta4<-out$solution
rho4<-gumbel_GAS_CL(theta4, data4, hessINFO4)
bic4_GAS <- log(l)*3 - 2 * (-gumbel_GAS_maxlik(theta4, data4, hessINFO4))
AIC4_GAS <- 2 * 3 - 2 * (-gumbel_GAS_maxlik(theta4, data4, hessINFO4))
lines(rho4, col = "red")

out_ARMA <- nloptr(x0 =c(w0,.1,.1), eval_f = gumbel_ARMA_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
                opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval=10000),
                data = data4)
theta4_ARMA <- out_ARMA$solution
rho4_ARMA <- gumbel_ARMA_CL(theta4_ARMA, data4)
BIC4_ARMA <- log(l)*3 - 2 * (-gumbel_ARMA_maxlik(theta4_ARMA, data4))
AIC4_ARMA <- 2 * 3 - 2 * (-gumbel_ARMA_maxlik(theta4_ARMA, data4))
lines(rho4_ARMA, col = "green")

MSE_gStation <- sum((par4 - rho0)^2)/l
MSE_gGAS <- sum((par4 - rho4)^2)/l
MSE_gARMA <- sum((par4 - rho4_ARMA)^2)/l

############################################ frank, 5############################################
par5 <- c(1 + runif(l/2), 4 + rnorm(l/2))
data5 <- matrix(NA, nrow = l, ncol = 2)
plot(par5, type = "l")
for(i in 1:l){
  cop <- frankCopula(par5[i],2)
  if(i <= l/2){
    data5[i,] <- rCopula(1, cop)
  }else{
    data5[i,] <- rCopula(1, cop)
  }
}

fit5 <- BiCopSelect(data5[,1], data5[,2], familyset = 5,indeptest = T, rotations = F)
abline(h = fit5$par, col = "blue")
rho0 <- fit5$par
BIC5_S <- fit5$BIC
AIC5_S <- fit5$AIC
w0<- 1/rho0 - 1
out <- nloptr(x0 =c(w0,.1,.1), eval_f = frank_GAS_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
              opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval=10000),
              data = data5, hessINFO = hessINFO5)
theta5<-out$solution
rho5 <- frank_GAS_CL(theta5, data5, hessINFO5)
BIC5_GAS <- log(l)*3 - 2 * (-frank_GAS_maxlik(theta5, data5, hessINFO5))
AIC5_GAS <- 2 * 3 - 2 * (-frank_GAS_maxlik(theta5, data5, hessINFO5))
lines(rho5, col = "red")

out_ARMA <- nloptr(x0 =c(w0,.1,.1), eval_f = frank_ARMA_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
                opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval=10000),
                data = data5)
theta5_ARMA <- out_ARMA$solution
rho5_ARMA <- frank_ARMA_CL(theta5_ARMA, data5)
BIC5_ARMA <- log(l)*3 - 2 * (-frank_ARMA_maxlik(theta5_ARMA, data5))
AIC5_ARMA <- 2 * 3 - 2 * (-frank_ARMA_maxlik(theta5_ARMA, data5))
lines(rho5_ARMA, col = "green")

MSE_fStation <- sum((par5 - rho0)^2)/l
MSE_fGAS <- sum((par5 - rho5)^2)/l
MSE_fARMA <- sum((par5 - rho5_ARMA)^2)/l


######################################### joe, 6#################################################
par6 <- c(2 + runif(l/2), 3 + runif(l/2))
data6 <- matrix(NA, nrow = l, ncol = 2)
plot(par6, type = "l")
for(i in 1:l){
  cop <- joeCopula(par6[i],2)
  if(i <= l/2){
    data6[i,] <- rCopula(1, cop)
  }else{
    data6[i,] <- rCopula(1, cop)
  }
}

fit6 <- BiCopSelect(data6[,1], data6[,2], familyset = 6,indeptest = T, rotations = F)
abline(h = fit6$par, col = "blue")
rho0 <- fit6$par
BIC6_S <- fit6$BIC
AIC6_S <- fit6$AIC
w0<-log(rho0-1)
out <- nloptr(x0 =c(w0,.1,.1), eval_f = joe_GAS_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
              opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval=10000),
              data = data6, hessINFO = hessINFO6)
theta6<-out$solution
rho6<-joe_GAS_CL(theta6, data6, hessINFO6)
BIC6_GAS <- log(l)*3 - 2 * (-joe_GAS_maxlik(theta6, data6, hessINFO6))
AIC6_GAS <- 2 * 3 - 2 * (-joe_GAS_maxlik(theta6, data6, hessINFO6))
lines(rho6, col = "red")

out_ARMA <- nloptr(x0 =c(w0,.1,.1), eval_f = joe_ARMA_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
                opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval=10000),
                data = data6)
theta6_ARMA <- out_ARMA$solution
rho6_ARMA <- joe_ARMA_CL(theta6_ARMA, data6)
BIC6_ARMA <- log(l)*3 - 2 * (-joe_ARMA_maxlik(theta6_ARMA, data6))
AIC6_ARMA <- 2 * 3 - 2 * (-joe_ARMA_maxlik(theta6_ARMA, data6))
lines(rho6_ARMA, col = "green")

MSE_jStation <- sum((par6 - rho0)^2)/l
MSE_jGAS <- sum((par6 - rho6)^2)/l
MSE_jARMA <- sum((par6 - rho6_ARMA)^2)/l


###########################plot##############################################
m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)

layout(mat = m, heights = c(0.4,0.4,0.2))
# normal copula
par(mar = c(2, 2, 1, 1))
plot(par1, type = "l", xlab = "", ylab = "",main = "nomal copula",
     ylim = c(min(rho1_ARMA) - .02, max(par1) + .05))
abline(h = fit1$par, col = "blue")
lines(rho1, col = "red")
lines(rho1_ARMA, col = "green")
# tcopula
plot(par2, type = "l", xlab = "", ylab = "", main = "t copula",
     ylim = c(min(rho2_ARMA) - .01, max(rho2_ARMA) + .2))
abline(h = fit2$par, col = "blue")
lines(rho2, col = "red")
lines(rho2_ARMA, col = "green")
# clayton copula
plot(par3, type = "l", xlab = "", ylab = "", main = "clayton copula",
     ylim = c(min(rho3_ARMA) - .05, max(rho3_ARMA) + .05))
abline(h = fit3$par, col = "blue")
lines(rho3, col = "red")
lines(rho3_ARMA, col = "green")
# gumbel copula
plot(par4, type = "l", xlab = "", ylab = "", main = "gumbel copula",
     ylim = c(min(rho4) - .5, max(par4) + .05))
abline(h = fit4$par, col = "blue")
lines(rho4, col = "red")
lines(rho4_ARMA, col = "green")
# frank copula
plot(par5, type = "l", xlab = "", ylab = "", main = "frank copula", 
     ylim = c(min(rho5_ARMA) - .3, max(par5) + .05))
abline(h = fit5$par, col = "blue")
lines(rho5, col = "red")
lines(rho5_ARMA, col = "green")
# joe copula
plot(par6, type = "l", xlab = "", ylab = "", main = "joe copula",
     ylim = c(min(rho6_ARMA) - .1, max(par6) + .05))
abline(h = fit6$par, col = "blue")
lines(rho6, col = "red")
lines(rho6_ARMA, col = "green")

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,
       legend = c("实际参数", "静态参数", "GAS时变参数", "ARMA时变参数"), 
       col = c("black", "blue", "red", "green"), lwd = 1, cex = .8, horiz = TRUE)
