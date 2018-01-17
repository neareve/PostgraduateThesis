library(copula)
library(VineCopula)
library(nloptr)
library(xts)
library(Rcpp)
source("GAS_clayton.R") #3
source("GAS_frank.R") #5
source("GAS_gumbel.R") #4
source("GAS_joe.R") #6
source("GAS_normal.R") #1
source("GAS_t.R") #2
source("GAS_rotclayton.R") #13
source("GAS_rotgumbel.R") #14
source("GAS_rotjoe.R") #16
source("IndepTest.R")
sourceCpp("hinv.cpp") #包含独立性检验，原假设为独立copula
sourceCpp("MlSpline.cpp")
#############################################################

# 静态选择
CDVine::CDVineCopSelect(udata, familyset = c(1:6, 13, 14, 16),
                        type = 2, indeptest = T)

#cop <- tCopula(.3,dim = 5)
#udata <- rCopula(100,cop)
# initialize
V <- udata
d <- ncol(V)
l <- nrow(V)
family <- c(1:6, 13, 14, 16)
theta_end <- array(NA, c(3, (d - 1), (d - 1)))
rho_end <- array(NA, c(l, (d - 1), (d - 1)))
fam_end <- matrix(NA, nrow = d - 1, ncol = d - 1)
nu_end <- matrix(NA, nrow = d - 1, ncol = d - 1)
theta <- matrix(0, nrow = 3, ncol = 9)
rho <- matrix(0, nrow = l, ncol = 9)
hessM <- array(0, c(100, 2, 4, d - 1))
aic <- 0
aic_s <- 0

# select optimum pair copula
for (i in 1:(d - 1)) {
  for (j in 1:(d - i)) {
    print(paste("Level:", i, ",", "Pair:", j, "  ", Sys.time(), sep = ""))
    data <- cbind(V[, j], V[, j + 1])
    p <- indep_test(data[, 1], data[, 2])
    if (p > .05) {
      theta_end[, j, i] <- 0
      rho_end[, j, i] <- 0
      fam_end[i, j] <- 0
    } else {
      # 1.normal copula
      print(paste("1:", "normalcopula", " ", Sys.time()))
      bicop <- BiCopSelect(data[, 1], data[, 2], familyset = 1)
      rho0 <- bicop$par
      aic_s[1] <- bicop$AIC
      w0 <- -log((1 - rho0) / (1 + rho0))
      out <- nloptr(x0 = c(w0, .1, .1), eval_f = ncopula_GAS_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
                    opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval = 10000),
                    data = data, hessINFO = hessINFO1)
      theta[, 1] <- out$solution
      rho[, 1] <- ncopula_GAS_CL(theta[, 1], data, hessINFO1)
      aic[1] <- 2 * 3 - 2 * (-ncopula_GAS_maxlik(theta[, 1], data, hessINFO1))
      if (aic_s[1] < aic[1]) {
        theta[, 1] <- c(w0, 0, 0)
        rho[, 1] <- rep(rho0, l)
        aic[1] <- aic_s[1]
      }
      # 2.tcopula
      print(paste("2:", "tcopula", " ", Sys.time()))
      bicop <- BiCopSelect(data[, 1], data[, 2], familyset = 2)
      nu_end[i, j] <- nu0 <- bicop$par2
      rho0 <- bicop$par
      aic_s[2] <- bicop$AIC
      w0 <- -log((1 - rho0) / (1 + rho0))
      # 生成hessen矩阵
      reps <- 1e4
      R <- matrix(c(seq(-.99, -.5, length = 40), seq(-.49, .5, length = 20), seq(.51, .99, length = 40)), ncol = 1)
      HESSstudt <- matrix(NA, length(R), ncol = 1)
      for (k in 1:length(R)) {
        type <- tCopula(R[k], df = nu0, dim = 2)
        Usim <- rCopula(reps, type)
        scoreSIM <- LLgrad_t(R[k], Usim)
        HESSstudt[k] <- mean(scoreSIM ^ 2, na.rm = T)
      }
      hessM[,, j, i] <- hessINFO2 <- cbind(R, HESSstudt[, 1]) #唯一确定t的随机矩阵位置
      out <- nloptr(x0 = c(w0, .1, .1), eval_f = tcopula_GAS_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
                    opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval = 10000),
                    data = data, hessINFO = hessINFO2)
      theta[, 2] <- out$solution
      rho[, 2] <- tcopula_GAS_CL(theta[, 2], data, hessINFO2)
      aic[2] <- 2 * 3 - 2 * (-tcopula_GAS_maxlik(theta[, 2], data, hessINFO2))
      if (aic_s[2] < aic[2]) {
        theta[, 2] <- c(w0, 0, 0)
        rho[, 2] <- rep(rho0, l)
        aic[2] <- aic_s[2]
      }
      # 3.clayton copula
      print(paste("3:", "claytoncopula", " ", Sys.time()))
      bicop <- BiCopSelect(data[, 1], data[, 2], familyset = 3)
      rho0 <- bicop$par
      if (rho0 > -1) {
        aic_s[3] <- bicop$AIC
        w0 <- log(1 + rho0)
        out <- nloptr(x0 = c(w0, .1, .1), eval_f = clayton_GAS_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
                      opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval = 10000),
                      data = data, hessINFO = hessINFO3)
        theta[, 3] <- out$solution
        rho[, 3] <- clayton_GAS_CL(theta[, 3], data, hessINFO3)
        aic[3] <- 2 * 3 - 2 * (-clayton_GAS_maxlik(theta[, 3], data, hessINFO3))
        if (aic_s[3] < aic[3]) {
          theta[, 3] <- c(w0, 0, 0)
          rho[, 3] <- rep(rho0, l)
          aic[3] <- aic_s[3]
        }
      } else {
        theta[, 3] <- 0
        rho[, 3] <- 0
        aic[3] <- 1e8
      }
      # 4.gumbel copula
      print(paste("4:", "gumbelcopula", " ", Sys.time()))
      bicop <- BiCopSelect(data[, 1], data[, 2], familyset = 4)
      rho0 <- bicop$par
      if (rho0 > 1) {
        aic_s[4] <- bicop$AIC
        w0 <- log(rho0 - 1)
        out <- nloptr(x0 = c(w0, .1, .1), eval_f = gumbel_GAS_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
                      opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval = 10000),
                      data = data, hessINFO = hessINFO4)
        theta[, 4] <- out$solution
        rho[, 4] <- gumbel_GAS_CL(theta[, 4], data, hessINFO4)
        aic[4] <- 2 * 3 - 2 * (-gumbel_GAS_maxlik(theta[, 4], data, hessINFO4))
        if (aic_s[4] < aic[4]) {
          theta[, 4] <- c(w0, 0, 0)
          rho[, 4] <- rep(rho0, l)
          aic[4] <- aic_s[4]
        }
      } else {
        theta[, 4] <- 0
        rho[, 4] <- 0
        aic[4] <- 1e8
      }
      # 5.frank copula
      print(paste("5:", "frankcopula", " ", Sys.time()))
      bicop <- BiCopSelect(data[, 1], data[, 2], familyset = 5)
      rho0 <- bicop$par
      aic_s[5] <- bicop$AIC
      w0 <- 1 / rho0 - 1
      out <- nloptr(x0 = c(w0, .1, .1), eval_f = frank_GAS_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
                    opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval = 10000),
                    data = data, hessINFO = hessINFO5)
      theta[, 5] <- out$solution
      rho[, 5] <- frank_GAS_CL(theta[, 5], data, hessINFO5)
      aic[5] <- 2 * 3 - 2 * (-frank_GAS_maxlik(theta[, 5], data, hessINFO5))
      if (aic_s[5] < aic[5]) {
        theta[, 5] <- c(w0, 0, 0)
        rho[, 5] <- rep(rho0, l)
        aic[5] <- aic_s[5]
      }
      # 6.joe copula
      print(paste("6:", "joecopula", " ", Sys.time()))
      bicop <- BiCopSelect(data[, 1], data[, 2], familyset = 6)
      rho0 <- bicop$par
      if (rho0 > 1) {
        aic_s[6] <- bicop$AIC
        w0 <- log(rho0 - 1)
        out <- nloptr(x0 = c(w0, .1, .1), eval_f = joe_GAS_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
                      opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval = 10000),
                      data = data, hessINFO = hessINFO6)
        theta[, 6] <- out$solution
        rho[, 6] <- joe_GAS_CL(theta[, 6], data, hessINFO6)
        aic[6] <- 2 - 2 * (-joe_GAS_maxlik(theta[, 6], data, hessINFO6))
        if (aic_s[6] < aic[6]) {
          theta[, 6] <- c(w0, 0, 0)
          rho[, 6] <- rep(rho0, l)
          aic[6] <- aic_s[6]
        }
      } else {
        theta[, 6] <- 0
        rho[, 6] <- 0
        aic[6] <- 1e8
      }
      # 13 rotclayton copula
      print(paste("7:", "rotclaytoncopula", " ", Sys.time()))
      bicop <- BiCopSelect(data[, 1], data[, 2], familyset = 13)
      rho0 <- bicop$par
      if (rho0 > -1) {
        aic_s[7] <- bicop$AIC
        w0 <- log(1 + rho0)
        out <- nloptr(x0 = c(w0, .1, .1), eval_f = rotclayton_GAS_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
                      opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval = 10000),
                      data = data, hessINFO = hessINFO13)
        theta[, 7] <- out$solution
        rho[, 7] <- rotclayton_GAS_CL(theta[, 7], data, hessINFO13)
        aic[7] <- 2 * 3 - 2 * (-rotclayton_GAS_maxlik(theta[, 7], data, hessINFO13))
        if (aic_s[7] < aic[7]) {
          theta[, 7] <- c(w0, 0, 0)
          rho[, 7] <- rep(rho0, l)
          aic[7] <- aic_s[7]
        }
      } else {
        theta[, 7] <- 0
        rho[, 7] <- 0
        aic[7] <- 1e8
      }
      #14 rotgumbel copula
      print(paste("8:", "rotgumbelcopula", " ", Sys.time()))
      bicop <- BiCopSelect(data[, 1], data[, 2], familyset = 14)
      rho0 <- bicop$par
      if (rho0 > 1) {
        aic_s[8] <- bicop$AIC
        w0 <- log(rho0 - 1)
        out <- nloptr(x0 = c(w0, .1, .1), eval_f = rotgumbel_GAS_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
                      opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval = 10000),
                      data = data, hessINFO = hessINFO14)
        theta[, 8] <- out$solution
        rho[, 8] <- rotgumbel_GAS_CL(theta[, 8], data, hessINFO14)
        aic[8] <- 2 * 3 - 2 * (-rotgumbel_GAS_maxlik(theta[, 8], data, hessINFO14))
        if (aic_s[8] < aic[8]) {
          theta[, 8] <- c(w0, 0, 0)
          rho[, 8] <- rep(rho0, l)
          aic[8] <- aic_s[8]
        }
      } else {
        theta[, 8] <- 0
        rho[, 8] <- 0
        aic[8] <- 1e8
      }
      # 16. rotjoe copula
      print(paste("9:", "rotjoecopula", " ", Sys.time()))
      bicop <- BiCopSelect(data[, 1], data[, 2], familyset = 16)
      rho0 <- bicop$par
      if (rho0 > 1) {
        aic_s[9] <- bicop$AIC
        w0 <- log(rho0 - 1)
        out <- nloptr(x0 = c(w0, .1, .1), eval_f = rotjoe_GAS_maxlik, lb = rep(-10, 3), ub = rep(10, 3),
                      opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval = 10000),
                      data = data, hessINFO = hessINFO16)
        theta[, 9] <- out$solution
        rho[, 9] <- rotjoe_GAS_CL(theta[, 9], data, hessINFO16)
        aic[9] <- 2 * 3 - 2 * (-rotjoe_GAS_maxlik(theta[, 9], data, hessINFO16))
        if (aic_s[9] < aic[9]) {
          theta[, 9] <- c(w0, 0, 0)
          rho[, 9] <- rep(rho0, l)
          aic[9] <- aic_s[9]
        }
      } else {
        theta[, 9] <- 0
        rho[, 9] <- 0
        aic[9] <- 1e8
      }
      # 选择BIC最小时变参数
      n <- which.min(aic)
      theta_end[, j, i] <- theta[, n]
      rho_end[, j, i] <- rho[, n]
      fam_end[i, j] <- family[n]
    }
  }
  temp <- matrix(0, nrow = l, ncol = (d - i))
  for (k in 1:(d - i)) {
    temp[, k] <- Hfunc(fam_end[i, k], V[, k], V[, k + 1], rho_end[, k, i], nu_end[i, k])
  }
  V <- temp
}
##########################################plot###############################
day <- as.Date(mdata$index[-1])
tree1 <- rho_end[,, 1]
tree1 <- xts(tree1, order.by = day)
par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))
plot(tree1[, 1], main = "BVSP-SP500", yaxis.right = F)
plot(tree1[, 2], main = "SP500-FT100", yaxis.right = F)
plot(tree1[, 4], main = "N225-SSEC", yaxis.right = F)
