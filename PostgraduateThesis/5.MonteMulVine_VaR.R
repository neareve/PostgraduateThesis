library(copula)
source("VaR.R")
#############################################################
# Muti t-copula Estimation
(normal_m <- fitCopula(normalCopula(dim = 5), udata))
(t_m <- fitCopula(tCopula(dim = 5), udata)) # 选择多元t分布
tcop <- tCopula(param = t_m@copula@parameters[1], dim = 5, df = t_m@copula@parameters[2])
(clayton_m <- fitCopula(claytonCopula(dim = 5), udata))
(gumbel_m <- fitCopula(gumbelCopula(dim = 5), udata))
(frank_m <- fitCopula(frankCopula(dim = 5), udata))
(joe_m <- fitCopula(joeCopula(dim = 5), udata))

# Simulation
Monte_t_VaR <- function(start, end, len = 1000, alpha = .05, W = rep(.2, 5)) {
  Daily_VaR <- Daily_ES <- NA
  Rt <- matrix(NA, nrow = len, ncol = 5)
  for (t in start:end) {
    print(paste("Counting Day", t, " ", Sys.time(), sep = ""))
    x <- rCopula(len, tcop)
    # trace back log return
    for (i in 1:5) {
      Rt[, i] <- qdist("sstd", p = x[, i], mu = garch_m[t, i], sigma = garch_s[t, i],
                       skew = garch_par[1, i], shape = garch_par[2, i],
                       lambda = garch_par[3, i])
    }
    Daily_VaR[t - start + 1] <- risk(W, Rt, alpha)
  }
  return(Daily_VaR)
}

try <- Monte_t_VaR(1, l)
try2 <- Monte_t_VaR(1, l, alpha = .01)
try3 <- Monte_t_VaR(1, l, alpha = .1)
plot(day, return_seq, xlab = "",
     ylim = c(min(c(return_seq, try, try2, try3) - .02),
              max(c(return_seq + 0.02, try, try2, try3))),
     ylab = "", main = "静态多元Copula", cex = 0.7, xaxt = "n")
axis(1, c(day[1], day[500], day[1000]), c(day[1], day[500], day[1000]))
lines(day, return_seq)
lines(day, try, col = "blue") #5%
lines(day, try2, col = "red")
lines(day, try3, col = "green")
legend("bottomright", c("组合回报率", "1%VaR", "5%VaR", "10%VaR"), text.width = 200,
       y.intersp = .5,
       col = c("black", "red", "blue", "green"), lty = c(1, 1, 1, 1), pch = c(1, NA, NA, NA), cex = 0.5, inset = 0)

VaRTest(.05, return_seq, try)
VaRTest(.01, return_seq, try2)
VaRTest(.1, return_seq, try3)

##################################Vine############################################
# static Vine Copula Estimation
VineCop <- CDVine::CDVineCopSelect(udata, type = 2, familyset = c(1:6, 13, 14, 16), indeptest = T)
par <- VineCop$par
par2 <- VineCop$par2
family_S <- VineCop$family
l <- nrow(udata)
Monte_sVine_VaR <- function(start, end, len = 1000, alpha = .05, W = rep(.2, 5)) {
  Daily_VaR <- NA
  Rt <- matrix(NA, nrow = len, ncol = 5)
  for (t in start:end) {
    print(paste("Counting Day", t, " ", Sys.time(), sep = ""))
    x <- CDVine::CDVineSim(len, family = family_S, par = par, par2 = par2, type = 2)
    for (i in 1:5) {
      Rt[, i] <- qdist("sstd", p = x[, i], mu = garch_m[t, i], sigma = garch_s[t, i],
                       skew = garch_par[1, i], shape = garch_par[2, i],
                       lambda = garch_par[3, i])
    }
    # trace back log return
    Daily_VaR[t - start + 1] <- risk(W, Rt, alpha)
  }
  return(Daily_VaR)
}

return_seq <- .2 * apply(data_ret, 1, sum)
try4 <- Monte_sVine_VaR(1, l)
try5 <- Monte_sVine_VaR(1, l, alpha = .01)
try6 <- Monte_sVine_VaR(1, l, alpha = .1)

plot(day, return_seq, xlab = "",
     ylim = c(min(c(return_seq, try4, try5, try6)),
              max(c(return_seq + 0.02, try4, try5, try6))),
     ylab = "", main = "静态D藤Copula", cex = 0.7, xaxt = "n")
axis(1, c(day[1], day[500], day[1000]), c(day[1], day[500], day[1000]))
lines(day, return_seq)
lines(day, try4, col = "blue") #5%
lines(day, try5, col = "red")
lines(day, try6, col = "green")
legend("bottomright", c("组合回报率", "1%VaR", "5%VaR", "10%VaR"), text.width = 200,
       y.intersp = .5,
       col = c("black", "red", "blue", "green"), lty = c(1, 1, 1, 1), pch = c(1, NA, NA, NA), cex = 0.5, inset = 0)

VaRTest(.05, return_seq, try4)
VaRTest(.01, return_seq, try5)
VaRTest(.1, return_seq, try6)
