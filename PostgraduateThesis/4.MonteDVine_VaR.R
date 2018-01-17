source("VaR.R")
######################################################################
# calculate daily VaR
Monte_VaR <- function(rho, fam, nu = NULL, start, end, alpha = .05, len = 1000,
                      W = rep(.2, 5)) {
  Daily_VaR <- NA
  Rt <- matrix(NA, nrow = len, ncol = 5)
  for (t in start:end) {
    print(paste("Counting Day", t, " ", Sys.time(), sep = ""))
    w <- matrix(runif(len * 5), nrow = len, ncol = 5)
    v <- array(NA, c(len, 6, 6))
    x <- matrix(NA, nrow = len, ncol = 5)
    x[, 1] <- v[, 1, 1] <- w[, 1]
    x[, 2] <- v[, 1, 2] <- Hinv(fam[1, 1], w[, 2], v[, 1, 1],
                                rep(rho[t, 1, 1], len), nu[1, 1]) # first could be i
    v[, 2, 2] <- Hfunc(fam[1, 1], v[, 1, 1], v[, 1, 2],
                       rep(rho[t, 1, 1], len), nu[1, 1])
    for (i in 3:5) {
      v[, 1, i] <- w[, i]
      for (k in (i - 1):2) {
        v[, 1, i] <- Hinv(fam[k, i - k], v[, 1, i], v[, 2 * k - 2, i - 1],
                          rep(rho[t, i - k, k], len), nu[k, i - k])
      }
      v[, 1, i] <- Hinv(fam[1, i - 1], v[, 1, i], v[, 1, i - 1],
                        rep(rho[t, i - 1, 1], len), nu[1, i - 1])
      x[, i] <- v[, 1, i]
      if (i == 5) { break }
      v[, 2, i] <- Hfunc(fam[1, i - 1], v[, 1, i - 1], v[, 1, i],
                         rep(rho[t, i - 1, 1], len), nu[1, i - 1])
      v[, 3, i] <- Hfunc(fam[1, i - 1], v[, 1, i], v[, 1, i - 1],
                         rep(rho[t, i - 1, 1], len), nu[1, i - 1])
      if (i > 3) {
        for (j in 2:(i - 2)) {
          v[, 2 * j, i] <- Hfunc(fam[j, i - j], v[, 2 * j - 2, i - 1], v[, 2 * j - 1, i],
                                 rep(rho[t, i - j, j], len), nu[j, i - j])
          v[, 2 * j + 1, i] <- Hfunc(fam[j, i - j], v[, 2 * j - 1, i], v[, 2 * j - 2, i - 1],
                                     rep(rho[t, i - j, j], len), nu[j, i - j])
        }
      }
      v[, 2 * i - 2, i] <- Hfunc(fam[i - 1, 1], v[, 2 * i - 4, i - 1], v[, 2 * i - 3, i],
                                 rep(rho[t, 1, i - 1], len), nu[i - 1, 1])
    }
    # 返回对数收益率
    for (i in 1:5) {
      Rt[, i] <- qdist("sstd", p = x[, i], mu = garch_m[t, i], sigma = garch_s[t, i],
                       skew = garch_par[1, i], shape = garch_par[2, i],
                       lambda = garch_par[3, i])
    }
    # couting VaR
    Daily_VaR[t - start + 1] <- risk(W, Rt, alpha)
  }
  return(Daily_VaR)
}
######################################BackTest##########################################
# backtest
return_seq <- .2 * apply(data_ret, 1, sum) # use weights equal to 1/5 to calculate returns
GAS_VaRm <- Monte_VaR(rho_end, fam_end, nu_end, start = 1, end = l, alpha = .05) # alpha=.05
GAS_VaRl <- Monte_VaR(rho_end, fam_end, nu_end, start = 1, end = l, alpha = .01) # alpha=.01
GAS_VaRu <- Monte_VaR(rho_end, fam_end, nu_end, start = 1, end = l, alpha = .1) # alpha=.1
# fail test
VaRTest(.05, return_seq, GAS_VaRm)
VaRTest(.01, return_seq, GAS_VaRl)
VaRTest(.1, return_seq, GAS_VaRu)

# 全区间
plot(day, return_seq, cex = .7, ylim = c(min(return_seq) - .02, max(return_seq) + .02),
     xlab = "", ylab = "", main = "动态D藤copula", xaxt = "n")
axis(1, c(day[1], day[500], day[1000]), c(day[1], day[500], day[1000]))
lines(day, return_seq)
lines(day, GAS_VaRm, col = "blue")
lines(day, GAS_VaRl, col = "red")
lines(day, GAS_VaRu, col = "green")
legend("bottomright", c("组合回报率", "1%VaR", "5%VaR", "10%VaR"), text.width = 200,
       y.intersp = .5,
       col = c("black", "red", "blue", "green"), lty = c(1, 1, 1, 1), pch = c(1, NA, NA, NA), cex = 0.5, inset = 0)
