#################################dynamic portfolio######################################
source("One_Step.R")
rho_onestep <- matrix(NA, nrow = 4, ncol = 4)
rho_onestep[1, 1] <- onestep(theta_end[, 1, 1], rho_end[, 1, 1], udata[, 1:2], hessM[,, 1, 1], nu_end[1, 1])
rho_onestep[1, 2] <- onestep(theta_end[, 2, 1], rho_end[, 2, 1], udata[, 2:3], hessM[,, 2, 1], nu_end[1, 2])
rho_onestep[1, 3] <- t(rho_end[1, 3, 1])
rho_onestep[1, 4] <- onestep(theta_end[, 4, 1], rho_end[, 4, 1], udata[, 4:5], hessM[,, 4, 1], nu_end[1, 4])
rho_onestep[2:4,] <- t(rho_end[1,, 2:4])
#
opt_portfolio <- function(W = rep(.2, 5), alpha = .05, garch.obj, len = 1000,
                          fam, rho, nu = NULL) {
  Rt <- matrix(NA, nrow = len, ncol = 5)
  weight <- matrix(NA, nrow = 100, ncol = 5)
  obj <- mu <- sigma <- NA
  for (i in 1:5) {
    ahead <- ugarchforecast(garch.obj[[i]], n.ahead = 1)
    mu[i] <- fitted(ahead);
    sigma[i] <- sigma(ahead)
  }
  for (t in 1:100) {
    print(paste("Counting : ", t, sep = ""))
    w <- matrix(runif(len * 5), nrow = len, ncol = 5)
    v <- array(NA, c(len, 6, 6))
    x <- matrix(NA, nrow = len, ncol = 5)
    x[, 1] <- v[, 1, 1] <- w[, 1]
    x[, 2] <- v[, 1, 2] <- Hinv(fam[1, 1], w[, 2], v[, 1, 1],
                                rep(rho[1, 1], len), nu[1, 1]) # first could be i
    v[, 2, 2] <- Hfunc(fam[1, 1], v[, 1, 1], v[, 1, 2],
                       rep(rho[1, 1], len), nu[1, 1])
    for (i in 3:5) {
      v[, 1, i] <- w[, i]
      for (k in (i - 1):2) {
        v[, 1, i] <- Hinv(fam[k, i - k], v[, 1, i], v[, 2 * k - 2, i - 1],
                          rep(rho[k, i - k], len), nu[k, i - k])
      }
      v[, 1, i] <- Hinv(fam[1, i - 1], v[, 1, i], v[, 1, i - 1],
                        rep(rho[1, i - 1], len), nu[1, i - 1])
      x[, i] <- v[, 1, i]
      if (i == 5) { break }
      v[, 2, i] <- Hfunc(fam[1, i - 1], v[, 1, i - 1], v[, 1, i],
                         rep(rho[1, i - 1], len), nu[1, i - 1])
      v[, 3, i] <- Hfunc(fam[1, i - 1], v[, 1, i], v[, 1, i - 1],
                         rep(rho[1, i - 1], len), nu[1, i - 1])
      if (i > 3) {
        for (j in 2:(i - 2)) {
          v[, 2 * j, i] <- Hfunc(fam[j, i - j], v[, 2 * j - 2, i - 1], v[, 2 * j - 1, i],
                                 rep(rho[j, i - j], len), nu[j, i - j])
          v[, 2 * j + 1, i] <- Hfunc(fam[j, i - j], v[, 2 * j - 1, i], v[, 2 * j - 2, i - 1],
                                     rep(rho[j, i - j], len), nu[j, i - j])
        }
      }
      v[, 2 * i - 2, i] <- Hfunc(fam[i - 1, 1], v[, 2 * i - 4, i - 1], v[, 2 * i - 3, i],
                                 rep(rho[i - 1, 1], len), nu[i - 1, 1])
    }
    # 返回对数收益率
    for (i in 1:5) {
      Rt[, i] <- qdist("sstd", p = x[, i], mu = mu[i], sigma = sigma[i],
                       skew = garch_par[1, i], shape = garch_par[2, i],
                       lambda = garch_par[3, i])
    }
    # couting VaR
    nl <- nloptr(W, toOpt, lb = rep(0, 5), ub = rep(1, 5), eval_g_eq = eqCon,
                 opts = opts, data = Rt, alpha = alpha)
    weight[t,] <- nl$solution
    obj[t] <- nl$objective
  }
  opt_weight <- apply(weight, 2, mean)
  opt_obj <- mean(obj)
  return(list(opt_weight, opt_obj))
}
port_GAS <- opt_portfolio(garch.obj = fit, fam = fam_end, rho = rho_onestep, nu = nu_end)
port_GAS_l <- opt_portfolio(garch.obj = fit, fam = fam_end, rho = rho_onestep, nu = nu_end, alpha = .01)
port_GAS_u <- opt_portfolio(garch.obj = fit, fam = fam_end, rho = rho_onestep, nu = nu_end, alpha = .1)
####################################static vine portfolio###################################
opt_portfolio_s <- function(W = rep(.2, 5), alpha = .05, garch.obj, len = 1000,
                            par, par2, family) {
  Rt <- matrix(NA, nrow = len, ncol = 5)
  weight <- matrix(NA, nrow = 100, ncol = 5)
  obj <- mu <- sigma <- NA
  for (i in 1:5) {
    ahead <- ugarchforecast(garch.obj[[i]], n.ahead = 1)
    mu[i] <- fitted(ahead);
    sigma[i] <- sigma(ahead)
  }
  for (t in 1:100) {
    print(paste("Counting : ", t, sep = ""))
    x <- CDVine::CDVineSim(len, family = family, par = par, par2 = par2, type = 2)
    for (i in 1:5) {
      Rt[, i] <- qdist("sstd", p = x[, i], mu = mu, sigma = sigma,
                       skew = garch_par[1, i], shape = garch_par[2, i],
                       lambda = garch_par[3, i])
    }
    nl <- nloptr(W, toOpt, lb = rep(0, 5), ub = rep(1, 5), eval_g_eq = eqCon,
                 opts = opts, data = Rt, alpha = alpha)
    weight[t,] <- nl$solution
    obj[t] <- nl$objective
  }
  opt_weight <- apply(weight, 2, mean)
  opt_obj <- mean(obj)
  return(list(opt_weight, opt_obj))
}
port_sVine <- opt_portfolio_s(garch.obj = fit, par = par, par2 = par2, family = family_S)
port_sVine_l <- opt_portfolio_s(garch.obj = fit, par = par, par2 = par2, family = family_S, alpha = .01)
port_sVine_u <- opt_portfolio_s(garch.obj = fit, par = par, par2 = par2, family = family_S, alpha = .1)
############################################multiple portfolio#############################

opt_portfolio_t <- function(len = 1000, alpha = .05, W = rep(.2, 5), cop, garch.obj) {
  Rt <- matrix(NA, nrow = len, ncol = 5)
  weight <- matrix(NA, nrow = 100, ncol = 5)
  obj <- mu <- sigma <- NA
  for (i in 1:5) {
    ahead <- ugarchforecast(garch.obj[[i]], n.ahead = 1)
    mu[i] <- fitted(ahead);
    sigma[i] <- sigma(ahead)
  }
  for (t in 1:100) {
    print(paste("Counting : ", t, sep = ""))
    x <- rCopula(len, cop)
    # trace back log return
    for (i in 1:5) {
      Rt[, i] <- qdist("sstd", p = x[, i], mu = mu, sigma = sigma,
                       skew = garch_par[1, i], shape = garch_par[2, i],
                       lambda = garch_par[3, i])
    }
    nl <- nloptr(W, toOpt, lb = rep(0, 5), ub = rep(1, 5), eval_g_eq = eqCon,
                 opts = opts, data = Rt, alpha = alpha)
    weight[t,] <- nl$solution
    obj[t] <- nl$objective
  }
  opt_weight <- apply(weight, 2, mean)
  opt_obj <- mean(obj)
  return(list(opt_weight, opt_obj))
}
port_t <- opt_portfolio_t(cop = tcop, garch.obj = fit)
port_t_l <- opt_portfolio_t(cop = tcop, garch.obj = fit, alpha = .01)
port_t_u <- opt_portfolio_t(cop = tcop, garch.obj = fit, alpha = .1)
