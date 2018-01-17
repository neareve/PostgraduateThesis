onestep <- function(theta,rho,data,hessINFO,nu=NULL){
  # theta = [w,a,b,nuinv], the parameters of the GAS specification, and the *inverse* DoF parameter
  # lag_rho, a scalar, the lagged rho
  # lag_U, a 1x2 vector, the lagged value of [U1,U2]
  # rho0, a scalar, the value of the correlation parameter to use as the starting value (perhaps the estimate from a constant version of this model)
  # RR, a k1x1 vector of values of rho at which the hessian was computed
  # NN, a k2x1 vector of values of nu at which the hessian was computed
  # HESSstudt, a k1 x k2 x 2 x 2 matrix, containing the 2x2 hessian for each combination of values of [rho,nu]
  w<-theta[1]
  a<-theta[2]
  b<-theta[3]
  h<-1e-4
  lag_rho<-rho
  RBAR <- .9999
  lag_f <- log((RBAR+lag_rho)/(RBAR-lag_rho))
  It<-approx(hessINFO[,1],hessINFO[,2],lag_rho)$y
  DELTAt<-(-tcopulaCL(lag_rho+h,data)--tcopulaCL(lag_rho,data))/h
  drhodf<-2*RBAR*exp(-lag_f)/((1+exp(-lag_f))^2)
  Itildat<-It/(drhodf^2)
  DELTAtildat <- DELTAt/drhodf
  Itildat <- max(Itildat,1e-6)
  DELTAtildat <- max(min(DELTAtildat,1e4),-1e4)
  ft <- w + b*lag_f + a*DELTAtildat/sqrt(Itildat)
  ft <- max(min(ft,100),-100)
  rho_new <- RBAR*(1-exp(-ft))/(1+exp(-ft))
  return(rho_new)
}

