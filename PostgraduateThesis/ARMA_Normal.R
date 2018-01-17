ncopula_ARMA_CL<-function(theta,data){
  t <- nrow(data)
  u <- data[,1]
  v <- data[,2]
  w <- theta[1]
  a <- theta[2]
  b <- theta[3]
  ft <- rep(0, t)
  rho <- rep(0, t)
  rho[1] <- rho0
  ft[1] <- log((1 + rho[1])/(1 - rho[1]))
  for(i in 2:t){
    if(i <= 10){
    ft[i] <- w + a * ft[i-1] + b * mean(abs(u[1:i-1] - v[1:i-1]))
  }else{
    ft[i] <- w + a * ft[i-1] + b * mean(abs(u[(i-10):i-1] - v[(i-10):i-1]))
  }
  rho[i] <- (1 - exp(-ft[i]))/(1 + exp(-ft[i]))
  }
  return(rho)
}

ncopula_ARMA_maxlik <- function(theta, data){
  rho <- ncopula_ARMA_CL(theta, data)
  CL <- ncopulaCL(rho, data)
}
