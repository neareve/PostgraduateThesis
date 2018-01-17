gumbel_ARMA_CL<-function(theta,data){
  t <- nrow(data)
  w <- theta[1]
  a <- theta[2]
  b <- theta[3]
  u <- data[,1]
  v <- data[,2]
  ft <- rep(NA,t)
  rho <- rep(NA,t) 
  rho[1] <- rho0
  ft[1]<- log(rho[1]-1)
  for (i in 2:t){
    if(i <= 10){
      ft[i] <- w + a * ft[i-1] + b * mean(abs(u[1:i-1] - v[1:i-1]))
    }
    else{
      ft[i] <- w + a * ft[i-1] + b * mean(abs(u[(i-10):i-1] - v[(i-10):i-1]))
    }
    rho[i] <- 1 + exp(ft[i])
  }
  return(rho)
} 

gumbel_ARMA_maxlik <- function(theta, data){
  rho <- gumbel_ARMA_CL(theta, data)
  CL <- gumbelCL(rho, data)
}