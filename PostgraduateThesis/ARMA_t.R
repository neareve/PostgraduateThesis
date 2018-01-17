tcopula_ARMA_CL<-function(theta,data){
  RBAR <- .9999
  t<-nrow(data)
  u <- data[,1]
  v <- data[,2]
  w<-theta[1]
  a<-theta[2]
  b<-theta[3]
  nu <- nu0
  ft<-rep(0,t)
  rhot<-rep(0,t)
  rhot[1]<-rho0
  ft[1]<-log((RBAR+rhot[1])/(RBAR-rhot[1]))
  for(i in 2:t){
    if(i <= 10){
      ft[i] <- w + a * ft[i-1] + b * mean(abs(u[1:i-1] - v[1:i-1]))
    }else{
      ft[i] <- w + a * ft[i-1] + b * mean(abs(u[(i-10):i-1] - v[(i-10):i-1]))
    }
    rhot[i]<-RBAR*(1-exp(-ft[i]))/(1+exp(-ft[i]))
  }
  return(rhot)
}

# 极大似然函数
tcopula_ARMA_maxlik<-function(theta,data){
  rhot<-tcopula_ARMA_CL(theta,data)
  CL<-tcopulaCL(rhot,data)
}