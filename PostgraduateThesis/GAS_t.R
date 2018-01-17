#t-copula pdf
tcopula_pdf<-function(rhot,data){
  nu <- nu0
  x1<-qt(data[,1],nu)
  x2<-qt(data[,2],nu)
  t<-length(data[,1])
  nu<-matrix(1,t,1)*nu
  rhot<-matrix(1,t,1)*rhot
  pdf<--999.99*matrix(1,t,1)
  for(i in 1:t){
    pdf[i]<-gamma((nu[i]+2)/2)/gamma(nu[i]/2)*((nu[i]*pi)^(-1))/sqrt(1-rhot[i]^2)*
      (1+(x1[i]^2+x2[i]^2-2*rhot[i]*x1[i]*x2[i])/(nu[i]*(1-rhot[i]^2)))^(-(nu[i]+2)/2)
    pdf[i]<-pdf[i]/(dt(x1[i],nu[i])*dt(x2[i],nu[i]))
  }
  return(pdf)
}

# t-copula的似然函数
tcopulaCL<-function(rhot,data){
  nu <- nu0
  if(nu>=100){
    x<-qnorm(data[,1])
    y<-qnorm(data[,2])
  }
  else{
    x<-qt(data[,1],nu)
    y<-qt(data[,2],nu)
  }
  CL<-lgamma((nu+2)/2)+lgamma(nu/2)-2*lgamma((nu+1)/2)-0.5*log(1-rhot^2)
  CL<-CL-(nu+2)/2*log(1+(x^2+y^2-2*rhot*x*y)/(nu*(1-rhot^2)))
  CL<-CL+(nu+1)/2*log(1+x^2/nu)+(nu+1)/2*log(1+y^2/nu)
  CL<-sum(CL)
  if(is.na(CL)){
    CL<--1e8
  }
  return(CL=-CL)
}

# t-copula的似然函数
tcopulaCLa<-function(rhot,data){
  nu <- nu0
  if(nu>=100){
    x<-qnorm(data[,1])
    y<-qnorm(data[,2])
  }
  else{
    x<-qt(data[,1],nu)
    y<-qt(data[,2],nu)
  }
  CL<-lgamma((nu+2)/2)+lgamma(nu/2)-2*lgamma((nu+1)/2)-0.5*log(1-rhot^2)
  CL<-CL-(nu+2)/2*log(1+(x^2+y^2-2*rhot*x*y)/(nu*(1-rhot^2)))
  CL<-CL+(nu+1)/2*log(1+x^2/nu)+(nu+1)/2*log(1+y^2/nu)
  return(CLa=-CL)
}

# 定义样条插值函数
LLgrad_t<-function(theta,data){
  l<-length(data[,1])
  B<-matrix(NA,nrow = l,ncol = 2)
  B[,1] <- tcopulaCLa(theta,data)
  h <- 2.2204e-16^(1/3)*max(abs(theta),0.01)#相当于导数定义的取微小值
  thetad <- theta - h
  B[,2] <- tcopulaCLa(thetad,data)
  B2 <- (B[,1]-B[,2])/h
  return(grad=B2)
}

# time-varing
tcopula_GAS_CL<-function(theta,data,hessINFO){
  RBAR<-0.9999
  t<-nrow(data)
  w<-theta[1]
  a<-theta[2]
  b<-theta[3]
  nu <- nu0
  h<-0.0001 #计算得分所用步长
  # 生成rhot时间序列
  ft<-rep(0,t)
  rhot<-rep(0,t)
  rhot[1]<-rho0
  ft[1]<-log((RBAR+rhot[1])/(RBAR-rhot[1]))
  for(i in 2:t){
    It<-ml_spline(hessINFO[,1],hessINFO[,2],rhot[i-1])
    DELTAt<-(-tcopulaCL(rhot[i-1]+h,t(data[i-1,]))--tcopulaCL(rhot[i-1],t(data[i-1,])))/h
    drhodf<-2*RBAR*exp(-ft[i-1])/((1+exp(-ft[i-1]))^2)
    Itildat<-It/(drhodf^2)
    DELTAtildat<-DELTAt/drhodf
    Itildat<-max(Itildat,1e-6)
    DELTAtildat<-max(min(DELTAtildat,1e4),-1e4)
    ft[i]<-w+b*ft[i-1]+a*DELTAtildat/sqrt(Itildat)
    ft[i]<-max(min(ft[i],100),-100)
    rhot[i]<-RBAR*(1-exp(-ft[i]))/(1+exp(-ft[i]))
  }
  return(rhot)
}

# 极大似然函数
tcopula_GAS_maxlik<-function(theta,data,hessINFO){
  rhot<-tcopula_GAS_CL(theta,data,hessINFO)
  CL<-tcopulaCL(rhot,data)
}

