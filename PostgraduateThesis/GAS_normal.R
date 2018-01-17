# normal copula pdf
ncopula_pdf<-function(rho,data){
  x<-qnorm(data[,1])
  y<-qnorm(data[,2])
  pdf<-(1-rho^2)^(-0.5)*exp(-0.5*(1-rho^2)^(-1)*(x^2+y^2-2*rho*x*y))*exp(0.5*(x^2+y^2))
  return(pdf)
}

# ncopula 似然函数
ncopulaCL<-function(theta,data){
  x<-qnorm(data[,1])
  y<-qnorm(data[,2])
  out1<-(1-theta^2)^(-0.5)
  out2<-exp(-0.5*(1-theta^2)^(-1)*(x^2+y^2-2*theta*x*y))
  out3<-exp(0.5*(x^2+y^2))
  out<-out1*out2*out3
  CL<-sum(log(out))
  if(is.na(CL)){
    CL<--1e8
  }
  return(CL=-CL)
}

# ncopula CLa
ncopulaCLa<-function(theta,data){
  x<-qnorm(data[,1])
  y<-qnorm(data[,2])
  CLa<--1*(2*(1-theta^2))^(-1)*(x^2+y^2-2*theta*x*y)
  CLa<-CLa+0.5*(x^2+y^2)
  CLa<-CLa-0.5*log(1-theta^2)
  return(CLa=-CLa)
}

# 样条差值函数
LLgrad_n<-function(theta,data){
  l<-length(data[,1])
  B<-matrix(NA,nrow = l,ncol = 2)
  B[,1] <- ncopulaCLa(theta,data)
  h <- 2.2204e-16^(1/3)*max(abs(theta),0.01)#相当于导数定义的取微小值
  thetad <- theta - h
  B[,2] <- ncopulaCLa(thetad,data)
  B2 <- (B[,1]-B[,2])/h
  return(grad=B2)
}

# Hessian matrix
reps<-1e4 # 生成hessian矩阵
K<-matrix(c(seq(-.99,-0.5,length=40),seq(-0.51,0.5,length=20),seq(0.51,.99,length=40)),ncol = 1)
HESSnormal<-matrix(NA,nrow=length(K),ncol=1)
for(i in 1:length(K)){
  type<-normalCopula(K[i],dim=2)
  Usim<-rCopula(reps,type)
  scoreSIM<-LLgrad_n(K[i],Usim)#相当于求delta
  HESSnormal[i]<-mean(scoreSIM^2,na.rm = T) #相当于求delta平方
}
hessINFO1<-cbind(K,HESSnormal)

# 时变函数
ncopula_GAS_CL<-function(theta,data,hessINFO){
  RBAR<-0.9999
  t<-nrow(data)
  w<-theta[1]
  a<-theta[2]
  b<-theta[3]
  h<-0.0001 #计算得分所用步长
  # 生成rhot时间序列
  ft<-rep(0,t)
  rho<-rep(0,t)
  rho[1]<-rho0
  ft[1]<-log((1+rho[1])/(1-rho[1]))
  for(i in 2:t){
    It<-ml_spline(hessINFO[,1],hessINFO[,2],rho[i-1])
    DELTAt<-(-ncopulaCL(rho[i-1]+h,t(data[i-1,]))--ncopulaCL(rho[i-1],t(data[i-1,])))/h
    drhodf<-2*exp(-ft[i-1])/((1+exp(-ft[i-1]))^2)
    Itildat<-It/(drhodf^2)
    DELTAtildat<-DELTAt/drhodf
    Itildat<-max(Itildat,1e-6)
    DELTAtildat<-max(min(DELTAtildat,1e4),-1e4)
    ft[i]<-w+b*ft[i-1]+a*DELTAtildat/sqrt(Itildat)
    rho[i]<-(1-exp(-ft[i]))/(1+exp(-ft[i]))
  }
  return(rho)
}

# 极大似然函数
ncopula_GAS_maxlik<-function(theta,data,hessINFO){
  rho<-ncopula_GAS_CL(theta,data,hessINFO)
  CL<-ncopulaCL(rho,data)
}
