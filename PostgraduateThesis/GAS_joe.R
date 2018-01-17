joe_pdf<-function(theta,data){
  ub<-1-data[,1]
  vb<-1-data[,2]
  pdf<-(ub^theta+vb^theta-(ub*vb)^theta)^(-2+1/theta)
  pdf<-pdf*ub^(theta-1)*vb^(theta-1)*(theta-1+ub^theta+vb^theta-(ub*vb)^theta)
  return(pdf)
}

joeCL<-function(theta,data){
  CL<-sum(log(joe_pdf(theta,data)))
  if(is.na(CL)){
    CL<--1e8
  }
  return(CL=-CL)
}

joeCLa<-function(theta,data){
  CLa<-log(joe_pdf(theta,data))
  return(CLa=-CLa)
}

# 样条差值函数
LLgrad_j<-function(theta,data){
  l<-length(data[,1])
  B<-matrix(NA,nrow = l,ncol = 2)
  B[,1] <- joeCLa(theta,data)
  h <- 2.2204e-16^(1/3)*max(abs(theta),0.01)#相当于导数定义的取微小值
  thetad <- theta - h
  B[,2] <- joeCLa(thetad,data)
  B2 <- (B[,1]-B[,2])/h
  return(grad=B2)
}

################################################
reps<-1e4 # 生成hessian矩阵
K<-matrix(seq(1.01,20,length.out =100),ncol = 1)
HESSjoe<-matrix(NA,nrow=length(K),ncol=1)
for(i in 1:length(K)){
  type<-joeCopula(K[i],dim=2)
  Usim<-rCopula(reps,type)
  scoreSIM<-LLgrad_j(K[i],Usim)#相当于求delta
  HESSjoe[i]<-mean(scoreSIM^2,na.rm = T) #相当于求delta平方
}
hessINFO6<-cbind(K,HESSjoe)
# 时变函数
joe_GAS_CL<-function(theta,data,hessINFO){
  t<-nrow(data)
  w<-theta[1]
  a<-theta[2]
  b<-theta[3]
  h<-0.0001 #计算得分所用步长
  # 生成rhot时间序列
  ft<-rep(0,t)
  rho<-rep(0,t)
  rho[1]<-rho0
  ft[1]<- log(rho[1]-1)
  for(i in 2:t){
    It<-ml_spline(hessINFO[,1],hessINFO[,2],rho[i-1])
    DELTAt<-(-joeCL(rho[i-1]+h,t(data[i-1,]))--joeCL(rho[i-1],t(data[i-1,])))/h
    drhodf<- exp(ft[i-1])
    Itildat<-It/(drhodf^2)
    DELTAtildat<-DELTAt/drhodf
    ft[i]<-w+b*ft[i-1]+a*DELTAtildat/sqrt(Itildat)
    rho[i]<-1+exp(ft[i])
  }
  return(rho)
}

# 极大似然函数
joe_GAS_maxlik<-function(theta,data,hessINFO){
  rho<-joe_GAS_CL(theta,data,hessINFO)
  CL<-joeCL(rho,data)
}