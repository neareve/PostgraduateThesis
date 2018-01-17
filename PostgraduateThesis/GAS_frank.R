# para definition (0,Inf)
frank_pdf<-function(theta,data){
  u<-data[,1]
  v<-data[,2]
  pdf<--theta*(exp(-theta)-1)*exp(-theta*(u+v))
  pdf<-pdf/((exp(-theta)-1)+(exp(-theta*u)-1)*(exp(-theta*v)-1))^2
  return(pdf)
}

frankCL<-function(theta,data){
  CL<-sum(log(frank_pdf(theta,data)))
  if(is.na(CL)){
    CL<--1e8
  }
  return(CL=-CL)
}

frankCLa<-function(theta,data){
  CLa<-log(frank_pdf(theta,data))
  return(CLa=-CLa)
}

# ������ֵ����
LLgrad_f<-function(theta,data){
  l<-length(data[,1])
  B<-matrix(NA,nrow = l,ncol = 2)
  B[,1] <- frankCLa(theta,data)
  h <- 2.2204e-16^(1/3)*max(abs(theta),0.01)#�൱�ڵ��������ȡ΢Сֵ
  thetad <- theta - h
  B[,2] <- frankCLa(thetad,data)
  B2 <- (B[,1]-B[,2])/h
  return(grad=B2)
}

######################################################
reps<-1e4 # ����hessian����
K<-matrix(c(seq(-10,-.01,length.out = 50),seq(.01,10,length.out = 50)),ncol = 1)
HESSfrank<-matrix(NA,nrow=length(K),ncol=1)
for(i in 1:length(K)){
  type<-frankCopula(K[i],dim=2)
  Usim<-rCopula(reps,type)
  scoreSIM<-LLgrad_f(K[i],Usim)#�൱����delta
  HESSfrank[i]<-mean(scoreSIM^2,na.rm = T) #�൱����deltaƽ��
}
hessINFO5<-cbind(K,HESSfrank)
# ʱ�亯��
frank_GAS_CL<-function(theta,data,hessINFO){
  t<-nrow(data)
  w<-theta[1]
  a<-theta[2]
  b<-theta[3]
  h<-.0001 #����÷����ò���
  # ����rhotʱ������
  ft<-rep(0,t)
  rho<-rep(0,t)
  rho[1]<-rho0
  ft[1] <- 1/rho[1]-1
  for(i in 2:t){
    It<-ml_spline(hessINFO[,1],hessINFO[,2],rho[i-1])
    DELTAt<-(-frankCL(rho[i-1]+h,t(data[i-1,]))--frankCL(rho[i-1],t(data[i-1,])))/h
    drhodf<- -1/(ft[i] + 1)^2
    Itildat<-It/(drhodf^2)
    DELTAtildat<-DELTAt/drhodf
    ft[i]<-w+b*ft[i-1]+a*DELTAtildat/sqrt(Itildat)
    rho[i]<-1/(ft[i] + 1)
  }
  return(rho)
}
# ������Ȼ����
frank_GAS_maxlik<-function(theta,data,hessINFO){
  rho<-frank_GAS_CL(theta,data,hessINFO)
  CL<-frankCL(rho,data)
}