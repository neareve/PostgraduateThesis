reps<-1e4 # ����hessian����
K <- matrix(c(seq(-.99, -.01, length.out = 20), seq(.01, 20, length.out = 80)), ncol = 1)
HESSclayton<-matrix(NA,nrow=length(K),ncol=1)
for(i in 1:length(K)){
  type<-claytonCopula(K[i],dim=2)
  Usim<-rCopula(reps,type)
  Usim<-1-Usim
  scoreSIM<-LLgrad_c(K[i],Usim)#�൱����delta
  HESSclayton[i]<-mean(scoreSIM^2,na.rm = T) #�൱����deltaƽ��
}
hessINFO13<-cbind(K,HESSclayton)
# ʱ�亯��
rotclayton_GAS_CL<-function(theta,data,hessINFO){
  RBAR<-0.9999
  t<-nrow(data)
  w<-theta[1]
  a<-theta[2]
  b<-theta[3]
  h<-0.0001 #����÷����ò���
  # ����rhotʱ������
  ft<-rep(0,t)
  rho<-rep(0,t)
  rho[1]<-rho0
  ft[1]<-log(1 + rho[1])
  for(i in 2:t){
    It<-ml_spline(hessINFO[,1],hessINFO[,2],rho[i-1])
    DELTAt<-(-claytonCL(rho[i-1]+h,t(data[i-1,]))--claytonCL(rho[i-1],t(data[i-1,])))/h
    drhodf<-exp(ft[i-1])
    Itildat<-It/(drhodf^2)
    DELTAtildat<-DELTAt/drhodf
    ft[i]<-w+b*ft[i-1]+a*DELTAtildat/sqrt(Itildat)
    rho[i]<- -1 + exp(ft[i])
  }
  return(rho)
}
# ������Ȼ����
rotclayton_GAS_maxlik<-function(theta,data,hessINFO){
  data<-1-data
  rho<-clayton_GAS_CL(theta,data,hessINFO)
  CL<-claytonCL(rho,data)
}