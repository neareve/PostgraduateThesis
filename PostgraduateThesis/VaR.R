# ��ʷ������VaR
risk<-function(W,data=NA,alpha=.05){
  port <- data%*%W
  rVaR <- quantile(port, alpha)
}

# ��ʽԼ������
eval_g0<-function(W,data=NA,alpha=NA){
  return(sum(W)-1)
}
# ��ֵ�ƽ��ݶ�
grad<-function(W,data=NA,alpha=.05){
  n<-length(W)
  out<-W
  for(i in 0:n){
    up<-W
    dn<-W
    up[i]<-up[i]+.0001
    dn[i]<-dn[i]-.0001
    out[i]<-(risk(up,data,alpha)-risk(dn,data,alpha))/.0002
  }
  return(out)
}
# �Ż��������ݶȺ���
toOpt <- function(W,data=NA,alpha=.05){
  list(objective=risk(W,data=data,alpha=alpha),gradient=grad(W,data=data,alpha=alpha))
}
# ��ʽԼ�����ſ˱Ⱦ�������б���Ϊ1
eqCon <- function(W,data=NA,alpha=.05){
  list(constraints=eval_g0(W,data=NA,alpha=.05),jacobian=rep(1,5))
}
# �������
opts <- list(algorithm = "NLOPT_LD_SLSQP",
             xtol_rel = 1.0e-7,
             maxeval = 10000)

