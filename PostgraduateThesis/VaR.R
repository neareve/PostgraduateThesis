# 历史法计算VaR
risk<-function(W,data=NA,alpha=.05){
  port <- data%*%W
  rVaR <- quantile(port, alpha)
}

# 等式约束条件
eval_g0<-function(W,data=NA,alpha=NA){
  return(sum(W)-1)
}
# 数值逼近梯度
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
# 优化函数及梯度函数
toOpt <- function(W,data=NA,alpha=.05){
  list(objective=risk(W,data=data,alpha=alpha),gradient=grad(W,data=data,alpha=alpha))
}
# 等式约束，雅克比矩阵对所有变量为1
eqCon <- function(W,data=NA,alpha=.05){
  list(constraints=eval_g0(W,data=NA,alpha=.05),jacobian=rep(1,5))
}
# 求解设置
opts <- list(algorithm = "NLOPT_LD_SLSQP",
             xtol_rel = 1.0e-7,
             maxeval = 10000)


