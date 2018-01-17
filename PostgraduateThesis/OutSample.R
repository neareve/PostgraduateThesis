library(quantmod)
library(data.table)
library(stochvol)
source("One_Step.R")
#提取数据
ssec<-getSymbols("^SSEC",from ="2017-01-01",to = "2017-06-30",auto.assign = F)
sp500<-getSymbols("^GSPC",from ="2017-01-01",to = "2017-06-30",auto.assign = F)
ft100<-getSymbols("^FTSE",from ="2017-01-01",to = "2017-06-30",auto.assign = F)
n225<-getSymbols("^N225",from ="2017-01-01",to = "2017-06-30",auto.assign = F)
bvsp<-getSymbols("^BVSP",from ="2017-01-01",to = "2017-06-30",auto.assign = F)

ssec<-as.data.table(ssec[,6])
sp500<-as.data.table(sp500[,6])
ft100<-as.data.table(ft100[,6])
n225<-as.data.table(n225[,6])
bvsp<-as.data.table(bvsp[,6])

data_handle <- function(x,y){
  lx<-nrow(x)
  ly<-nrow(y)
  if(lx>ly){
    x_y<-x[index%in%y$index]
    temp<-y[index%in%x_y$index]
    x_y<-cbind(x_y,temp[,-1])
  }
  else{
    x_y<-y[index%in%x$index]
    temp<-x[index%in%x_y$index]
    x_y<-cbind(x_y,temp[,-1])
  }
}
a<-data_handle(sp500,ssec)
b<-data_handle(ft100,a)
c<-data_handle(n225,b)
new_data<-data_handle(bvsp,c)
new_data <- as.data.frame(na.omit(new_data))
day_out <- as.Date(new_data[-1,1])
write.csv(new_data, "OutSampleData.csv")

bvsp_ret<-logret(new_data$BVSP.Adjusted)
n225_ret<-logret(new_data$N225.Adjusted)
ft100_ret<-logret(new_data$FTSE.Adjusted)
sp500_ret<-logret(new_data$GSPC.Adjusted)
ssec_ret<-logret(new_data$SSEC.Adjusted)
new_ret<-data.frame(BVSP_ret=bvsp_ret,SP500_ret=sp500_ret,FT100_ret=ft100_ret,
                    N225_ret=n225_ret,SSEC_ret=ssec_ret)
write.csv(new_ret, "OutSample.csv")
new_ret <- read.csv("OutSample.csv")[,-1]

# marginal model
#garch_par_out <- matrix(NA, nrow = 3, ncol = 5, 
#                    dimnames = list(c("skew", "shape", "lambda"),
#                                    c("bvsp", "sp500", "ft100", "n225", "ssec")))
udata_out <- garch_m_out <- garch_s_out <- matrix(NA, nrow = nrow(new_ret), ncol = ncol(new_ret))
#for(i in 1:5){
# fit_out <- ugarchfit(spec, new_ret[,i])
#  udata_out[,i] <- as.numeric(pit(fit_out))
#  garch_m_out[,i] <- as.numeric(fitted(fit_out))
#  garch_s_out[,i] <- as.numeric(sigma(fit_out))
#  garch_par_out[1,i] <- fit_out@model$pars["skew",1]
#  garch_par_out[2,i] <- fit_out@model$pars["shape",1]
#  garch_par_out[3,i] <- fit_out@model$pars["lambda",1]
#}
ks.test(udata_out[,1], "punif")
ks.test(udata_out[,2], "punif")
ks.test(udata_out[,3], "punif")
ks.test(udata_out[,4], "punif")
ks.test(udata_out[,5], "punif")
udata_out <- as.data.frame(udata_out)

for(i in 1:5){
  spec2 <- spec
  setfixed(spec2) <- as.list(coef(fit[[i]]))
  filt <- ugarchfilter(spec2, new_ret[,i])
  udata_out[,i] <- pit(filt)
  garch_m_out[,i] <- fitted(filt)
  garch_s_out[,i] <- sigma(filt)
}


# initial step, t+1 day
t_out <- nrow(udata_out)
rho_new <- array(NA, c(t_out, 4, 4))
rho_new[1,1,1] <- onestep(theta_end[,1,1],rho_end[l,1,1],udata[l,1:2],hessM[,,1,1],nu_end[1,1])
rho_new[1,2,1] <- onestep(theta_end[,2,1],rho_end[l,2,1],udata[l,2:3],hessM[,,2,1],nu_end[1,2])

# iteration
for(j in 2:t_out){
  rho_new[j,1,1] <- onestep(theta_end[,1,1], rho_new[j-1,1,1], udata_out[j-1,1:2],
                            hessM[,,1,1], nu_end[1,1])
}
for(j in 2:t_out){
  rho_new[j,2,1] <- onestep(theta_end[,2,1], rho_new[j-1,2,1], udata_out[j-1,2:3],
                            hessM[,,2,1], nu_end[1,2])
}

rho_new[,3:4,1] <- rho_end[1:t_out,3:4,1]
rho_new[,,2:4] <- rho_end[1:t_out,,2:4]

######################################################################################  
library(copula)
source("VaR.R")
#############################################################
############################# Muti t-copula Estimation##########################
(normal_m<-fitCopula(normalCopula(dim=5),udata))
(t_m<-fitCopula(tCopula(dim=5),udata_out)) # 选择多元t分布
tcop<-tCopula(param = t_m@copula@parameters[1],dim=5,df=t_m@copula@parameters[2])
(clayton_m<-fitCopula(claytonCopula(dim=5),udata))
(gumbel_m<-fitCopula(gumbelCopula(dim=5),udata))
(frank_m<-fitCopula(frankCopula(dim=5),udata))
(joe_m<-fitCopula(joeCopula(dim=5),udata))

# Simulation
Monte_t_VaR<-function(start,end,len=1000,alpha=.05,W=rep(.2,5)){
  Daily_VaR <- NA
  Rt<-matrix(NA,nrow = len,ncol = 5)
  for(t in start:end){
    print(paste("Counting Day",t," ",Sys.time(),sep = ""))
    x<-rCopula(len,tcop)
    # trace back log return
    for(i in 1:5){
      Rt[,i] <- qdist("sstd", p = x[,i], mu = garch_m_out[t,i], sigma = garch_s_out[t,i], 
                      skew = garch_par[1,i], shape = garch_par[2,i], 
                      lambda = garch_par[3,i])
    }
    Daily_VaR[t-start+1]<-risk(W,Rt,alpha)
  }
  return(Daily_VaR)
}

try<-Monte_t_VaR(1,t_out)
try2<-Monte_t_VaR(1,t_out,alpha = .01)
try3<-Monte_t_VaR(1,t_out,alpha = .1)
plot(day_out, return_seq_out,xlab="",
     ylim=c(min(c(return_seq_out,try,try2,try3)-.02),
            max(c(return_seq_out+0.02,try,try2,try3))),
     ylab = "",main="静态多元Copula",cex=0.7, xaxt = "n")
axis(1, c(day_out[1], day_out[50], day_out[100]), c(day_out[1], day_out[50], day_out[100]))
lines(day_out, return_seq_out)
lines(day_out, try, col="blue") #5%
lines(day_out, try2, col = "red")
lines(day_out, try3, col = "green")
legend("bottomright",c("组合回报率","1%VaR","5%VaR","10%VaR"),  text.width = 40,
       y.intersp = .5, 
       col=c("black","red","blue","green"),lty=c(1,1,1,1),pch=c(1,NA,NA,NA),cex=0.5, inset = 0)
VaRTest(.05,return_seq_out,try) 
VaRTest(.01,return_seq_out,try2) 
VaRTest(.1,return_seq_out,try3) 

################################### static Vine Copula Estimation##################
VineCop<-CDVine::CDVineCopSelect(udata_out, type = 2, familyset = c(1:6, 13, 14, 16),
                                 selectioncrit = "BIC",indeptest = T)
par<-VineCop$par
par2<-VineCop$par2
family_S<-VineCop$family
l <- nrow(udata_out)
Monte_sVine_VaR<-function(start,end,len = 1000,alpha=.05,W=rep(.2,5)){
  Daily_VaR<-NA
  Rt<-matrix(NA,nrow = len,ncol = 5)
  for(t in start:end){
    print(paste("Counting Day",t," ",Sys.time(),sep = ""))
    x <- CDVine::CDVineSim(len,family = family_S,par = par,par2 = par2,type = 2)
    for(i in 1:5){
      Rt[,i] <- qdist("sstd", p = x[,i], mu = garch_m_out[t,i], sigma = garch_s_out[t,i], 
                      skew = garch_par[1,i], shape = garch_par[2,i], 
                      lambda = garch_par[3,i])
    }
    # trace back log return
    Daily_VaR[t-start+1]<-risk(W,Rt,alpha)
  }
  return(Daily_VaR)
}

return_seq_out <- .2 * apply(new_ret, 1, sum)
plot(day_out ,return_seq,xlab="",
     ylim=c(min(c(return_seq_out,try4,try5,try6)),
            max(c(return_seq+0.02,try4,try5,try6))),
     ylab="",main="静态D藤Copula",cex=0.7, xaxt = "n")
axis(1, c(day_out[1], day_out[50], day_out[100]), c(day_out[1], day_out[50], day_out[100]))
lines(day_out, return_seq_out)
lines(day,try4, col="blue") #5%
lines(day,try5, col = "red")
lines(day,try6, col = "green")
legend("bottomright",c("组合回报率","1%VaR","5%VaR","10%VaR"),  text.width = 40,
       y.intersp = .5, 
       col=c("black","red","blue","green"),lty=c(1,1,1,1),pch=c(1,NA,NA,NA),cex=0.5, inset = 0)
try4<-Monte_sVine_VaR(1,t_out)
try5<-Monte_sVine_VaR(1,t_out,alpha = .01)
try6<-Monte_sVine_VaR(1,t_out,alpha = .1)
VaRTest(.05, return_seq_out, try4)
VaRTest(.01, return_seq_out, try5)
VaRTest(.1, return_seq_out, try6)
#########################################################################################  
library(rugarch)
source("VaR.R")
source("HInv.R") 
######################################################################
# calculate daily VaR
Monte_VaR<-function(rho,fam,nu=NULL,start,end,alpha=.05,len=1000,
                    W=rep(.2,5)){
  Daily_VaR<-NA
  Rt<-matrix(NA,nrow = len,ncol = 5)
  for(t in start:end){
    print(paste("Counting Day",t," ",Sys.time(),sep = ""))
    w<-matrix(runif(len*5),nrow = len,ncol = 5)
    v<-array(NA,c(len,6,6))
    x<-matrix(NA,nrow = len,ncol = 5)
    x[,1]<-v[,1,1]<-w[,1]
    x[,2]<-v[,1,2]<-Hinv(fam[1,1],w[,2],v[,1,1],
                         rep(rho[t,1,1],len),nu[1,1]) # first could be i
    v[,2,2]<-Hfunc1(fam[1,1],v[,1,1],v[,1,2],
                    rep(rho[t,1,1],len),nu[1,1])
    for(i in 3:5){
      v[,1,i]<-w[,i]
      for(k in (i-1):2){
        v[,1,i]<-Hinv(fam[k,i-k],v[,1,i],v[,2*k-2,i-1],
                      rep(rho[t,i-k,k],len),nu[k,i-k])
      }
      v[,1,i]<-Hinv(fam[1,i-1],v[,1,i],v[,1,i-1],
                    rep(rho[t,i-1,1],len),nu[1,i-1])
      x[,i]<-v[,1,i]
      if(i==5){break}
      v[,2,i]<-Hfunc1(fam[1,i-1],v[,1,i-1],v[,1,i],
                      rep(rho[t,i-1,1],len),nu[1,i-1])
      v[,3,i]<-Hfunc1(fam[1,i-1],v[,1,i],v[,1,i-1],
                      rep(rho[t,i-1,1],len),nu[1,i-1])
      if(i>3){
        for(j in 2:(i-2)){
          v[,2*j,i]<-Hfunc1(fam[j,i-j],v[,2*j-2,i-1],v[,2*j-1,i],
                            rep(rho[t,i-j,j],len),nu[j,i-j])
          v[,2*j+1,i]<-Hfunc1(fam[j,i-j],v[,2*j-1,i],v[,2*j-2,i-1],
                              rep(rho[t,i-j,j],len),nu[j,i-j])
        }
      }
      v[,2*i-2,i]<-Hfunc1(fam[i-1,1],v[,2*i-4,i-1],v[,2*i-3,i],
                          rep(rho[t,1,i-1],len),nu[i-1,1])
    }
    # 返回对数收益率
    for(i in 1:5){
      Rt[,i] <- qdist("sstd", p = x[,i], mu = garch_m_out[t,i], sigma = garch_s_out[t,i], 
                      skew = garch_par[1,i], shape = garch_par[2,i], 
                      lambda = garch_par[3,i])
    }
    # couting VaR
    Daily_VaR[t-start+1]<-risk(W,Rt,alpha)
  }
  return(Daily_VaR)
}

# backtest
return_seq_out <- .2*apply(new_ret,1,sum) # use weights equal to 1/5 to calculate returns
GAS_VaRm <- Monte_VaR(rho_new,fam_end,nu_end,start = 1,end = t_out,alpha = .05) # alpha=.05
GAS_VaRl <- Monte_VaR(rho_new,fam_end,nu_end,start = 1,end = t_out,alpha = .01) # alpha=.01
GAS_VaRu <- Monte_VaR(rho_new,fam_end,nu_end,start = 1,end = t_out,alpha = .1) # alpha=.1
# fail test
VaRTest(.05,return_seq_out,GAS_VaRm) 
VaRTest(.01,return_seq_out,GAS_VaRl) 
VaRTest(.1,return_seq_out,GAS_VaRu)
