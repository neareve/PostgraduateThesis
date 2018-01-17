library(rugarch)
source("VaR.R")

MonteCVine_VaR<-function(rho, fam, nu=NULL, start, end, alpha=.05, len=1000,
                    W=rep(.2,5)){
  Daily_VaR <- NA
  Rt <- matrix(NA,nrow = len,ncol = 5)
  for(t in start:end){
    print(paste("Counting Day",t," ",Sys.time(),sep = ""))
    w <- matrix(runif(len * 5), nrow = len, ncol = 5)
    v <- array(NA, c(len,6,6))
    x <- matrix(NA, nrow = len, ncol = 5)
    x[,1] <- v[,1,1] <- w[,1]
    for(i in 2:5){
      v[,1,i] <- w[,i]
      for(k in (i-1):1){
        v[,1,i] <- Hinv(fam[k,i-k], v[,1,i], v[,k,k], 
                        rep(rho[t,i-k,k], len), nu[k,i-k])
      }
      x[,i] <- v[,1,i]
      if(i == 5){break}
      for(j in 1:(i-1)){
        v[,j+1,i] <- Hfunc(fam[j,i-j], v[,j,i], v[,j,j],
                            rep(rho[t,i-j,j], len), nu[j,i-j])
      }
    }
    # 返回对数收益率
    for(i in 1:5){
      Rt[,i] <- qdist("sstd", p = x[,i], mu = garch_m[t,i], sigma = garch_s[t,i], 
                      skew = garch_par[1,i], shape = garch_par[2,i], 
                      lambda = garch_par[3,i])
    }
    # couting VaR
    Daily_VaR[t-start+1]<-risk(W,Rt,alpha)
  }
  return(Daily_VaR)
}

#######################################test########################
return_seq <- .2*apply(data_ret,1,sum) 
CVine95_VaR <- MonteCVine_VaR(rho_end,fam_end,nu_end,start = 1,end = l,alpha = .05) 
CVine99_VaR <- MonteCVine_VaR(rho_end,fam_end,nu_end,start = 1,end = l,alpha = .01) # alpha=.01
CVine90_VaR <- MonteCVine_VaR(rho_end,fam_end,nu_end,start = 1,end = l,alpha = .1) # alpha=.1

VaRTest(.05,return_seq, CVine95_VaR) 
VaRTest(.01,return_seq, CVine99_VaR) 
VaRTest(.1,return_seq, CVine90_VaR) 

par(mar = c(2,2,2,1))
plot(day, return_seq, cex=.7,ylim=c(min(return_seq)-.02,max(return_seq)+.02),
     xlab = "", ylab = "", main = "动态C藤copula", xaxt = "n")
axis(1, c(day[1], day[500], day[1000]), c(day[1], day[500], day[1000])) 
lines(day, return_seq)
lines(day, CVine95_VaR, col = "blue")
lines(day, CVine99_VaR, col = "red")
lines(day, CVine90_VaR, col = "green")
legend("bottomright",c("组合回报率","1%VaR","5%VaR","10%VaR"),  text.width = 200,
       y.intersp = .5,
       col=c("black","red","blue","green"),lty=c(1,1,1,1),pch=c(1,NA,NA,NA),cex=0.5, inset = 0)
