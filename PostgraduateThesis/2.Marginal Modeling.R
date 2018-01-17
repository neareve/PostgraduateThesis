library(stochvol)
library(rugarch)
library(psych)
library(xts)
#############################################################
# 生成对数收益率
mdata<-read.csv("匹配指数.csv")
mdata<-mdata[,-1]

bvsp_ret<-logret(mdata$BVSP.Adjusted)
n225_ret<-logret(mdata$N225.Adjusted)
ft100_ret<-logret(mdata$FTSE.Adjusted)
sp500_ret<-logret(mdata$GSPC.Adjusted)
ssec_ret<-logret(mdata$SSEC.Adjusted)
# 可重新使用 匹配指数.csv中扣除首日日期匹配
data_ret<-data.frame(BVSP_ret=bvsp_ret,SP500_ret=sp500_ret,FT100_ret=ft100_ret,
                     N225_ret=n225_ret,SSEC_ret=ssec_ret)
rownames(data_ret) <- mdata[-1,1]
data_rets <- as.xts(data_ret)

#var(data_ret[,5])
#skew(data_ret[,5])
#kurtosi(data_ret[,5])
#ks.test(data_ret[,4],"pnorm")
#tseries::adf.test(data_ret[,5])
#查看收益率秩相关系数
# BVSP-SP500-FT100-N225-SSEC
cor_mat<-cor(data_ret,method = "kendall")
pairs.panels(data_ret,method = "kendall")
# write.csv(cor_mat,"收益率相关矩阵.csv")
# 绘制对数收益率图
par(mfrow=c(2,3))
BVSP <- data_rets[,"BVSP_ret"]
SP500 <- data_rets[,"SP500_ret"]
FT100 <- data_rets[,"FT100_ret"]
N225 <- data_rets[,"N225_ret"]
SSEC <- data_rets[,"SSEC_ret"]
plot(BVSP, yaxis.right = F)
plot(SP500, yaxis.right = F)
plot(FT100, yaxis.right = F)
plot(N225, yaxis.right = F)
plot(SSEC, yaxis.right = F)

# marginal distribution modeling
spec <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
                   mean.model = list(armaOrder = c(1,0)), 
                   distribution.model = "sstd")
garch_par <- matrix(NA, nrow = 3, ncol = 5, 
                    dimnames = list(c("skew", "shape", "lambda"),
                                    c("bvsp", "sp500", "ft100", "n225", "ssec")))
udata <- garch_m <- garch_s <- matrix(NA, nrow = nrow(data_ret), ncol = ncol(data_ret))
fit <- list()
for(i in 1:5){
  fit[[i]] <- ugarchfit(spec, data_ret[,i])
  udata[,i] <- as.numeric(pit(fit[[i]]))
  garch_m[,i] <- as.numeric(fitted(fit[[i]]))
  garch_s[,i] <- as.numeric(sigma(fit[[i]]))
  garch_par[1,i] <- fit[[i]]@model$pars["skew",1]
  garch_par[2,i] <- fit[[i]]@model$pars["shape",1]
  garch_par[3,i] <- fit[[i]]@model$pars["lambda",1]
}
ks.test(udata[,1], "punif")
ks.test(udata[,2], "punif")
ks.test(udata[,3], "punif")
ks.test(udata[,4], "punif")
ks.test(udata[,5], "punif")

#write.csv(udata,"copula建模数据.csv")