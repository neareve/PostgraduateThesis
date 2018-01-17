library(quantmod)
library(data.table)
##############################################################################
#提取数据
ssec<-getSymbols("^SSEC",from ="2012-01-01",to = "2016-12-31",auto.assign = F)
sp500<-getSymbols("^GSPC",from ="2012-01-01",to = "2016-12-31",auto.assign = F)
ft100<-getSymbols("^FTSE",from ="2012-01-01",to = "2016-12-31",auto.assign = F)
n225<-getSymbols("^N225",from ="2012-01-01",to = "2016-12-31",auto.assign = F)
bvsp<-getSymbols("^BVSP",from ="2012-01-01",to = "2016-12-31",auto.assign = F)

ssec<-as.data.table(ssec[,6])
sp500<-as.data.table(sp500[,6])
ft100<-as.data.table(ft100[,6])
n225<-as.data.table(n225[,6])
bvsp<-as.data.table(bvsp[,6])

# 绘图
par(mfrow=c(2,3))
plot(ssec,type="l",xlab="",ylab="",main="上证综指")
plot(sp500,type="l",xlab="",ylab="",main="标准普尔500指数")
plot(ft100,type="l",xlab="",ylab="",main="富时100指数")
plot(n225,type="l",xlab="",ylab="",main="日经225指数")
plot(bvsp,type="l",xlab="",ylab="",main="巴西圣保罗指数")

# 整理数据，取相同的交易日期
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
data<-data_handle(bvsp,c)
write.csv(data,"匹配指数.csv")
