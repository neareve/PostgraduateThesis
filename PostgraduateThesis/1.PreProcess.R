library(quantmod)
library(data.table)
##############################################################################
#��ȡ����
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

# ��ͼ
par(mfrow=c(2,3))
plot(ssec,type="l",xlab="",ylab="",main="��֤��ָ")
plot(sp500,type="l",xlab="",ylab="",main="��׼�ն�500ָ��")
plot(ft100,type="l",xlab="",ylab="",main="��ʱ100ָ��")
plot(n225,type="l",xlab="",ylab="",main="�վ�225ָ��")
plot(bvsp,type="l",xlab="",ylab="",main="����ʥ����ָ��")

# �������ݣ�ȡ��ͬ�Ľ�������
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
write.csv(data,"ƥ��ָ��.csv")