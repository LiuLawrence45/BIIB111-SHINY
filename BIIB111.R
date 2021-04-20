rm(list=ls())

library(VCA)
library(readxl)
library(plotrix)

path = c("/Users/lawrenceliu/Desktop/stat4ward/biib111")
setwd(path)

data = data.frame(read_excel('NSR-REP-02_Covance_ELISPOT_PROD_02JUL2020_v1.xlsx'))

# exclude cv.s.t<80%
data.0=data[data$Analyte.ID=='% Cell Viability',]
data.0=data.0[data.0$Lab.Result.Value>80,]
data.0=data.0[,c(2,13)]
colnames(data.0)=c('subject','time.point')

data$assay=substr(data$Analyte.ID,1,nchar(data$Analyte.ID)-6)
data=data[,c(2,13,17,19,23)]
colnames(data)=c('subject','time.point','analyte','value','assay')
data=data[substr(data$assay,1,3)!='CD3',]
data=data[(grep('Rep', data$analyte, fixed = T)),]
data=data[-grep('Human',data$analyte,fixed=T),]
data$value=as.numeric(data$value)
data=data[!is.na(data$value),]

data=merge(data,data.0,by=c('subject','time.point'))



############### varplot for each subject
n=length(unique(data$subject))
id=unique(data$subject)
pdf("BIIB111 Varplots for all assays.pdf", width = 12, height = 6)
for (i in 1:n){
  temp=data[data$subject==id[i],]
  temp=temp[!is.na(temp$value),]
  varPlot(value~time.point/assay,temp,Title=list(main=paste0(id[i])),
          Boxplot=list(jitter=1),VarLab = list(las=2, col='black'),
          Points=list(pch=16, cex=0.5), YLabel = list(text='Count',cex=0.8)
  )
}
dev.off()

#### data for each remained assays
IDs=c("AAV2","REP1 Pool 1","REP1 Pool 2","REP1 Pool 3","Neg Control")
dat=data.frame()
for (i in 1:length(IDs)){
  temp=data[data$assay==IDs[i],]
  dat=rbind(dat,temp)
}

data=dat
data=data[order(data$subject),]

########
data.neg.his=data[data$assay=='Neg Control',]
hist(data.neg.his$value ,xlab='Spot Counts',ylab='Frequency',main='Histgram of NC Spot Counts', breaks=20)
box()


############### plot for each remained subject
n=length(unique(data$subject))
id=unique(data$subject)

pdf("BIIB111 Varplots for all remained assays.pdf", width = 12, height = 6)
for (i in 1:n){
  temp=data[data$subject==id[i],]
  temp=temp[!is.na(temp$value),]
  varPlot(value~time.point/assay,temp,Title=list(main=paste0(id[i])),
          Boxplot=list(jitter=1),VarLab = list(las=2, col='black'),
          Points=list(pch=16, cex=0.5), YLabel = list(text='Count',cex=0.8)
  )
}
dev.off()






##########plot mean.cv
n=length(unique(data$assay))
id=unique(data$assay)
id.color=c('Green','pink','chocolate1','blue','black')
cvs=data.frame()

for (i in 1:n){
  temp=data[data$assay==id[i],]
  temp$ID=paste0(temp$subject, temp$time.point,temp$assay)
  mean.temp = tapply(temp$value, temp$ID,mean,na.rm=T)
  sd.temp= tapply(temp$value, temp$ID,sd,na.rm=T)
  criteria.temp= sd.temp^2/mean.temp
  cv.temp=sd.temp/mean.temp*100
  data.temp=data.frame(mean=mean.temp,cv=cv.temp,criteria=criteria.temp)
  data.temp$color=id.color[i]
  cvs=rbind(data.temp,cvs)
}
cvs$id=row.names(cvs)



# cvs=cvs[!is.na(cvs$cv),]

####### outliers ######
crit1=qchisq(0.99,df=2)
crit2=qchisq(0.999,df=2)
length(cvs[cvs$criteria>crit1,]$criteria)/length(cvs$criteria)
length(cvs[cvs$criteria>crit2,]$criteria)/length(cvs$criteria)

########## plot Variance/Mean vs mean #######
par(mar=c(5,5,5,5))
k=14
rgx=range(cvs$mean)
rg1=rgx[1]-(rgx[2]-rgx[1])/k
rg2=rgx[2]+(rgx[2]-rgx[1])/k
rgy=range(cvs$criteria)
rg3=rgy[1]-(rgy[2]-rgy[1])/k
rg4=rgy[2]+(rgy[2]-rgy[1])/k
plot(cvs$mean, cvs$criteria, xlim=c(0.1, rg2), ylim=c(0.01, rg4), xlab="Mean",
     ylab = 'Variance/Mean', log='xy',col=cvs$color, pch=19, cex=0.2)
ablineclip(h=crit1,col='black',lty=2)
ablineclip(h=crit2,col='red',lty=2)
title("BIIB111")
legend('bottomleft', id, col=id.color, pch=19, cex=0.5)
legend('topleft', c('α=0.01','α=0.001'),lty=c(2,2), col=c('black','red'), pch=19, cex=0.5)



###### plot mean vs cv #####
par(mar=c(5,5,5,5))

k=14
rgx=range(cvs$mean)
rg1=rgx[1]-(rgx[2]-rgx[1])/k
rg2=rgx[2]+(rgx[2]-rgx[1])/k
rgy=range(cvs$cv)
rg3=rgy[1]-(rgy[2]-rgy[1])/k
rg4=rgy[2]+(rgy[2]-rgy[1])/k
plot(cvs$mean, cvs$cv, xlim=c(rg1, rg2), ylim=c(rg3, rg4), xlab="Mean",
     ylab = 'CV', col=cvs$color, pch=19, cex=0.2)
lines(smooth.spline(cvs$mean, cvs$cv,df=4), col='blue', lty=1, lwd=2)
abline(h=50,lty=2)
title("BIIB111")
legend('topright', id, col=id.color, pch=19, cex=0.5)
length(cvs[cvs$cv<50,]$cv)/length(cvs$cv)
############# cut data
cvs.cut=cvs[cvs$mean<500,]

plot(cvs.cut$mean, cvs.cut$cv, xlab="Mean",xlim=c(-50,550),
     ylab = 'CV', col=cvs$color, pch=19, cex=0.2)
lines(smooth.spline(cvs$mean, cvs$cv,df=4), col='blue', lty=1, lwd=2)
abline(h=50,lty=2)
lines(smooth.spline(cvs$mean, cvs$cv,df=4), col='blue', lty=1, lwd=2)
title("BIIB111")
legend('topright', id, col=id.color, pch=19, cex=0.5)


########### and mean vs cv for neg control #######
data.neg=cvs[substr(cvs$id,nchar(cvs$id)-6,nchar(cvs$id))=='Control',]

par(mar=c(5,5,5,5))
k=14
rgx=range(data.neg$mean)
rg1=rgx[1]-(rgx[2]-rgx[1])/k
rg2=rgx[2]+(rgx[2]-rgx[1])/k
rgy=range(data.neg$cv)
rg3=rgy[1]-(rgy[2]-rgy[1])/k
rg4=rgy[2]+(rgy[2]-rgy[1])/k
plot(data.neg$mean, data.neg$cv, xlim=c(rg1, rg2), ylim=c(rg3, rg4), xlab="Mean",
     ylab = 'CV', pch=19, cex=0.2)
lines(smooth.spline(data.neg$mean, data.neg$cv,df=4), col='blue', lty=1, lwd=2)
abline(h=50,lty=2)
title("Mean vs CV Plot for NC")







############## varplot by sample
pdf("BIIB111 Varplots by sample.pdf", width = 12, height = 6)
id=unique(data$assay)
n=length(id)

for (i in 1:n){
  temp=data[data$assay==id[i],]
  varPlot(value~subject/time.point,temp,Title=list(main=paste(id[i])),
          Boxplot=list(jitter=1),VarLab = list(las=2, col='black'),
          Points=list(pch=16, cex=0.5), YLabel = list(text='Value',cex=0.8)
          )
}
dev.off()


######### add id column in data
data$id=paste0(data$subject, data$time.point,data$assay)
data = merge(data,cvs,by='id')

####### ####### varplot all subject
data$id.s.t=paste0(data$subject, data$time.point)
id=unique(data$id.s.t)
n=length(id)
data.n=NULL

for (i in 1:n){
  temp=data[data$id.s.t==id[i],]
  if ('Neg Control' %in% temp$assay) {
    temp$cut=temp$value/unique(temp[temp$assay=='Neg Control',]$mean)
  }
  else {temp=NULL}
  data.n=rbind(data.n,temp)
}


varPlot(cut~subject/time.point/assay,data.n,Title=list(main=paste("AAV2,","REP1 Pool 1,","REP1 Pool 2,","REP1 Pool 3")),
        Boxplot=list(jitter=1),VarLab = list(las=2, col='black'),
        Points=list(pch=16, cex=0.5), YLabel = list(text='Mean(Trt)/Mean(NC)',cex=0.8),
        ylim=c(0,1.1*max(data.n$cut))
)
abline(h=3, lty=2, col='purple')



####### ####### plot for each remained subject 
n=length(unique(data.n$subject))
id=unique(data.n$subject)

pdf("ex by sample.pdf", width = 12, height = 6)
for (i in 1:n){
  temp=data.n[data.n$subject==id[i],]
  temp=temp[!is.na(temp$value),]
  varPlot(value~time.point/assay,temp,Title=list(main=paste0(id[i])),
          Boxplot=list(jitter=1),VarLab = list(las=2, col='black'),
          Points=list(pch=16, cex=0.5), YLabel = list(text='Count',cex=0.8)
  )
  # abline(h=temp[temp$assay=="Neg Control",]$mean[1]*2.2, lty=2, col='green')
}
dev.off()

####### ####### example for plot for each remained subject
k.x=0.85
k.y=1.15

dat.exp=data.n[data.n$subject==id[1],]
n=length(unique(dat.exp$time.point))
id=unique(dat.exp$time.point)

pdf("ex seprate by sample.pdf", width = 12, height = 6)
for (i in 1:n){
  temp=dat.exp[dat.exp$time.point==id[i],]
  temp=temp[!is.na(temp$value),]
  h=temp[temp$assay=="Neg Control",]$mean*3
  varPlot(value~assay,temp,
          ylim=c(k.x*min(temp$value),k.y*h),
          Title=list(main=paste0(unique(temp$subject),' ',id[i])),
          Boxplot=list(jitter=1),VarLab = list(las=2, col='black'),
          Points=list(pch=16, cex=0.5), YLabel = list(text='Count',cex=0.8)
  )
  abline(h=h, lty=2, col='green')
}
dev.off()