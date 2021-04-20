rm(list=ls())

library(VCA)
library(readxl)
library(plotrix)

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
data$value=ifelse(substr(data$value,1,1)=='>',substr(data$value,2,nchar(data$value)),data$value)
data$value=as.numeric(data$value)
data$group=ifelse(substr(data$assay,1,1)=='C','Positive','Negative')
# data$value2=2*sqrt(data$value)
data=data[!is.na(data$value),]
data=data[order(data$subject),]

data=merge(data,data.0,by=c('subject','time.point'))



data$id.s.t=paste0(data$subject, data$time.point)
id=unique(data$id.s.t)
n=length(id)
data.n=NULL

for (i in 1:n){
  temp=data[data$id.s.t==id[i],]
  if ('Neg Control' %in% temp$assay) {
    temp$cut1=temp$value/mean(temp[temp$assay=='Neg Control',]$value)
    # temp$cut2=temp$value2-mean(temp[temp$assay=='Neg Control',]$value2)
  }
  else {temp=NULL}
  data.n=rbind(data.n,temp)
}

cutoff=3

varPlot(cut1~group/subject/time.point/assay,data.n,Title=list(main=paste('Negative vs Positive')),
        Boxplot=list(jitter=1),VarLab = list(las=2, col='black'),
        Points=list(pch=16, cex=0.5), YLabel = list(text='Mean(Trt)/Mean(NC)',cex=0.8),
        ylim=c(0,1.1*max(data.n$cut1))
)
abline(h=cutoff, lty=2, col='purple')

length(data.n[data.n$group=='Negative'&data.n$cut1>cutoff,]$value)/length(data.n[data.n$group=='Negative',]$value)
length(data.n[data.n$group=='Negative'&data.n$cut1<=cutoff,]$value)/length(data.n[data.n$group=='Negative',]$value)

length(data.n[data.n$group=='Positive'&data.n$cut1>cutoff,]$value)/length(data.n[data.n$group=='Positive',]$value)
length(data.n[data.n$group=='Positive'&data.n$cut1<=cutoff,]$value)/length(data.n[data.n$group=='Positive',]$value)


# varPlot(cut2~group/subject/time.point/assay,data.n,Title=list(main=paste('Negative vs Positive')),
#         Boxplot=list(jitter=1),VarLab = list(las=2, col='black'),
#         Points=list(pch=16, cex=0.5), YLabel = list(text='Mean(Drug)-Mean(NC)',cex=0.8),
#         ylim=c(0,1.1*max(data.n$cut2))
# )
# abline(h=16.8, lty=2, col='purple')



