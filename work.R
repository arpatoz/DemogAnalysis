rm(list=ls(all=TRUE))
setwd("/Users/Arpat/Documents/Work/Rotifer/0908_Shane")

load("work.Rda")

# survival analysis
library(survival)

survdata<-subset(work,surv==0)      ##subset of the data containing the days where individuals died
survdata$status=1                   ##adding a status column as there's no censored data


########################################################################


with(survdata,plot(survfit(Surv(day,status)~1),ylab="Survivorship",xlab="Days"))    ##plot of survivorship for all of the data
 
##plotting each species separately
with(survdata[survdata$sp=='M',],plot(survfit(Surv(day,status)~1),ylab="Survivorship",xlab="Days",col="blue",xlim=c(0,85)))
par(new=T)
with(survdata[survdata$sp=='A',],plot(survfit(Surv(day,status)~1),ylab="Survivorship",xlab="Days",col="red",xlim=c(0,85)))

model1<-with(survdata,(survfit(Surv(day,status)~sp)))
summary(model1)
plot(model1,lty=c(1,3),ylab="Survivorship",xlab="day")              ##plots the same graph as above but without confidence intervals



#cox proportional hazards:

##macrotrachela
summary(m1<-coxph(with(subset(survdata,sp=='M'),(Surv(day,status)~occ*trt))))
summary(m2<-coxph(with(subset(survdata,sp=='M'),(Surv(day,status)~occ+trt))))  
summary(m3<-coxph(with(subset(survdata,sp=='M'),(Surv(day,status)~occ))))           
summary(m4<-coxph(with(subset(survdata,sp=='M'),(Surv(day,status)~trt))))         
summary(m5<-coxph(with(subset(survdata,sp=='M'),(Surv(day,status)~1)))) 
anova(m1,m2)
anova(m1,m3)
anova(m1,m4)

##adineta
summary(a1<-coxph(with(subset(survdata,sp=='A'),(Surv(day,status)~occ*trt))))
summary(a2<-coxph(with(subset(survdata,sp=='A'),(Surv(day,status)~occ+trt))))  
summary(a3<-coxph(with(subset(survdata,sp=='A'),(Surv(day,status)~occ))))           
summary(a4<-coxph(with(subset(survdata,sp=='A'),(Surv(day,status)~trt))))         
summary(a5<-coxph(with(subset(survdata,sp=='A'),(Surv(day,status)~1)))) 
anova(a1,a2)
anova(a1,a3)
anova(a1,a4)

coxpoly<-function(model,newdat,br=1,xlab="",ylab="",...){
	xx<-survfit(model, newdata=newdat)
	plot(xx$surv~xx$time,type='l',axes=F,col=ifelse(br==1,1,2),xlim=c(0,80),ylim=c(0,1),xlab=xlab,ylab=ylab,...)
	polygon(c(xx$time,rev(xx$time)),c(xx$lower,rev(xx$upper)),col=ifelse(br==1,'#BFBFBF99','#FF000050'),border=NA)
	}

quartz('',12,6)
layout(matrix(c(0,2,3,4,1,2,3,4,1,5,6,7,0,5,6,7), 4, 4, byrow = TRUE))

coxpoly(a3,data.frame(occ="0"),1,xlab='days',ylab='survival',main='month 0')
par(new=T)
coxpoly(m3,data.frame(occ="0"),1)
text(c(20,60),c(0.2,0.8),c('adineta','macrotrachela'))
axis(1)
axis(2)

coxpoly(m1,data.frame(occ="4",trt='D'),1,main='month 4')
par(new=T)
coxpoly(m1,data.frame(occ="4",trt='H'),2)
mtext('macrotrachela', side=2, line=3, font=2)
legend("bottomleft",c('hydrated','dehydrated'),bty='n',lty=1,col=1:2)
axis(1)
axis(2)

coxpoly(m1,data.frame(occ="8",trt='D'),1,main='month 8')
par(new=T)
coxpoly(m1,data.frame(occ="8",trt='H'),2)
axis(1)

coxpoly(m1,data.frame(occ="12",trt='D'),1,main='month 12')
par(new=T)
coxpoly(m1,data.frame(occ="12",trt='H'),2)
axis(1)


coxpoly(a1,data.frame(occ="4",trt='D'),1)
par(new=T)
coxpoly(a1,data.frame(occ="4",trt='H'),2)
mtext('adineta', side=2, line=3, font=2)
axis(1)
axis(2)

coxpoly(a1,data.frame(occ="8",trt='D'),1)
par(new=T)
coxpoly(a1,data.frame(occ="8",trt='H'),2)
axis(1)

coxpoly(a1,data.frame(occ="12",trt='D'),1)
par(new=T)
coxpoly(a1,data.frame(occ="12",trt='H'),2)
axis(1)


summary(M00C<-coxph(with(subset(survdata,sp=='M' & occ=="0"),(Surv(day,status)~1)))) 
summary(M04D<-coxph(with(subset(survdata,sp=='M' & occ=="4"&trt=="D"),(Surv(day,status)~1)))) 
summary(M04H<-coxph(with(subset(survdata,sp=='M' & occ=="4"&trt=="H"),(Surv(day,status)~1)))) 
summary(M08D<-coxph(with(subset(survdata,sp=='M' & occ=="8"&trt=="D"),(Surv(day,status)~1)))) 
summary(M08H<-coxph(with(subset(survdata,sp=='M' & occ=="8"&trt=="H"),(Surv(day,status)~1)))) 
summary(M12D<-coxph(with(subset(survdata,sp=='M' & occ=="12"&trt=="D"),(Surv(day,status)~1)))) 
summary(M12H<-coxph(with(subset(survdata,sp=='M' & occ=="12"&trt=="H"),(Surv(day,status)~1)))) 
summary(A00C<-coxph(with(subset(survdata,sp=='A' & occ=="0"),(Surv(day,status)~1)))) 
summary(A04D<-coxph(with(subset(survdata,sp=='A' & occ=="4"&trt=="D"),(Surv(day,status)~1)))) 
summary(A04H<-coxph(with(subset(survdata,sp=='A' & occ=="4"&trt=="H"),(Surv(day,status)~1)))) 
summary(A08D<-coxph(with(subset(survdata,sp=='A' & occ=="8"&trt=="D"),(Surv(day,status)~1)))) 
summary(A08H<-coxph(with(subset(survdata,sp=='A' & occ=="8"&trt=="H"),(Surv(day,status)~1)))) 
summary(A12D<-coxph(with(subset(survdata,sp=='A' & occ=="12"&trt=="D"),(Surv(day,status)~1)))) 
summary(A12H<-coxph(with(subset(survdata,sp=='A' & occ=="12"&trt=="H"),(Surv(day,status)~1)))) 


quartz('',12,6)
layout(matrix(c(0,2,3,4,1,2,3,4,1,5,6,7,0,5,6,7), 4, 4, byrow = TRUE))

coxpoly(M00C,br=1,xlab='days',ylab='survival',main='month 0')
par(new=T)
coxpoly(A00C,br=1)
text(c(20,60),c(0.2,0.8),c('adineta','macrotrachela'))
axis(1)
axis(2)

coxpoly(M04D,br=1,main='month 4')
par(new=T)
coxpoly(M04H,br=2)
mtext('macrotrachela', side=2, line=3, font=2)
legend("bottomleft",c('hydrated','dehydrated'),bty='n',lty=1,col=1:2)
axis(1)
axis(2)

coxpoly(M08D,br=1,main='month 8')
par(new=T)
coxpoly(M08H,br=2)
axis(1)

coxpoly(M12D,br=1,main='month 12')
par(new=T)
coxpoly(M12H,br=2)
axis(1)


coxpoly(A04D,br=1)
par(new=T)
coxpoly(A04H,br=2)
mtext('adineta', side=2, line=3, font=2)
axis(1)
axis(2)

coxpoly(A08D,br=1)
par(new=T)
coxpoly(A08H,br=2)
axis(1)

coxpoly(A12D,br=1)
par(new=T)
coxpoly(A12H,br=2)
axis(1)

###################################

library(MASS)
library(lme4)
library(gplots)

work$stg2 <- work$stg
levels(work$stg2)<-list(A=c("A"),S=c("S"))
workA=subset(work,sp=='A')
workM=subset(work,sp=='M')

survA0<-glm(surv~(stg+occ+trt)^2,family=binomial,data=workA)
survA1<-stepAIC(survA0)      
survM0<-glm(surv~(stg+occ+trt)^2,family=binomial,data=workM)
survM1<-stepAIC(survM0)      

predictors<-expand.grid(trt=c(NA,levels(work$trt)),occ=levels(work$occ),stg=levels(work$stg))
predictors<-predictors[!(predictors$occ==0 & (predictors$trt=="H" | is.na(predictors$trt))),] 
predsA<-predict(survA0,predictors,interval="confidence",se.fit=T)
predsM<-predict(survM0,predictors,interval="confidence",se.fit=T)

survA<-data.frame(
sp = "A",
surv = plogis(predsA$fit),
surv.upper = plogis(predsA$fit + 1.96*predsA$se.fit),
surv.lower = plogis(predsA$fit - 1.96*predsA$se.fit))
survA[survA==0]<-1

survM<-data.frame(
sp = "M",
surv = plogis(predsM$fit),
surv.upper = plogis(predsM$fit + 1.96*predsM$se.fit),
surv.lower = plogis(predsM$fit - 1.96*predsM$se.fit))
survM[survM==0]<-1

Surv<-cbind(rbind(predictors,predictors),rbind(survA,survM))
rownames(Surv)<-NULL

plotCI(Surv$surv, ui=Surv$surv.upper, li=Surv$surv.lower, axes=F, ylab="Survival", xlab="", pch=16, gap=0, ylim=c(0.7,1), col=rep(c(1,1,2,3,2,3,2,3),3))

plot.surv <- function(X,Y, ylab="", xlab="", ...){
	with(subset(Surv,sp==X & stg==Y), plotCI(surv, ui=surv.upper, li=surv.lower, axes=F, ylab=ylab, xlab=xlab, pch=16, gap=0, xlim=c(0,10), ylim=c(0.85,1), col=rep(c(1,2,3,2,3,2,3),3), ...))
	axis(1, at = c(1,3.5,6.5,9.5), labels= c('month0','month4','month8','month12'))
	axis(2)
	}


quartz('',12,6)
layout(matrix(1:6, 2, 3, byrow = TRUE))

plot.surv('A','J',ylab='daily survival')
plot.surv('A','A')
plot.surv('A','S')
plot.surv('M','J',ylab='daily survival')
plot.surv('M','A')
plot.surv('M','S')


