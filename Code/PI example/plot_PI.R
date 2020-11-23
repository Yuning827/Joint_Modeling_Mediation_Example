library(survival)
res1=res2=res3=sd1=sd2=sd3=matrix(data=0,nrow=9,ncol=2)
sd1t=sd2t=sd3t=rep(0,9)
for (i in 1:8){
t=i*0.25
tmpb1=read.csv(paste("est1wo_boot_t=",as.character(t),".csv",sep=""))
tmpb2=read.csv(paste("est2wo_boot_t=",as.character(t),".csv",sep=""))
tmpb3=read.csv(paste("est1wo_s_boot_t=",as.character(t),".csv",sep=""))
tmp=read.csv(paste("param_boot_t=",as.character(t),"_res.csv",sep=""))
tmp_s=read.csv(paste("param_s_boot_t=",as.character(t),"_res.csv",sep=""))
res1[i+1,]=as.numeric(tmp[1,c(2,4)])
res2[i+1,]=as.numeric(tmp[2,c(2,4)])
res3[i+1,]=as.numeric(tmp_s[1,c(2,4)])
sd1[i+1,]=as.numeric(tmp[1,c(3,5)])
sd2[i+1,]=as.numeric(tmp[2,c(3,5)])
sd3[i+1,]=(apply(tmpb3[,c(2:3)],2,quantile,probs=0.75,na.rm=TRUE)-apply(tmpb3[,c(2:3)],2,quantile,probs=0.25,na.rm=TRUE))/(qnorm(0.75)-qnorm(0.25))
sd1t[i+1]=sd(tmpb1[,2]+tmpb1[,3],na.rm=TRUE)
sd2t[i+1]=sd(tmpb2[,2]+tmpb2[,3],na.rm=TRUE)
sd3t[i+1]=quantile(tmpb3[,2]+tmpb3[,3],probs=0.75,na.rm=TRUE)-quantile(tmpb3[,2]+tmpb3[,3],probs=0.25,na.rm=TRUE)
}

LL1<-res1-1.96*sd1
LL2<-res2-1.96*sd2
LL3<-res3-1.96*sd3
UU1<-res1+1.96*sd1
UU2<-res2+1.96*sd2
UU3<-res3+1.96*sd3
####read in data to get total effect
data<-read.csv("pro.csv")
sdata<-data[which(data$Time==data$stop),]
#fit1=survfit(Surv(Time,death)~1,data=sdata[which(sdata$trt==1),])
#fit0=survfit(Surv(Time,death)~1,data=sdata[which(sdata$trt==0),])
#te=summary(fit1,time=seq(0,2,by=0.25))$surv-summary(fit0,time=seq(0,2,by=0.25))$surv
####plot
fit<-coxph(Surv(Time,death)~trt,data=data)
   bh=basehaz(fit,centered=FALSE)
    S0=exp(-bh[,1])
    S1=exp(-bh[,1])^(exp(fit$coef))
   stepS0=stepfun(bh$time, c(1, S0))
    stepS1=stepfun(bh$time, c(1, S1))

tseq<-seq(0,2,by=0.25)
te=stepS1(tseq)-stepS0(tseq)

pdf("PI.pdf")
par(mfrow=c(2,3))
plot(res1[,1]~tseq,ylim=2*c(min(cbind(LL1,LL2,LL3)),max(cbind(UU1,UU2,UU3))),type="n",xlab="Time",ylab="Effect",main="Model 1")
lines(res1[,1]~tseq,col="red",lwd=3)
lines(res1[,2]~tseq,col="blue",lwd=3)
lines(LL1[,1]~tseq,col="red",lty=2,lwd=3)
lines(UU1[,1]~tseq,col="red",lty=2,lwd=3)
lines(LL1[,2]~tseq,col="blue",lty=2,lwd=3)
lines(UU1[,2]~tseq,col="blue",lty=2,lwd=3)
lines(te~tseq,col="black",lwd=3)
lines(res1[,1]+res1[,2]~tseq,col="green",lwd=3)
lines(res1[,1]+res1[,2]+1.96*sd1t~tseq,col="green",lty=2,lwd=3)
lines(res1[,1]+res1[,2]-1.96*sd1t~tseq,col="green",lty=2,lwd=3)
legend("bottomleft",c("NDE","NIE","NIE+NDE","TE"),col=c("red","blue","green","black"),lty=1,lwd=3)
abline(h=0)

plot(res1[,1]~tseq,ylim=2*c(min(cbind(LL1,LL2,LL3)),max(cbind(UU1,UU2,UU3))),type="n",xlab="Time",ylab="Effect",main="Model 2")
lines(res2[,1]~tseq,col="red",lwd=3)
lines(res2[,2]~tseq,col="blue",lwd=3)
lines(LL2[,1]~tseq,col="red",lty=2,lwd=3)
lines(UU2[,1]~tseq,col="red",lty=2,lwd=3)
lines(LL2[,2]~tseq,col="blue",lty=2,lwd=3)
lines(UU2[,2]~tseq,col="blue",lty=2,lwd=3)
lines(te~tseq,col="black",lwd=3)
lines(res2[,1]+res2[,2]~tseq,col="green",lwd=3)
lines(res2[,1]+res2[,2]+1.96*sd2t~tseq,col="green",lty=2,lwd=3)
lines(res2[,1]+res2[,2]-1.96*sd2t~tseq,col="green",lty=2,lwd=3)
legend("bottomleft",c("NDE","NIE","NIE+NDE","TE"),col=c("red","blue","green","black"),lty=1,lwd=3)
abline(h=0)


plot(res1[,1]~tseq,ylim=2*c(min(cbind(LL1,LL2,LL3)),max(cbind(UU1,UU2,UU3))),type="n",xlab="Time",ylab="Effect",main="Separate")
lines(res3[,1]~tseq,col="red",lwd=3)
lines(res3[,2]~tseq,col="blue",lwd=3)
lines(LL3[,1]~tseq,col="red",lty=2,lwd=3)
lines(UU3[,1]~tseq,col="red",lty=2,lwd=3)
lines(UU3[,2]~tseq,col="blue",lty=2,lwd=3)
lines(LL3[,2]~tseq,col="blue",lty=2,lwd=3)
lines(te~tseq,col="black",lwd=3)
lines(res3[,1]+res3[,2]~tseq,col="green",lwd=3)
lines(res3[,1]+res3[,2]+1.96*sd3t~tseq,col="green",lty=2,lwd=3)
lines(res3[,1]+res3[,2]-1.96*sd3t~tseq,col="green",lty=2,lwd=3)
legend("bottomleft",c("NDE","NIE","NIE+NDE","TE"),col=c("red","blue","green","black"),lty=1,lwd=3)
abline(h=0)
dev.off()
