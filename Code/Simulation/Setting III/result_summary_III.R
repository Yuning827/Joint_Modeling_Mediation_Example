library(MASS)
library(doParallel)
BB=100
KK=20
MM=50

setwd("C:/Users/cheng.zheng/Desktop/Lei")
result5<-read.csv("data5_res.csv")
names(result5)[10:13]<-c("h01","h02","h03","h04")
setwd("C:/Users/cheng.zheng/Desktop/Lei/data5")
truepar<-c(0.1,0.1,0.1,0.1,1,0,1,0.2,1,0.5,-0.5,1,1)
cr<-matrix(data=NA,nrow=200,ncol=length(truepar))

qq<-c(0,2,4,6,8)
nr<-200
tlist=c(2,4,6,8)


med5<-function(qq,par,K,M,t,nr){
    bh=c(par[1:4],0)
    alpha1=par[5]
    eta=par[10]
    beta=par[6:8]
    Sigma_a=matrix(data=c(par[12]),nrow=1,ncol=1)
    vare=par[13]
    gamma=par[9]
    delta=par[11]
    NDE<-NIE<-rep(NA,nr)
    for (xi in 1:nr){
        lam<-stepfun(qq[2:(length(qq)-1)],as.numeric(bh[-length(bh)]))
        ###generate potential mediator
        t_vec=(1:K)*(t/K)
        xt=t_vec%o%c(1)
        mu1<-cbind(1,1,xt)%*%beta
        mu0<-cbind(1,0,xt)%*%beta
        Sigma<-diag(rep(vare,K))
        e_vec0<-mvrnorm(M,mu0,Sigma)
        e_vec1<-mvrnorm(M,mu1,Sigma)
        
        #samp_a=mvrnorm(M,c(0),Sigma_a)
	samp_a=matrix(log(rgamma(M, shape=1.5, scale=1)),nrow=M,ncol=1)
        a<-samp_a
        m_vec1=e_vec1+a%*%rep(1,K)
        m_vec0=e_vec0+a%*%rep(1,K)
        #        tdif<-t_vec-c(0,t_vec[-length(t_vec)])
        ccum0<-ccum1<-rep(0,M)
        m_a=c(samp_a%*%delta)%o%rep(1,K)
        for (i in 1:M){
            ee0<-stepfun(t_vec,c(exp(eta*m_vec0+m_a*m_vec0)[i,],0))
            ee1<-stepfun(t_vec,c(exp(eta*m_vec1+m_a*m_vec1)[i,],0))
            myf0<-function(u){lam(u)*ee0(u)}
            myf1<-function(u){lam(u)*ee1(u)}
            ccum0[i]<-integrate(myf0,0,t,subdivisions=2000)$value
            ccum1[i]<-integrate(myf1,0,t,subdivisions=2000)$value
        }
        surv_a=samp_a%*%gamma
        S00<-exp(-diag(ccum0)%*%exp(c(alpha1*0)+surv_a))
        S01<-exp(-diag(ccum1)%*%exp(c(alpha1*0)+surv_a))
        S10<-exp(-diag(ccum0)%*%exp(c(alpha1*1)+surv_a))
        S11<-exp(-diag(ccum1)%*%exp(c(alpha1*1)+surv_a))
        NDE[xi]<-mean(S10-S00,na.rm=TRUE)
        NIE[xi]<-mean(S11-S10,na.rm=TRUE)
    }
    mNDE<-mean(NDE)
    mNIE<-mean(NIE)
    return(c(mNDE,mNIE))
}


true_effects=foreach (t=tlist,.combine=c)%do%{
	set.seed(1111)
	med5(qq,as.numeric(truepar),K=KK,M=MM*10,t,nr)
}

setwd("C:/Users/cheng.zheng/Desktop/Lei/data5")
estmat=matrix(data=NA,nrow=200,ncol=length(true_effects))
semat=matrix(data=NA,nrow=200,ncol=length(true_effects))

for (i in 1:200){
	tmp=read.csv(paste("boot_res5_i=",as.character(i),".csv",sep=""))
	estmat[i,]=as.numeric(tmp[1,])
	semat[i,]=apply(tmp[-1,],2,sd,na.rm=TRUE)
}
Bias=apply(estmat,2,mean,na.rm=TRUE)-true_effects
SD=apply(estmat,2,sd,na.rm=TRUE)
MeSE=apply(semat,2,median,na.rm=TRUE)
CR=apply((abs(estmat-rep(1,200)%o%true_effects)/semat<=1.96),2,mean,na.rm=TRUE)

result=rbind(true_effects,Bias,SD,MeSE,CR)
write.csv(result,"result5.csv")






