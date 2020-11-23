
library(MASS)

library("parallel")
library("doParallel")
args = commandArgs(trailingOnly=TRUE)
i=as.numeric(args[1])
registerDoParallel(cores=16)

v1wo_asy<-as.matrix(read.csv("cov2s.csv",na.strings=".",header=FALSE)[,-1])
v2wo_asy<-as.matrix(read.csv("cov1s.csv",na.strings=".",header=FALSE)[,-1])
par1<-read.csv("res2s.csv",header=FALSE)
par2<-read.csv("res1s.csv",header=FALSE)

qq<-c(0.0109517,0.21082028,0.42437949,0.71459862,1.18552185,2.054813275,2.89398752,3.92892345,5.41835505,7.08575183,10.32745692)
nr<-2968

tlist=c(0.25*i)
for (t in tlist){

####Separate model
par<-par1[,2]
med1wo<-function(qq,par,K=20,M=50,t,nr){
    bh=c(par[1:10],0)
    alpha1=par[11]
    eta=par[12]
    beta=par[13:16]
    Sigma_a=matrix(data=c(par[17],0,0,par[18]),nrow=2,ncol=2)
    vare=par[19]
    NDE<-NIE<-rep(NA,nr)
    for (xi in 1:nr){
        lam<-stepfun(qq[2:(length(qq)-1)],as.numeric(bh[-length(bh)]))
        ###generate potential mediator
        t_vec=(1:K)*(t/K)
        xt=t_vec%o%c(1)
        mu1<-cbind(1,1,xt,xt)%*%beta
        mu0<-cbind(1,0,xt,0)%*%beta
        Sigma<-diag(rep(vare,K))
        e_vec0<-mvrnorm(M,mu0,Sigma)
        e_vec1<-mvrnorm(M,mu1,Sigma)
        
        samp_a=mvrnorm(M,c(0,0),Sigma_a)
        a<-samp_a%*%t(cbind(1,xt))
        m_vec1=e_vec1+a
        m_vec0=e_vec0+a
        ccum0<-ccum1<-rep(0,M)
        for (i in 1:M){
            ee0<-stepfun(t_vec,c(exp(eta*m_vec0)[i,],0))
            ee1<-stepfun(t_vec,c(exp(eta*m_vec1)[i,],0))
            myf0<-function(u){lam(u)*ee0(u)}
            myf1<-function(u){lam(u)*ee1(u)}
            ccum0[i]<-integrate(myf0,0,t,subdivisions=2000)$value
            ccum1[i]<-integrate(myf1,0,t,subdivisions=2000)$value
        }
        surv_a=samp_a%*%c(0,0)
        S00<-exp(-diag(ccum0)%*%exp(c(alpha1*0)+surv_a))
        S01<-exp(-diag(ccum1)%*%exp(c(alpha1*0)+surv_a))
        S10<-exp(-diag(ccum0)%*%exp(c(alpha1*1)+surv_a))
        S11<-exp(-diag(ccum1)%*%exp(c(alpha1*1)+surv_a))
        NDE[xi]<-mean(S10-S00,na.rm=TRUE)
        NIE[xi]<-mean(S11-S10,na.rm=TRUE)
    }
    mNDE<-mean(NDE,na.rm=TRUE)
    mNIE<-mean(NIE,na.rm=TRUE)
    return(c(mNDE,mNIE))
}
est1wo=med1wo(qq,par,K=20,M=50,t,nr)
est1wo_boot=foreach (i=1:100,.combine=rbind)%dopar%{
	out=rep(NA, 2)
    try({set.seed(1111+i)
    bootpar=mvrnorm(1,par,v1wo_asy)
        out=med1wo(qq,bootpar,K=20,M=50,t,nr)})
	return(out)
}
write.csv(est1wo_boot,paste("est1wo_s_boot_t=",as.character(t),".csv",sep=""))
vmat1wo=var(est1wo_boot,na.rm=TRUE)
res=rbind(c(est1wo[1],sqrt(vmat1wo[1,1]),est1wo[2],sqrt(vmat1wo[2,2])),
c(est1wo[1],sqrt(vmat1wo[1,1]),est1wo[2],sqrt(vmat1wo[2,2])))
write.csv(res,paste("param_s_boot_t=",as.character(t),"_res.csv",sep=""))
}
