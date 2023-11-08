
library(MASS)

library("parallel")
library("doParallel")
args = commandArgs(trailingOnly=TRUE)
i=as.numeric(args[1])
registerDoParallel(cores=16)

v1wo_asy<-as.matrix(read.csv("pro1wo_asy.csv",na.strings=".")[,-1])
v2wo_asy<-as.matrix(read.csv("pro2wo_asy.csv",na.strings=".")[,-1])


qq<-c(0.0109517,0.21082028,0.42437949,0.71459862,1.18552185,2.054813275,2.89398752,3.92892345,5.41835505,7.08575183,10.32745692)
nr<-2968

tlist=c(0.25*i)
for (t in tlist){

####model 1
par<-c(2.5240,3.1699,2.8707,2.0632,1.3762,1.8061,1.7174,1.6021,2.1489,3.9043,
       -0.1098,-0.03396,76.8744,-7.3243,0.7252,0.3836,337.23,12.8548,301.08,0.000862,-0.05400)
med1wo<-function(qq,par,K=20,M=50,t,nr){
    bh=c(par[1:10],0)
    alpha1=par[11]
    eta=par[12]
    beta=par[13:16]
    Sigma_a=matrix(data=c(par[17],0,0,par[18]),nrow=2,ncol=2)
    vare=par[19]
    gamma=par[20:21]
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
        #        tdif<-t_vec-c(0,t_vec[-length(t_vec)])
        ccum0<-ccum1<-rep(0,M)
        for (i in 1:M){
            ee0<-stepfun(t_vec,c(exp(eta*m_vec0)[i,],0))
            ee1<-stepfun(t_vec,c(exp(eta*m_vec1)[i,],0))
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
write.csv(est1wo_boot,paste("est1wo_boot_t=",as.character(t),".csv",sep=""))
vmat1wo=var(est1wo_boot,na.rm=TRUE)


####model 2 
par<-c(2.5636,3.1847,2.8780,2.0655,1.3689,1.8025,1.7135,1.5934,2.1306,3.8875,-0.1265,-0.04031,76.8720,-7.3263,0.7688,0.3982,337.23,12.8740,301.08,-0.00428,-0.7838,0.000159,0.01023)
med2wo<-function(qq,par,K=20,M=50,t,nr){
    bh=c(par[1:10],0)
    alpha1=par[11]
    eta=par[12]
    beta=par[13:16]
    Sigma_a=matrix(data=c(par[17],0,0,par[18]),nrow=2,ncol=2)
    vare=par[19]
    gamma=par[20:21]
    delta=par[22:23]
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
    mNDE<-mean(NDE,na.rm=TRUE)
    mNIE<-mean(NIE,na.rm=TRUE)
    return(c(mNDE,mNIE))
}
est2wo=med2wo(qq,par,K=20,M=50,t,nr)
est2wo_boot=foreach (i=1:100,.combine=rbind)%dopar%{
	out=rep(NA, 2)
    try({set.seed(1111+i)
        bootpar=mvrnorm(1,par,v2wo_asy)
        out=med2wo(qq,bootpar,K=20,M=50,t,nr)})
	return(out)
}
write.csv(est2wo_boot,paste("est2wo_boot_t=",as.character(t),".csv",sep=""))
vmat2wo=var(est2wo_boot,na.rm=TRUE)

res=rbind(c(est1wo[1],sqrt(vmat1wo[1,1]),est1wo[2],sqrt(vmat1wo[2,2])),
c(est2wo[1],sqrt(vmat2wo[1,1]),est2wo[2],sqrt(vmat2wo[2,2])))
write.csv(res,paste("param_boot_t=",as.character(t),"_res.csv",sep=""))
}
