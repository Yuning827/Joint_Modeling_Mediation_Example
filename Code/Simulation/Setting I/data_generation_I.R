#Define baseline cumulative hazard at each time point using parameters
Lambda.0 = function(t,par){
    (t/par[1])^par[2]
}
#Define inverse function explicit form when available
inv.Lambda.0 = function (lam, par){
    lam^(1/par[2])*par[1]
}

#Data generation function
#Lambda.0 is our assumed baseline hazard, par are parameters related to the baseline hazard
data.gen=function(beta.0, alpha.0, eta.0, gamma.0, delta.0, var.a.0, var.e.0, n.id, Lambda.0, p=1/2, J, par=NULL) {
	z=ifelse(runif(n.id, 0, 1)> p, 1, 0)  # Treatment indicator
	
	x.time=2* (0:(J-1))   # Observational times for longitudinal markers: months 0, 2, 4, 6, 8, 10
	
	a=rnorm(n.id, 0, sqrt(var.a.0))   # random effect
	
	mean.1=  beta.0[1] * z + a     #m(t) at baseline point
    
    censor.time= 6 + runif(n.id, 0, 2)
    y=rep(0, n.id* J)
    z2=y
    id=y
    month=y
    fup=y
    event=y
    follow.time=rep(NA,n.id)
    death.time=rep(NA,n.id)
    status=rep(NA,n.id)
    for (i in 1:n.id)
    {
        
        mean.u=x.time *beta.0[2] + mean.1[i]    # generate mean of longitudinal data
        y.obs =  mean.u  + rnorm(J, 0, sqrt(var.e.0)) #generate longitudinal data
        
        #Compute the cumulative hazard function for each individual given the longitudinal data
        
        factors = exp(alpha.0 * z[i] + eta.0*y.obs + gamma.0 * a[i] + delta.0 * a[i] * y.obs)  #Compute the multiplicative factors
        Lambdas = rep(0,J)  #Compute cumulative hazard at each time points x.time
        for (j in 2:J){
            Lambdas[j] = Lambdas[j-1]+( Lambda.0(x.time[j],par)-Lambda.0(x.time[j-1],par) )*factors[j-1]
        }
        
        #generate survival time
        tmp = -log(runif(1))
        if (tmp >= Lambdas[J]){
            death.time[i]=inv.Lambda.0( (tmp-Lambdas[J])/factors[J]+Lambda.0(x.time[J], par) , par)
        }
        else {
            j = min(which(Lambdas>tmp))-1
            death.time[i]=inv.Lambda.0( (tmp-Lambdas[j])/factors[j]+Lambda.0(x.time[j], par) , par)
        }
        ceil.follow.time=ceiling(min(death.time[i], censor.time[i]))
        
        
        follow.time[i]=ifelse(censor.time[i]>=death.time[i], death.time[i], censor.time[i])
        status[i]=ifelse(censor.time[i]>=death.time[i], 2 , 0)
        
        y[((i-1) * J+1): (i*J) ] = ( (0:(J-1))*2 < ceil.follow.time ) * y.obs
        z2[((i-1) * J+1): (i*J) ] = z[i]
        id[((i-1) * J+1): (i*J) ] =i
        month[((i-1) * J+1): (i*J) ] =(0:(J-1))*2
        fup[((i-1) * J+1): (i*J) ]=follow.time[i]
        event[((i-1) * J+1): (i*J) ]=status[i]
    
    }

	data.longt=data.frame(y, month, z2, id, fup, event)
	data.longt=data.longt[y != 0, ]	


	list(data.longt=data.longt, status=status, follow.time=follow.time)

}

library(flexsurv)
# True values of parameters
  n.id=200
  beta.0=c(1, .2)
  alpha.0=1
  eta.0=0.5
  gamma.0=1
  delta.0=0
  par=c(10,1)
  var.a.0=1
  var.e.0=1	


  J=4
  #n.month=6

  p=1/2

  rand.num=1001:1200

  n.rep=length(rand.num)
  total.death=rep(0, n.rep)
  total.recur=rep(0, n.rep)
  total.y=rep(0, n.rep)
  mean.y=rep(0, n.rep)

# Initial values of parameters

  gamma.init=gamma.0
  alpha.init=alpha.0

time1=date()

	filedir.1="./Code/Simulation/Setting I/simulation data I/"
#Pathway to save simulated dataset


	library(MASS)

for (ii in 1:n.rep)
{
# Generate the data

 set.seed(rand.num[ii])


  data.all=data.gen(beta.0, alpha.0, eta.0, gamma.0, delta.0, var.a.0, var.e.0, n.id, Lambda.0, p=1/2, J,par=par)

  attach(data.all)
  attach(data.longt)


	filename.1=paste(filedir.1, ii, ".csv", sep="")	# Generate comma delimited file

	data.longt.2=as.matrix(data.longt)
	write.matrix(data.longt.2, filename.1, sep=",")


	total.death[ii]=sum(data.all$status==2)
	total.y[ii]=length(data.longt$id)
	mean.y[ii]=mean(data.longt$y)
}

