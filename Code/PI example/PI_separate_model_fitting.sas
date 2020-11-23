libname pro "C:\Users\cheng.zheng\Desktop\Lei\data2";
PROC IMPORT OUT= WORK.pro 
            DATAFILE= "C:\Users\cheng.zheng\Desktop\Lei\data2\pro.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;
proc means data=pro;
var pro;
run;
data pro;
set pro;
pro_c=pro-79; * Centered pro level;
run;
data a;
set pro;
by id;
fup=stop;
if last.id;
keep id fup;
run;
data one;
merge pro a;
by id;
aa=1;
run;
data two;
set one;
by id;
if last.id;
aa=2;
bb=1;
run;

* Get all the death event times;
data three;
set two;
if death=1;
run;
* Calculate the quantiles of the death time;
proc univariate data=three noprint;
var fup; 
output out=quant pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=q; 
run;
data quant;
set quant;
bb=1;
run;
* Merge data with the quantiles;
data four;
merge two quant;
by bb;
run;

* Calculate the duration in each quantile interval, together with the indicator of event in each interval;
data five;
set four;
array quant {11} q0 q10 q20 q30 q40 q50 q60 q70 q80 q90 q100;
array dur {10} dur1-dur10;

array event {10} event1-event10;

do i=1 to 10;
	dur{i}=0;
	event{i}=0;
end;

do i=2 to 11;
	if fup<=quant{i} then do;
		dur{i-1}=fup-quant{i-1};
		event{i-1}=death;
		i=11;
	end;
	else do;
		dur{i-1}=quant{i}-quant{i-1};
	end;
end;
run;
data five2;
set one five;
run;
proc sort data=five2;
by id aa;
run;
proc univariate data=one;
var pro;
run;
data quant2;
set quant;
aa=1;
run;
data six;
merge one quant2;
by aa;
pro_c=pro-79;
run;
* For each Pro observational time, we need to calculate the duration in each quantile interval of death time;
data six2;
set six;
array quant {11} q0 q10 q20 q30 q40 q50 q60 q70 q80 q90 q100;
array dur {10} dur1-dur10;

last_start=start;
do i=1 to 10;
	dur{i}=0;
end;

do i=2 to 11;
	if stop>quant{i} then do;
		if last_start<= quant{i} then do;
			dur{i-1}=quant{i}-last_start;
			last_start=quant{i};
		end;
	end;
	else do;
		dur{i-1}=stop - last_start;
		i=11;
	end;
end;
run;	

data six3;
set six2;
array quant {11} q0 q10 q20 q30 q40 q50 q60 q70 q80 q90 q100;
array event {10} event1-event10;
by id;
lastid=0;
if last.id then do;
	lastid=1;
	do i=1 to 10;
		event{i}=0;
	end;

	do i=2 to 11;
		if stop<=quant{i} then do;
			event{i-1}=death;
			i=11;
		end;
	end;
end;
run;
* Model I with random intercept and random slope;

* Initial value for survival model;
proc nlmixed data=six3 qpoints=5 cov corr;
parms 	h1=3 h2=3 h3=3 h4=2 h5=1.5 h6=2 h7=2 h8=2 h9=2.5 h10=2
		alpha1=-.14 eta=0;
bounds h1 h2 h3 h4 h5 h6 h7 h8 h9 h10  >=0;

array dur {10} dur1-dur10;
array baseh {10} h1-h10;


		sum2=0;
		do k=1 to 10;
			/* cumulative baseline hazard for time dependent measure */
			sum2=sum2 + baseh{k} * dur{k} * exp(eta * pro );
		end;

  	    mu3= alpha1 * trt ;
		loglik=-exp(mu3) * sum2;

		if lastid=1 then do;

			if death=1 then do; 				/* for death event */
				base_haz_d=h1 * event1 + h2 * event2 + h3 * event3 + h4 * event4 + h5 * event5 + 
					h6 * event6 + h7 * event7 + h8* event8 +h9 * event9 + h10 * event10;
				mu4= alpha1 * trt + eta * pro   ;	
				loglik=loglik + mu4 + log(base_haz_d);
			end;
		end;
	
model id ~ general(loglik);
ods output ParameterEstimates=est0a; 

run;

* Initial value for longitudinal model;

proc nlmixed data=six3 qpoints=5;
parms 	beta0=76 beta1=-6.8 beta2=1.5 beta3=0 covab=0 vara=329 varb=21;
bounds vara  varb>=0;

/* likelihood for repeated measures */
	mu1= beta0 +  beta1 * trt + beta2 * start + beta3* trt * start +  a + b * start;
	loglik=-.5*(pro-mu1)**2/vare-.5*log(vare) - .5 * log(2 * 3.1415926);

	
model id ~ general(loglik);
random a b ~ normal([0, 0], [vara, covab, varb]) subject=id;
ods output ParameterEstimates=est0b; 

run;

proc mixed data=six3 method=ml covtest ;
model pro = trt | start /s;
Random Intercept start /  subject=id type=un g;
ods output solutionf=est0b2;
ods output g=est0b3; 
run;

data est0d;
set est0a est0b;
run;
* Model I with covariance of random intercept and slope;
proc nlmixed data=six3 qpoints=5 corr cov;
	parms/data=est0d;

bounds h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 vara varb >=0;

array dur {10} dur1-dur10;
array baseh {10} h1-h10;

/* likelihood for repeated measures */
	mu1= beta0 +  beta1 * trt + beta2 * start + beta3* trt * start +  a + b * start;
	loglik1=-.5*(pro-mu1)**2/vare-.5*log(vare);

		sum2=0;
		do k=1 to 10;
			/* cumulative baseline hazard for time dependent measure */
			sum2=sum2 + baseh{k} * dur{k} * exp(eta * pro );
		end;

  	    mu3= alpha1 * trt ;
		loglik2=-exp(mu3) * sum2;

		loglik=loglik1 + loglik2;

		if lastid=1 then do;

			if death=1 then do; 				/* for death event */
				base_haz_d=h1 * event1 + h2 * event2 + h3 * event3 + h4 * event4 + h5 * event5 + 
					h6 * event6 + h7 * event7 + h8* event8 +h9 * event9 + h10 * event10;
				mu4= alpha1 * trt + eta * pro ;	
				loglik=loglik + mu4 + log(base_haz_d);
			end;
		end;
	
model id ~ general(loglik);
random a b ~ normal([0, 0], [vara, covab, varb]) subject=id;
ods output ParameterEstimates=est1 FitStatistics=fit1; 
estimate 'Mediation' beta1 * eta;
*estimate 'eta+gamma' eta + gamma1;
predict pro-a out=b;
predict a out=a;

run;

* Model II;


* Model II without covariance of random intercept and slope;
************************** Used in paper ****************************;


data est1b;
set est1;
if parameter="covab" then delete;
run;

proc nlmixed data=six3 qpoints=5 corr cov;
	parms/data=est1b;
bounds h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 vara varb >=0;

array dur {10} dur1-dur10;
array baseh {10} h1-h10;

/* likelihood for repeated measures */
	mu1= beta0 +  beta1 * trt + beta2 * start + beta3* trt * start + a + b * start;
	loglik1=-.5*(pro-mu1)**2/vare-.5*log(vare);

		sum2=0;
		do k=1 to 10;
			/* cumulative baseline hazard for time dependent measure */
			sum2=sum2 + baseh{k} * dur{k} * exp(eta * pro );
		end;

  	    mu3= alpha1 * trt ;
		loglik2=-exp(mu3) * sum2;

		loglik=loglik1 + loglik2;

		if lastid=1 then do;

			if death=1 then do; 				/* for death event */
				base_haz_d=h1 * event1 + h2 * event2 + h3 * event3 + h4 * event4 + h5 * event5 + 
					h6 * event6 + h7 * event7 + h8* event8 +h9 * event9 + h10 * event10;
				mu4= alpha1 * trt + eta * pro;	
				loglik=loglik + mu4 + log(base_haz_d);
			end;
		end;
	
model id ~ general(loglik);
random a b ~ normal([0, 0], [vara, 0, varb]) subject=id;
ods output ParameterEstimates=est4 FitStatistics=fit4; 
estimate 'Mediation' beta1 * eta;
*estimate 'eta+gamma' eta + gamma1;
run;


