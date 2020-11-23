libname ddiddc "Y:\liulei\collaboration\Zheng Cheng\";

* Result for correlated random effects;

data ddiddc;
set ddiddc.ddiddc;
stratum=stratum-1;
gender=gender-1;
hemobl=hemobl-12;
id=seq;
trt=randgrp-1;

run;
data one;
set ddiddc;
array cd4_all {11} cd4bl cd402 cd404 cd406 cd408 cd410 cd412 cd414 cd416 cd418 cd420;
aa=1;
stoptime=t2death/30;
status=death;
do i=1 to 11;
	cd4=cd4_all{i};
	month=2*(i-1);
	output;
end;
run;
data one2;
set one;
if cd4 ne .;
run;

data two;
set one2;
by id;
if first.id;
aa=2;
bb=1;
run;

* Get all the death event times;
data three;
set two;
if status=1;
run;

* Calculate the duration in each quantile interval, together with the indicator of event in each interval;
data five;
set two;
array quant {11} (0 2 4 6 8 10 12 14 16 18 20);
array dur {10} dur1-dur10;
array median {10} median1-median10;

array event {10} event1-event10;

do i=1 to 10;
	dur{i}=0;
	event{i}=0;
	median{i}=0;
end;

do i=2 to 11;
	if stoptime<=quant{i} then do;
		dur{i-1}=stoptime-quant{i-1};
		event{i-1}=status;
		median{i-1}=quant{i-1}+dur{i-1}/2; /* Get the median of each interval */
		lastmedian=median{i-1};
		i=11;
	end;
	else do;
		dur{i-1}=quant{i}-quant{i-1};
		median{i-1}=quant{i-1}+dur{i-1}/2; /* Get the median of each interval */
		lastmedian=median{i-1};
	end;
end;
run;
proc print data=five;
var stoptime status;
where stoptime>20;
run;
data six;
set one2 five;
log_cd4=log(cd4+1);
run;

data seven2;
set six;
if month> stoptime;
run;
proc univariate data=six;
var log_cd4;
run;

data seven;
set six;
if log_cd4=. then log_cd4=stoptime;
year=month/12;
run;
proc sort data=seven;
by id aa;
run;

data eight;
set seven;
array cdfour {10} cdfour1-cdfour10;
cdfour1=log(cd4bl+1);
cdfour2=log(cd402+1);
cdfour3=log(cd404+1);
cdfour4=log(cd406+1);
cdfour5=log(cd408+1);
cdfour6=log(cd410+1);
cdfour7=log(cd412+1);
cdfour8=log(cd414+1);
cdfour9=log(cd416+1);
cdfour10=log(cd418+1);
* Last value carried forward imputation;
* The baseline CD4 is always measured;
do i=2 to 10;
	if cdfour{i}=. then cdfour{i}=cdfour{i-1};
	last_cd4=cdfour{i};
end;
run;
title "Model I Results";

proc nlmixed data=eight qpoints=5 cov;

parms 	h1=0.02 h2=0.02 h3=0.02 h4=0.02 h5=0.02 h6=0.02 h7=0.02 h8=0.02 h9=0.02 h10=0.02
		alpha1=-.28 alpha2=.3 alpha3=1.5 alpha4=-.16 alpha5=-.5
		beta0=3.6 beta1=-.1 beta2=-.8 beta3=0 beta4=-1.3 beta5=.12 beta6=.23 beta7=0
		gamma1=-.1 gamma2=0 eta=0
		vara=1.25 varb=.6 covab=0 vare=.3 ;
bounds h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 vara varb  >=0;

array dur {10} dur1-dur10;
array baseh {10} h1-h10;
array cdfour {10} cdfour1-cdfour10;


if aa=1 then do;	/* likelihood for CD4 measurements */
	mu1= beta0 +  beta1 * trt + beta2 * year + beta3* gender + beta4*prevoi + beta5*stratum + beta6 * hemobl + beta7 * trt * year + a + b * year;
	loglik=-.5*(log_cd4-mu1)**2/vare-.5*log(vare);
end;

if aa=2 then do;
	
		sum2=0;
		do k=1 to 10;
			/* cumulative baseline hazard for time dependent CD4 measure */
			sum2=sum2 + baseh{k} * dur{k} * exp(eta * cdfour{k});
		end;

		mu3= alpha1 * trt + alpha2 * gender + alpha3* prevoi + alpha4* stratum + alpha5 * hemobl +gamma1 * a + gamma2 * b;

		loglik2=-exp(mu3) * sum2;

		if status=0 then loglik=loglik2;		/* for censoring event */
		if status=1 then do; 				/* for death event */
			base_haz_d=h1 * event1 + h2 * event2 + h3 * event3 + h4 * event4 + h5 * event5 + 
				h6 * event6 + h7 * event7 + h8* event8 +h9 * event9 + h10 * event10;
			cd4_current=cdfour1 * event1 + cdfour2 * event2 + cdfour3 * event3 + cdfour4 * event4 + cdfour5 * event5 + 
				cdfour6 * event6 + cdfour7 * event7 + cdfour8* event8 +cdfour9 * event9 + cdfour10 * event10;
			mu4= alpha1 * trt + alpha2 * gender + alpha3* prevoi + alpha4* stratum + alpha5 * hemobl + eta * cd4_current + gamma1 * a + gamma2 * b;	
			loglik=loglik2 + mu4 + log(base_haz_d);
		end;
end;

model id ~ general(loglik);
random a  b~ normal([0, 0], [vara, covab, varb]) subject=id;
ods output ParameterEstimates=est1 FitStatistics=fit1; 
estimate 'Mediation' beta1 * eta;
estimate 'eta+gamma' eta + gamma1;
predict (log_cd4-mu1)*(aa=1) out=ddiddc.cd4res3; /* Output the residual in the mixed model */
predict (aa=2) *((status=1)+loglik2) out=ddiddc.cd4mres3;	/* Output the martingale residual */

run;


* New Model IV, with interaction of  a_i and Y(t);
* Model (6) in the published paper;;
title "Model II Results";
proc nlmixed data=eight qpoints=5 cov;

parms 	h1=0.02 h2=0.02 h3=0.02 h4=0.02 h5=0.02 h6=0.02 h7=0.02 h8=0.02 h9=0.02 h10=0.02
		alpha1=-.28 alpha2=.3 alpha3=1.5 alpha4=-.16 alpha5=-.5
		beta0=3.6 beta1=-.1 beta2=-.8 beta3=0 beta4=-1.3 beta5=.12 beta6=.23 beta7=0
		gamma1=-.1 gamma2=0 eta=0 delta1=0 delta2=0
		vara=1.25 varb=.6 covab=0 vare=.3 ;
bounds h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 vara varb  >=0;

array dur {10} dur1-dur10;
array baseh {10} h1-h10;
array cdfour {10} cdfour1-cdfour10;


if aa=1 then do;	/* likelihood for CD4 measurements */
	mu1= beta0 +  beta1 * trt + beta2 * year + beta3* gender + beta4*prevoi + beta5*stratum + beta6 * hemobl + beta7 * trt * year + a + b * year ;
	loglik=-.5*(log_cd4-mu1)**2/vare-.5*log(vare);
end;

if aa=2 then do;
	
		sum2=0;
		do k=1 to 10;
			/* cumulative baseline hazard for time dependent CD4 measure */
			sum2=sum2 + baseh{k} * dur{k} * exp(eta * cdfour{k} +  delta1* cdfour{k} * a + delta2* cdfour{k} * b);
		end;

		mu3= alpha1 * trt + alpha2 * gender + alpha3* prevoi + alpha4* stratum + alpha5 * hemobl +gamma1 * a + gamma2 * b;

		loglik2=-exp(mu3) * sum2;

		if status=0 then loglik=loglik2;		/* for censoring event */
		if status=1 then do; 				/* for death event */
			base_haz_d=h1 * event1 + h2 * event2 + h3 * event3 + h4 * event4 + h5 * event5 + 
				h6 * event6 + h7 * event7 + h8* event8 +h9 * event9 + h10 * event10;
			cd4_current=cdfour1 * event1 + cdfour2 * event2 + cdfour3 * event3 + cdfour4 * event4 + cdfour5 * event5 + 
				cdfour6 * event6 + cdfour7 * event7 + cdfour8* event8 +cdfour9 * event9 + cdfour10 * event10;
			mu4= alpha1 * trt + alpha2 * gender + alpha3* prevoi + alpha4* stratum + alpha5 * hemobl 
				+ eta * cd4_current +   delta1* cd4_current *a + delta2* cd4_current *b + gamma1 * a + gamma2 * b;	
			loglik=loglik2 + mu4 + log(base_haz_d);
		end;
end;

model id ~ general(loglik);
random a  b~ normal([0, 0], [vara, covab, varb]) subject=id;
ods output ParameterEstimates=est1 FitStatistics=fit1; 
estimate 'Mediation' beta1 * eta;
estimate 'eta+gamma' eta + gamma1;
*predict (log_cd4-mu1)*(aa=1) out=ddiddc.cd4res3; /* Output the residual in the mixed model */
*predict (aa=2) *((status=1)+loglik2) out=ddiddc.cd4mres3;	/* Output the martingale residual */

run;


title "Separate Model Results";
proc nlmixed data=eight qpoints=5 cov;

parms 	h1=0.02 h2=0.02 h3=0.02 h4=0.02 h5=0.02 h6=0.02 h7=0.02 h8=0.02 h9=0.02 h10=0.02
		alpha1=-.28 alpha2=.3 alpha3=1.5 alpha4=-.16 alpha5=-.5
		beta0=3.6 beta1=-.1 beta2=-.8 beta3=0 beta4=-1.3 beta5=.12 beta6=.23 beta7=0 eta=0
		vara=1.25 varb=.6 covab=0 vare=.3 ;
bounds h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 vara varb  >=0;

array dur {10} dur1-dur10;
array baseh {10} h1-h10;
array cdfour {10} cdfour1-cdfour10;


if aa=1 then do;	/* likelihood for CD4 measurements */
	mu1= beta0 +  beta1 * trt + beta2 * year + beta3* gender + beta4*prevoi + beta5*stratum + beta6 * hemobl + beta7 * trt * year + a + b * year;
	loglik=-.5*(log_cd4-mu1)**2/vare-.5*log(vare);
end;

if aa=2 then do;
	
		sum2=0;
		do k=1 to 10;
			/* cumulative baseline hazard for time dependent CD4 measure */
			sum2=sum2 + baseh{k} * dur{k} * exp(eta * cdfour{k});
		end;

		mu3= alpha1 * trt + alpha2 * gender + alpha3* prevoi + alpha4* stratum + alpha5 * hemobl ;

		loglik2=-exp(mu3) * sum2;

		if status=0 then loglik=loglik2;		/* for censoring event */
		if status=1 then do; 				/* for death event */
			base_haz_d=h1 * event1 + h2 * event2 + h3 * event3 + h4 * event4 + h5 * event5 + 
				h6 * event6 + h7 * event7 + h8* event8 +h9 * event9 + h10 * event10;
			cd4_current=cdfour1 * event1 + cdfour2 * event2 + cdfour3 * event3 + cdfour4 * event4 + cdfour5 * event5 + 
				cdfour6 * event6 + cdfour7 * event7 + cdfour8* event8 +cdfour9 * event9 + cdfour10 * event10;
			mu4= alpha1 * trt + alpha2 * gender + alpha3* prevoi + alpha4* stratum + alpha5 * hemobl + eta * cd4_current;	
			loglik=loglik2 + mu4 + log(base_haz_d);
		end;
end;

model id ~ general(loglik);
random a  b~ normal([0, 0], [vara, covab, varb]) subject=id;
ods output ParameterEstimates=est1 FitStatistics=fit1; 
*estimate 'Mediation' beta1 * eta;
*estimate 'eta+gamma' eta + gamma1;
predict (log_cd4-mu1)*(aa=1) out=ddiddc.cd4res3; /* Output the residual in the mixed model */
predict (aa=2) *((status=1)+loglik2) out=ddiddc.cd4mres3;	/* Output the martingale residual */

run;
