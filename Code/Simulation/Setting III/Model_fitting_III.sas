libname cheng "y:\liulei\collaboration\zheng cheng\paper2\data5";
%let datadir=y:\liulei\collaboration\zheng cheng\paper2\data5;
* Model II is correct: delta =-.5;

* Result for correlated random effects;

PROC IMPORT OUT= WORK.long 
            DATAFILE= "&datadir\1.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=20; 
RUN;
data one;
set long;
retain y0 y2 y4 y6;
by id;
aa=2;
if first.id then do;
	y0=0; y2=0; y4=0; y6=0;
end;
if month=0 then y0=y;
if month=2 then y2=y;
if month=4 then y4=y;
if month=6 then y6=y;
if last.id then output;
run;
* Calculate the duration in each quantile interval, together with the indicator of event in each interval;
data four;
set one;
array quant_d {5} (0 2 4 6 8);

array dur_d {4} dur_d1-dur_d4;

array event_d {4} event_d1-event_d4;

do i=1 to 4;
	dur_d{i}=0;
end;

do i=1 to 4;
	event_d{i}=0;
end;

	do i=2 to 5;
		if fup<=quant_d{i} then do;
			event_d{i-1}=(event=2);
			dur_d{i-1}=max(0, fup-quant_d{i-1});
			i=5;
		end;
		else dur_d{i-1}=quant_d{i}-quant_d{i-1};
	end;

run;
data long2;
set long;
aa=1;
run;
data five;
set four long2;
run;
data six;
set five;
if fup=. then fup=y;
run;
proc sort data=six;
by id aa;
run;
* Model II;
proc nlmixed data=six qpoints=5 corr cov;

parms 	h01=0.025 h02=0.05 h03=0.075 h04=0.1 
		alpha1=1 beta0=0 beta1=1 beta2=.2 gamma=1 eta=0.5 delta=-.5 vara=1 vare=1 ;
bounds h01 h02 h03 h04 vara >=0;

array dur {4} dur_d1-dur_d4;
array event_d {4} event_d1-event_d4;
array baseh {4} h01-h04;
array yall {4} y0 y2 y4 y6;


if aa=1 then do;	/* likelihood for CD4 measurements */
	mu1= beta0 +  beta1 * z2 + beta2* month +  a ;
	loglik=-.5*(y-mu1)**2/vare-.5*log(2*3.14159*vare);
end;

if aa=2 then do;	/* likelihood for survival */

	sum2=0;
	do k=1 to 4;
			/* cumulative baseline hazard for time dependent Y */
		sum2=sum2 + baseh{k} * dur{k} * exp(eta * yall{k} +  delta * yall{k} * a);
	end;

	mu3= alpha1 * z2  + gamma * a ;	/* for death event */

	loglik2=-exp(mu3) * sum2;

	if event=0 then loglik= loglik2;							/*log likelihood for censoring */
	if event=2 then do;
		base_haz_d=h01 * event_d1 + h02 * event_d2 + h03 * event_d3 + h04 * event_d4  ;
		y_current=y0 * event_d1 + y2 * event_d2 + y4 * event_d3 + y6 * event_d4  ;
		mu4= alpha1 * z2  + gamma * a + eta * y_current + delta * y_current * a;
		
		loglik=log(base_haz_d) + mu4 + loglik2;	/*log likelihood for death */
	end;
end;

model fup ~ general(loglik);
random a ~ normal(0, vara) subject=id;
ods output ParameterEstimates=est3_all FitStatistics=fit3_all CorrMatParmEst=corr1 CovMatParmEst=cov1; 
run;


%macro jointmodel();

%do ii=2 %to 200;

title "Iteration &ii";

PROC IMPORT OUT= WORK.long 
            DATAFILE= "&datadir\&ii..csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=20; 
RUN;

data one;
set long;
retain y0 y2 y4 y6;
by id;
aa=2;
if first.id then do;
	y0=0; y2=0; y4=0; y6=0;
end;
if month=0 then y0=y;
if month=2 then y2=y;
if month=4 then y4=y;
if month=6 then y6=y;
if last.id then output;
run;
* Calculate the duration in each quantile interval, together with the indicator of event in each interval;
data four;
set one;
array quant_d {5} (0 2 4 6 8);

array dur_d {4} dur_d1-dur_d4;

array event_d {4} event_d1-event_d4;

do i=1 to 4;
	dur_d{i}=0;
end;

do i=1 to 4;
	event_d{i}=0;
end;

	do i=2 to 5;
		if fup<=quant_d{i} then do;
			event_d{i-1}=(event=2);
			dur_d{i-1}=max(0, fup-quant_d{i-1});
			i=5;
		end;
		else dur_d{i-1}=quant_d{i}-quant_d{i-1};
	end;

run;
data long2;
set long;
aa=1;
run;
data five;
set four long2;
run;
data six;
set five;
if fup=. then fup=y;
run;
proc sort data=six;
by id aa;
run;

proc nlmixed data=six qpoints=5 corr cov;

parms 	h01=0.025 h02=0.05 h03=0.075 h04=0.1 
		alpha1=1 beta0=0 beta1=1 beta2=.2 gamma=1 eta=0.5 delta=-.5 vara=1 vare=1;
bounds h01 h02 h03 h04 vara >=0;

array dur {4} dur_d1-dur_d4;
array event_d {4} event_d1-event_d4;
array baseh {4} h01-h04;
array yall {4} y0 y2 y4 y6;


if aa=1 then do;	/* likelihood for CD4 measurements */
	mu1= beta0 +  beta1 * z2 + beta2* month +  a ;
	loglik=-.5*(y-mu1)**2/vare-.5*log(2*3.14159*vare);
end;

if aa=2 then do;	/* likelihood for survival */

	sum2=0;
	do k=1 to 4;
			/* cumulative baseline hazard for time dependent Y */
		sum2=sum2 + baseh{k} * dur{k} * exp(eta * yall{k} +  delta * yall{k} * a);
	end;

	mu3= alpha1 * z2  + gamma * a ;	/* for death event */

	loglik2=-exp(mu3) * sum2;

	if event=0 then loglik= loglik2;							/*log likelihood for censoring */
	if event=2 then do;
		base_haz_d=h01 * event_d1 + h02 * event_d2 + h03 * event_d3 + h04 * event_d4  ;
		y_current=y0 * event_d1 + y2 * event_d2 + y4 * event_d3 + y6 * event_d4  ;
		mu4= alpha1 * z2  + gamma * a + eta * y_current + delta * y_current * a;
		
		loglik=log(base_haz_d) + mu4 + loglik2;	/*log likelihood for death */
	end;
end;

model fup ~ general(loglik);
random a ~ normal(0, vara) subject=id;
ods output ParameterEstimates=est3 FitStatistics=fit3 CorrMatParmEst=corr&ii CovMatParmEst=cov&ii; 
run;

data est3_all;
set est3_all est3;
run;

data fit3_all;
set fit3_all fit3;
run;
%end;

%mend;

%jointmodel();


data result;
set _null_;
run;
%macro organize_result(varname);
data two;
set est3_all;
if parameter="&varname";
&varname=estimate;
se_&varname=standarderror;
keep &varname se_&varname;
run;

data result;
merge result two;
run;
%mend;
%organize_result(alpha1);

%organize_result(beta0);
%organize_result(beta1);
%organize_result(beta2);
%organize_result(gamma);
%organize_result(eta);
%organize_result(delta);

%organize_result(vara);
%organize_result(vare);

%organize_result(h01);
%organize_result(h02);
%organize_result(h03);
%organize_result(h04);

data _null_;
set result;
file "y:\liulei\collaboration\zheng cheng\paper2\data5\AGQresult3.txt";
put alpha1 beta0 beta1 beta2 gamma eta delta vara  vare h01 h02 h03 h04 
se_alpha1 se_beta0 se_beta1 se_beta2  se_gamma se_eta se_delta se_vara  se_vare se_h01 se_h02 se_h03 se_h04  ; 
run;

%macro outcov();

%do ii=1 %to 200;
	data _null_;
	set cov&ii;
	file "y:\liulei\collaboration\zheng cheng\paper2\data5\cov&ii..txt";
	put parameter h01 h02 h03 h04  alpha1 beta0 beta1 beta2 gamma eta delta vara  vare ; 
	run;
%end;
%mend;

%outcov();
