
*===============================================================
 ### Analyse simulated data on vaccine study ###
===============================================================;



proc import datafile="C:\Users\eidedgri\Documents\GitHub\TB-Vacc\Output\spec100.csv"
			dbms=csv
			out=spec100;
run;

proc contents data=spec100; run;


* ### 100% specificity ###;

proc genmod data=spec100 descending;
	by rep;
	model obs_case = newvac / dist=binomial link=log lrci;
	estimate 'newvac' newvac 1;
	ods output estimates=est100;
run;

data est100_power;
	set est100;

	if meanuppercl > 1.25 then power=1;
	else power=0;

run;

proc freq data=est100_power;
	tables power;
run;

* # RD model;
proc genmod data=spec100 descending;
	by rep;
	model obs_case = newvac / dist=binomial link=identity lrci;
	estimate 'newvac' newvac 1;
	ods output estimates=est100_rd;
run;

data est100_rd_power;
	set est100_rd;

	if meanuppercl > 0.0125 then power=1;
	else power=0;

run;

proc freq data=est100_rd_power;
	tables power;
run;




* ### 95% specificity ###;

proc genmod data=spec100 descending;
	by rep;
	model obs_case95 = newvac / dist=binomial link=log lrci;
	estimate 'newvac' newvac 1;
	ods output estimates=est95;
run;

proc contents data=est95; run;

data est95_power;
	set est95;

	if meanuppercl > 1.25 then power=1;
	else power=0;

run;

proc freq data=est95_power;
	tables power;
run;



* ### RR 1.25 ###;

proc import datafile="C:\Users\eidedgri\Documents\GitHub\TB-Vacc\Output\rr125.csv"
			dbms=csv
			out=rr125;
run;

proc contents data=rr125; run;


* ### 95% specificity ###;

proc genmod data=rr125 descending;
	by rep;
	model obs_case95 = newvac / dist=binomial link=log lrci;
	estimate 'newvac' newvac 1;
	ods output estimates=est95_125;
run;

proc contents data=est95_125; run;

proc means data=rr125 n mean; class newvac; var case case_n obs_case95; run;

data est95_125_type1;
	set est95_125;

	if meanuppercl < 1.25 then type1=1;
	else type1=0;

run;

proc freq data=est95_125_type1;
	tables type1;
run;


* # RD model;
proc genmod data=rr125 descending;
	by rep;
	model obs_case95 = newvac / dist=binomial link=identity lrci;
	estimate 'newvac' newvac 1;
	ods output estimates=est95_rd_125;
run;

proc contents data=est95_rd_125; run;

data est95_rd_125_type1;
	set est95_rd_125;

	if meanuppercl < 0.0125 then type1=1;
	else type1=0;

run;

proc freq data=est95_rd_125_type1;
	tables type1;
run;
