log using "/Users/frankcobelens1/Documents/CPCD/BMGF CTVD/Non-inferiority trials/triallogprime", replace

*RR 1.5, incid BCG arm 0.05, N=5000 per arm cf 90% power for d=1.25, specificity varies, fraction 1

*first run the simulation for RR=1 to check the power (should be around 85%)
program drop trial
program define trial, eclass
        drop _all
        set obs 10000
        gen newvac=0
        replace newvac=1 if _n>5000
        gen cases=5000*0.05 
        gen rr=1
        replace rr=1.0 if newvac==1
        gen incid=0
        gen sens=1
        gen spec=1
        gen fraction=1	
        bysort newvac: replace incid=1 if _n<=rr*sens*cases*(1+rnormal())
        bysort newvac incid: replace incid=1 if incid==0 & _n<=(1-spec)*fraction*_N
        poisson incid newvac
end
simulate _b _se, reps(1000) saving("/Users/frankcobelens1/Documents/CPCD/BMGF CTVD/Non-inferiority trials/trialspec100frac100", replace): trial
sum
gen RR=exp(incid_b_newvac)
gen higher=exp(incid_b_newvac+1.96*incid_se_newvac)
sum RR, d
sum higher, d
gen higher125=higher
recode higher125 0/1.24999=0 1.25/1000=1
tab higher125

*now run the simulations for RR=1.5 by declining specificity
program drop trial
program define trial, eclass
        drop _all
        set obs 10000
        gen newvac=0
        replace newvac=1 if _n>5000
        gen cases=5000*0.05 
        gen rr=1
        replace rr=1.5 if newvac==1
        gen incid=0
        gen sens=1
        gen spec=1
        gen fraction=1	
        bysort newvac: replace incid=1 if _n<=rr*sens*cases*(1+rnormal())
        bysort newvac incid: replace incid=1 if incid==0 & _n<=(1-spec)*fraction*_N
        poisson incid newvac
end
simulate _b _se, reps(1000) saving("/Users/frankcobelens1/Documents/CPCD/BMGF CTVD/Non-inferiority trials/trialspec100frac100", replace): trial
sum
gen RR=exp(incid_b_newvac)
gen higher=exp(incid_b_newvac+1.96*incid_se_newvac)
sum RR, d
sum higher, d
gen higher125=higher
recode higher125 0/1.24999=0 1.25/1000=1
tab higher125

program drop trial
program define trial, eclass
        drop _all
        set obs 10000
        gen newvac=0
        replace newvac=1 if _n>5000
        gen cases=5000*0.05 
        gen rr=1
        replace rr=1.5 if newvac==1
        gen incid=0
        gen sens=1
        gen spec=0.98
        gen fraction=1	
        bysort newvac: replace incid=1 if _n<=rr*sens*cases*(1+rnormal())
        bysort newvac incid: replace incid=1 if incid==0 & _n<=(1-spec)*fraction*_N
        poisson incid newvac
end
simulate _b _se, reps(1000) saving("/Users/frankcobelens1/Documents/CPCD/BMGF CTVD/Non-inferiority trials/trialspec98frac100", replace): trial
sum
gen RR=exp(incid_b_newvac)
gen higher=exp(incid_b_newvac+1.96*incid_se_newvac)
sum RR, d
sum higher, d
gen higher125=higher
recode higher125 0/1.24999=0 1.25/1000=1
tab higher125

program drop trial
program define trial, eclass
        drop _all
        set obs 10000
        gen newvac=0
        replace newvac=1 if _n>5000
        gen cases=5000*0.05 
        gen rr=1
        replace rr=1.5 if newvac==1
        gen incid=0
        gen sens=1
        gen spec=0.96
        gen fraction=1	
        bysort newvac: replace incid=1 if _n<=rr*sens*cases*(1+rnormal())
        bysort newvac incid: replace incid=1 if incid==0 & _n<=(1-spec)*fraction*_N
        poisson incid newvac
end
simulate _b _se, reps(1000) saving("/Users/frankcobelens1/Documents/CPCD/BMGF CTVD/Non-inferiority trials/trialspec96frac100", replace): trial
sum
gen RR=exp(incid_b_newvac)
gen higher=exp(incid_b_newvac+1.96*incid_se_newvac)
sum RR, d
sum higher, d
gen higher125=higher
recode higher125 0/1.24999=0 1.25/1000=1
tab higher125

program drop trial
program define trial, eclass
        drop _all
        set obs 10000
        gen newvac=0
        replace newvac=1 if _n>5000
        gen cases=5000*0.05 
        gen rr=1
        replace rr=1.5 if newvac==1
        gen incid=0
        gen sens=1
        gen spec=0.94
        gen fraction=1	
        bysort newvac: replace incid=1 if _n<=rr*sens*cases*(1+rnormal())
        bysort newvac incid: replace incid=1 if incid==0 & _n<=(1-spec)*fraction*_N
        poisson incid newvac
end
simulate _b _se, reps(1000) saving("/Users/frankcobelens1/Documents/CPCD/BMGF CTVD/Non-inferiority trials/trialspec94frac100", replace): trial
sum
gen RR=exp(incid_b_newvac)
gen higher=exp(incid_b_newvac+1.96*incid_se_newvac)
sum RR, d
sum higher, d
gen higher125=higher
recode higher125 0/1.24999=0 1.25/1000=1
tab higher125

program drop trial
program define trial, eclass
        drop _all
        set obs 10000
        gen newvac=0
        replace newvac=1 if _n>5000
        gen cases=5000*0.05 
        gen rr=1
        replace rr=1.5 if newvac==1
        gen incid=0
        gen sens=1
        gen spec=0.92
        gen fraction=1	
        bysort newvac: replace incid=1 if _n<=rr*sens*cases*(1+rnormal())
        bysort newvac incid: replace incid=1 if incid==0 & _n<=(1-spec)*fraction*_N
        poisson incid newvac
end
simulate _b _se, reps(1000) saving("/Users/frankcobelens1/Documents/CPCD/BMGF CTVD/Non-inferiority trials/trialspec92frac100", replace): trial
sum
gen RR=exp(incid_b_newvac)
gen higher=exp(incid_b_newvac+1.96*incid_se_newvac)
sum RR, d
sum higher, d
gen higher125=higher
recode higher125 0/1.24999=0 1.25/1000=1
tab higher125

program drop trial
program define trial, eclass
        drop _all
        set obs 10000
        gen newvac=0
        replace newvac=1 if _n>5000
        gen cases=5000*0.05 
        gen rr=1
        replace rr=1.5 if newvac==1
        gen incid=0
        gen sens=1
        gen spec=0.90
        gen fraction=1	
        bysort newvac: replace incid=1 if _n<=rr*sens*cases*(1+rnormal())
        bysort newvac incid: replace incid=1 if incid==0 & _n<=(1-spec)*fraction*_N
        poisson incid newvac
end
simulate _b _se, reps(1000) saving("/Users/frankcobelens1/Documents/CPCD/BMGF CTVD/Non-inferiority trials/trialspec90frac100", replace): trial
sum
gen RR=exp(incid_b_newvac)
gen higher=exp(incid_b_newvac+1.96*incid_se_newvac)
sum RR, d
sum higher, d
gen higher125=higher
recode higher125 0/1.24999=0 1.25/1000=1
tab higher125

program drop trial
program define trial, eclass
        drop _all
        set obs 10000
        gen newvac=0
        replace newvac=1 if _n>5000
        gen cases=5000*0.05 
        gen rr=1
        replace rr=1.5 if newvac==1
        gen incid=0
        gen sens=1
        gen spec=0.88
        gen fraction=1	
        bysort newvac: replace incid=1 if _n<=rr*sens*cases*(1+rnormal())
        bysort newvac incid: replace incid=1 if incid==0 & _n<=(1-spec)*fraction*_N
        poisson incid newvac
end
simulate _b _se, reps(1000) saving("/Users/frankcobelens1/Documents/CPCD/BMGF CTVD/Non-inferiority trials/trialspec88frac100", replace): trial
sum
gen RR=exp(incid_b_newvac)
gen higher=exp(incid_b_newvac+1.96*incid_se_newvac)
sum RR, d
sum higher, d
gen higher125=higher
recode higher125 0/1.24999=0 1.25/1000=1
tab higher125

program drop trial
program define trial, eclass
        drop _all
        set obs 10000
        gen newvac=0
        replace newvac=1 if _n>5000
        gen cases=5000*0.05 
        gen rr=1
        replace rr=1.5 if newvac==1
        gen incid=0
        gen sens=1
        gen spec=0.86
        gen fraction=1	
        bysort newvac: replace incid=1 if _n<=rr*sens*cases*(1+rnormal())
        bysort newvac incid: replace incid=1 if incid==0 & _n<=(1-spec)*fraction*_N
        poisson incid newvac
end
simulate _b _se, reps(1000) saving("/Users/frankcobelens1/Documents/CPCD/BMGF CTVD/Non-inferiority trials/trialspec86frac100", replace): trial
sum
gen RR=exp(incid_b_newvac)
gen higher=exp(incid_b_newvac+1.96*incid_se_newvac)
sum RR, d
sum higher, d
gen higher125=higher
recode higher125 0/1.24999=0 1.25/1000=1
tab higher125

program drop trial
program define trial, eclass
        drop _all
        set obs 10000
        gen newvac=0
        replace newvac=1 if _n>5000
        gen cases=5000*0.05 
        gen rr=1
        replace rr=1.5 if newvac==1
        gen incid=0
        gen sens=1
        gen spec=0.84
        gen fraction=1	
        bysort newvac: replace incid=1 if _n<=rr*sens*cases*(1+rnormal())
        bysort newvac incid: replace incid=1 if incid==0 & _n<=(1-spec)*fraction*_N
        poisson incid newvac
end
simulate _b _se, reps(1000) saving("/Users/frankcobelens1/Documents/CPCD/BMGF CTVD/Non-inferiority trials/trialspec84frac100", replace): trial
sum
gen RR=exp(incid_b_newvac)
gen higher=exp(incid_b_newvac+1.96*incid_se_newvac)
sum RR, d
sum higher, d
gen higher125=higher
recode higher125 0/1.24999=0 1.25/1000=1
tab higher125

program drop trial
program define trial, eclass
        drop _all
        set obs 10000
        gen newvac=0
        replace newvac=1 if _n>5000
        gen cases=5000*0.05 
        gen rr=1
        replace rr=1.5 if newvac==1
        gen incid=0
        gen sens=1
        gen spec=0.82
        gen fraction=1	
        bysort newvac: replace incid=1 if _n<=rr*sens*cases*(1+rnormal())
        bysort newvac incid: replace incid=1 if incid==0 & _n<=(1-spec)*fraction*_N
        poisson incid newvac
end
simulate _b _se, reps(1000) saving("/Users/frankcobelens1/Documents/CPCD/BMGF CTVD/Non-inferiority trials/trialspec82frac100", replace): trial
sum
gen RR=exp(incid_b_newvac)
gen higher=exp(incid_b_newvac+1.96*incid_se_newvac)
sum RR, d
sum higher, d
gen higher125=higher
recode higher125 0/1.24999=0 1.25/1000=1
tab higher125

program drop trial
program define trial, eclass
        drop _all
        set obs 10000
        gen newvac=0
        replace newvac=1 if _n>5000
        gen cases=5000*0.05 
        gen rr=1
        replace rr=1.5 if newvac==1
        gen incid=0
        gen sens=1
        gen spec=0.80
        gen fraction=1	
        bysort newvac: replace incid=1 if _n<=rr*sens*cases*(1+rnormal())
        bysort newvac incid: replace incid=1 if incid==0 & _n<=(1-spec)*fraction*_N
        poisson incid newvac
end
simulate _b _se, reps(1000) saving("/Users/frankcobelens1/Documents/CPCD/BMGF CTVD/Non-inferiority trials/trialspec80frac100", replace): trial
sum
gen RR=exp(incid_b_newvac)
gen higher=exp(incid_b_newvac+1.96*incid_se_newvac)
sum RR, d
sum higher, d
gen higher125=higher
recode higher125 0/1.24999=0 1.25/1000=1
tab higher125

log close
