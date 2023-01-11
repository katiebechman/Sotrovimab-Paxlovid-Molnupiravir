********************************************************************************
*
*	Do-file:		cox.do
*
*	Project:		sotrovimab-and-Paxlovid
*
*	Programmed by:	Bang Zheng
*
*	Data used:		output/main.dta
*
*	Output:	        logs/cox.log  output/phtest.svg  output/phtest_psw.svg
*
********************************************************************************
*
*	Purpose: This do-file implements stratified Cox regression, propensity score
*   weighted Cox, and subgroup analyses.
*  
********************************************************************************

* Open a log file
cap log close
log using ./logs/cox_update, replace t
clear

use ./output/main_update.dta

*follow-up time and events*
stset end_date ,  origin(start_date) failure(failure==1)
keep if _st==1
tab _t,m
tab _t drug,m col
by drug, sort: sum _t ,de
tab _t drug if failure==1,m col
tab _t drug if failure==1&end_date==covid_hospitalisation_outcome_da&end_date!=death_with_covid_on_the_death_ce,m col
tab _t drug if failure==1&end_date==death_with_covid_on_the_death_ce,m col
tab failure drug,m col
*check censor reasons*
tab _t drug if failure==0&_t<28&end_date==death_date,m col
tab _t drug if failure==0&_t<28&end_date==dereg_date,m col
tab _t drug if failure==0&_t<28&end_date==covid_hosp_date_day_cases_mab,m col
tab _t drug if failure==0&_t<28&end_date==min(molnupiravir_covid_therapeutics,sotrovimab_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)&drug==1,m col
tab _t drug if failure==0&_t<28&end_date==min(paxlovid_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)&drug==0,m col

*main analysis*
stcox i.drug age i.sex, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
*PH test*
estat phtest,de
estat phtest, plot(1.drug)
graph export ./output/phtest_update.svg, as(svg) replace

do "analysis/ado/psmatch2.ado"
psmatch2 drug age i.sex i.region_nhs, logit
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug age i.sex i.region_nhs ) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug

psmatch2 drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug

psmatch2 drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug

psmatch2 drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, logit
histogram _pscore, by(drug, col(1))
graph export ./output/psgraph_update.svg, as(svg) replace
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug
estat phtest,de
estat phtest, plot(1.drug)
graph export ./output/phtest_psw_update.svg, as(svg) replace


*Poisson regression*
stset end_date ,  origin(start_date) failure(failure==1)
gen fu=28
replace fu=_t-_t0 if failure==0
*tab failure fu,m
poisson failure i.drug age i.sex i.region_nhs, exposure(fu) irr
poisson failure i.drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, exposure(fu) irr
poisson failure i.drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, exposure(fu) irr
poisson failure i.drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, exposure(fu) irr

*Bayesian Poisson regression*
*set seed 2022
*bayes, saving(poisson1,replace): poisson failure i.drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, exposure(fu)
*estimates store poisson1
*bayes, saving(poisson2,replace): poisson failure  age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, exposure(fu)
*estimates store poisson2
*bayesstats ic poisson2 poisson1


*un-stratified Cox, with covariate adjustment, complete case*
stset end_date ,  origin(start_date) failure(failure==1)
stcox i.drug
stcox i.drug age i.sex
stcox i.drug age i.sex i.region_nhs
stcox i.drug age i.sex i.stp
stcox i.drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro
stcox i.drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline* 
stcox i.drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline* b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease
stcox i.drug age i.sex i.stp downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new  rare_neuro
stcox i.drug age i.sex i.stp downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new  rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline* 
stcox i.drug age i.sex i.stp downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new  rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline*  b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease
*age: 5-year band*
stcox i.drug b7.age_5y_band i.sex
stcox i.drug b7.age_5y_band i.sex i.region_nhs
stcox i.drug b7.age_5y_band i.sex i.stp
stcox i.drug b7.age_5y_band i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro
stcox i.drug b7.age_5y_band i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline* 
stcox i.drug b7.age_5y_band i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline* b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease
stcox i.drug b7.age_5y_band i.sex i.stp downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new  rare_neuro
stcox i.drug b7.age_5y_band i.sex i.stp downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new  rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline* 
stcox i.drug b7.age_5y_band i.sex i.stp downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new  rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline*  b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease
*age: Restricted cubic spline*
stcox i.drug age_spline* i.sex
stcox i.drug age_spline* i.sex i.region_nhs
stcox i.drug age_spline* i.sex i.stp
stcox i.drug age_spline* i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro
stcox i.drug age_spline* i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline* 
stcox i.drug age_spline* i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline* b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease
stcox i.drug age_spline* i.sex i.stp downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new  rare_neuro
stcox i.drug age_spline* i.sex i.stp downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new  rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline* 
stcox i.drug age_spline* i.sex i.stp downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new  rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline*  b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease
*PH test*
estat phtest,de

*un-stratified Cox, missing values as a separate category*
stcox i.drug age i.sex i.stp downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* 
stcox i.drug age i.sex i.stp downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*  b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease
stcox i.drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* 
stcox i.drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*  b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease
*age: 5-year band*
stcox i.drug b7.age_5y_band i.sex i.stp downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* 
stcox i.drug b7.age_5y_band i.sex i.stp downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*  b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease
stcox i.drug b7.age_5y_band i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* 
stcox i.drug b7.age_5y_band i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*  b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease
*age: Restricted cubic spline*
stcox i.drug age_spline* i.sex i.stp downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* 
stcox i.drug age_spline* i.sex i.stp downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*  b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease
stcox i.drug age_spline* i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* 
stcox i.drug age_spline* i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*  b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease

*stratified Cox, complete case*
*stcox i.drug age i.sex, strata(region_nhs week_after_campaign)
*too few events to allow two-level stratification*
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline* , strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline*  b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex i.region_nhs, strata(week_after_campaign)
stcox i.drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro, strata(week_after_campaign)
stcox i.drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status , strata(week_after_campaign)
stcox i.drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status  b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(week_after_campaign)
*age: 5-year band*
stcox i.drug b7.age_5y_band i.sex, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline* , strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline*  b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex i.region_nhs, strata(week_after_campaign)
stcox i.drug b7.age_5y_band i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro, strata(week_after_campaign)
stcox i.drug b7.age_5y_band i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status , strata(week_after_campaign)
stcox i.drug b7.age_5y_band i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status  b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(week_after_campaign)
*age: Restricted cubic spline*
stcox i.drug age_spline* i.sex, strata(region_nhs)
stcox i.drug age_spline* i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro, strata(region_nhs)
stcox i.drug age_spline* i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline* , strata(region_nhs)
stcox i.drug age_spline* i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline*  b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age_spline* i.sex i.region_nhs, strata(week_after_campaign)
stcox i.drug age_spline* i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro, strata(week_after_campaign)
stcox i.drug age_spline* i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status , strata(week_after_campaign)
stcox i.drug age_spline* i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status  b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(week_after_campaign)

*stratified Cox, missing values as a separate category*
stcox i.drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status, strata(week_after_campaign)
stcox i.drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(week_after_campaign)
*age: 5-year band*
stcox i.drug b7.age_5y_band i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* , strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*  b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status, strata(week_after_campaign)
stcox i.drug b7.age_5y_band i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(week_after_campaign)
*age: Restricted cubic spline*
stcox i.drug age_spline* i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* , strata(region_nhs)
stcox i.drug age_spline* i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*  b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age_spline* i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status, strata(week_after_campaign)
stcox i.drug age_spline* i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(week_after_campaign)



*propensity score weighted Cox*
*do "analysis/ado/psmatch2.ado"
*age continuous, complete case*
psmatch2 drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline*, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline*) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug

psmatch2 drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline* b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline* b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug

*age: 5-year band, complete case*
psmatch2 drug b7.age_5y_band i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline*, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug b7.age_5y_band i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline*) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug

psmatch2 drug b7.age_5y_band i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline* b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug b7.age_5y_band i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline* b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug

*age: 5-year band, missing values as a separate categorye*
psmatch2 drug b7.age_5y_band i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug b7.age_5y_band i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug

psmatch2 drug b7.age_5y_band i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug b7.age_5y_band i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug




*secondary outcomes*
*all-cause hosp/death*
*follow-up time and events*
stset end_date_allcause ,  origin(start_date) failure(failure_allcause==1)
tab _t drug,m col
by drug, sort: sum _t ,de
tab _t drug if failure_allcause==1,m col
tab _t drug if failure_allcause==1&end_date_allcause==hospitalisation_outcome_date&end_date_allcause!=death_date,m col
tab _t drug if failure_allcause==1&end_date_allcause==death_date,m col
tab failure_allcause drug if _st==1,m col
*stratified Cox, missing values as a separate category*
stcox i.drug age i.sex, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
estat phtest,de
stcox i.drug age i.sex, strata(stp)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro, strata(stp)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, strata(stp)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(stp)
stcox i.drug b7.age_5y_band i.sex, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age_spline* i.sex, strata(region_nhs)
stcox i.drug age_spline* i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro, strata(region_nhs)
stcox i.drug age_spline* i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, strata(region_nhs)
stcox i.drug age_spline* i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
*PSW*
psmatch2 drug age i.sex i.region_nhs  if _st==1, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure_allcause) (drug age i.sex i.region_nhs ) if _pscore!=.
tebalance summarize
stset end_date_allcause [pwei=psweight],  origin(start_date) failure(failure_allcause==1)
stcox i.drug

psmatch2 drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra  if _st==1, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure_allcause) (drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra) if _pscore!=.
tebalance summarize
stset end_date_allcause [pwei=psweight],  origin(start_date) failure(failure_allcause==1)
stcox i.drug

psmatch2 drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*  if _st==1, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure_allcause) (drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*) if _pscore!=.
tebalance summarize
stset end_date_allcause [pwei=psweight],  origin(start_date) failure(failure_allcause==1)
stcox i.drug

psmatch2 drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease  if _st==1, logit
histogram _pscore, by(drug, col(1))
graph export ./output/psgraph.svg, as(svg) replace
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure_allcause) (drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease) if _pscore!=.
tebalance summarize
stset end_date_allcause [pwei=psweight],  origin(start_date) failure(failure_allcause==1)
stcox i.drug
estat phtest,de

*2m covid hosp/death*
*follow-up time and events*
stset end_date_2m if start_date<=mdy(9,1,2022),  origin(start_date) failure(failure_2m==1)
tab _t,m
tab _t drug,m col
by drug, sort: sum _t ,de
tab _t drug if failure_2m==1,m col
tab _t drug if failure_2m==1&end_date_2m==covid_hospitalisation_outcome_da&end_date_2m!=death_with_covid_on_the_death_ce,m col
tab _t drug if failure_2m==1&end_date_2m==death_with_covid_on_the_death_ce,m col
tab failure_2m drug if _st==1,m col
*stratified Cox, missing values as a separate category*
stcox i.drug age i.sex, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
estat phtest,de
stcox i.drug age i.sex, strata(stp)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro, strata(stp)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, strata(stp)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(stp)
stcox i.drug b7.age_5y_band i.sex, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age_spline* i.sex, strata(region_nhs)
stcox i.drug age_spline* i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro, strata(region_nhs)
stcox i.drug age_spline* i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, strata(region_nhs)
stcox i.drug age_spline* i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
*PSW*
psmatch2 drug age i.sex i.region_nhs if _st==1, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure_2m) (drug age i.sex i.region_nhs ) if _pscore!=.
tebalance summarize
stset end_date_2m [pwei=psweight],  origin(start_date) failure(failure_2m==1)
stcox i.drug

psmatch2 drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if _st==1, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure_2m) (drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra) if _pscore!=.
tebalance summarize
stset end_date_2m [pwei=psweight],  origin(start_date) failure(failure_2m==1)
stcox i.drug

psmatch2 drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if _st==1, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure_2m) (drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*) if _pscore!=.
tebalance summarize
stset end_date_2m [pwei=psweight],  origin(start_date) failure(failure_2m==1)
stcox i.drug

psmatch2 drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if _st==1, logit
histogram _pscore, by(drug, col(1))
graph export ./output/psgraph.svg, as(svg) replace
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure_2m) (drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease) if _pscore!=.
tebalance summarize
stset end_date_2m [pwei=psweight],  origin(start_date) failure(failure_2m==1)
stcox i.drug
estat phtest,de



*subgroup analysis*
stset end_date ,  origin(start_date) failure(failure==1)
*main analysis*
stcox i.drug##i.sex age downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if sex==0, strata(region_nhs)
stcox i.drug age downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if sex==1, strata(region_nhs)

stcox i.drug##i.age_group3 i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if age_group3==0, strata(region_nhs)
stcox i.drug i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if age_group3==1, strata(region_nhs)
stcox i.drug i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if age_group3==2, strata(region_nhs)

stcox i.drug##i.age_55 i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if age_55==0, strata(region_nhs)
stcox i.drug i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if age_55==1, strata(region_nhs)

stcox i.drug##i.age_60 i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if age_60==0, strata(region_nhs)
stcox i.drug i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if age_60==1, strata(region_nhs)

stcox i.drug##i.White age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if White==1, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if White==0, strata(region_nhs)

stcox i.drug##i.downs_syndrome  solid_cancer_new age i.sex  haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease , strata(region_nhs)
stcox i.drug  solid_cancer_new age i.sex  haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if downs_syndrome==0, strata(region_nhs)
stcox i.drug##i.solid_cancer_new age i.sex downs_syndrome  haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease , strata(region_nhs)
stcox i.drug age i.sex downs_syndrome  haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if solid_cancer_new==1, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome  haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if solid_cancer_new==0, strata(region_nhs)
stcox i.drug##i.haema_disease age i.sex downs_syndrome solid_cancer_new    imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new    imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if haema_disease==1, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new    imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if haema_disease==0, strata(region_nhs)
stcox i.drug##i.imid age i.sex downs_syndrome solid_cancer_new haema_disease    immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease , strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease    immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if imid==1, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease    immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if imid==0, strata(region_nhs)
stcox i.drug##i.immunosupression_new age i.sex downs_syndrome solid_cancer_new haema_disease   imid   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease , strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if immunosupression_new==1, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if immunosupression_new==0, strata(region_nhs)
stcox i.drug##i.rare_neuro age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new  drugs_consider_risk_contra   b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease , strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   drugs_consider_risk_contra  b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if rare_neuro==0, strata(region_nhs)

stcox i.drug##i.bmi_g3 age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*  diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*  diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if bmi_g3==1, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if bmi_g3==2, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if bmi_g3==3, strata(region_nhs)

stcox i.drug##i.bmi_25 age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*  diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*  diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if bmi_25==0, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if bmi_25==1, strata(region_nhs)

stcox i.drug##i.bmi_30 age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*  diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*  diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if bmi_30==0, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if bmi_30==1, strata(region_nhs)

stcox i.drug##i.diabetes age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing  chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing chronic_cardiac_disease hypertension chronic_respiratory_disease if diabetes ==0, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing chronic_cardiac_disease hypertension chronic_respiratory_disease if diabetes ==1, strata(region_nhs)

stcox i.drug##i.chronic_cardiac_disease age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes  hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes  hypertension chronic_respiratory_disease if chronic_cardiac_disease==0, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes  hypertension chronic_respiratory_disease if chronic_cardiac_disease==1, strata(region_nhs)

stcox i.drug##i.hypertension age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease  chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease  chronic_respiratory_disease if hypertension==0, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease  chronic_respiratory_disease if hypertension==1, strata(region_nhs)

stcox i.drug##i.chronic_respiratory_disease age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension , strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension if chronic_respiratory_disease==0, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension if chronic_respiratory_disease==1, strata(region_nhs)

stcox i.drug##i.vaccination_g3 age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing  calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing  calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if vaccination_g3==1, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing  calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if vaccination_g3==2, strata(region_nhs)
*stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing  calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if vaccination_3==0, strata(region_nhs)
*stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing  calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if vaccination_status==0, strata(region_nhs)

stcox i.drug##i.d_postest_treat_g2 age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
*stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if d_postest_treat_g2==0, strata(region_nhs)
*stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if d_postest_treat_g2==1, strata(region_nhs)

gen may_31_after=(start_date>=mdy(5,31,2022))
tab drug may_31_after,row chi
stcox i.drug##i.may_31_after i.sex age downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status  b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status i.week_after_campaign b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if may_31_after==0, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status i.week_after_campaign b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if may_31_after==1, strata(region_nhs)

*use minimal-adjusted model*
stcox i.drug##i.sex age , strata(region_nhs)
stcox i.drug age if sex==0, strata(region_nhs)
stcox i.drug age if sex==1, strata(region_nhs)

stcox i.drug##i.age_group3 i.sex , strata(region_nhs)
stcox i.drug i.sex  if age_group3==0, strata(region_nhs)
stcox i.drug i.sex if age_group3==1, strata(region_nhs)
stcox i.drug i.sex if age_group3==2, strata(region_nhs)

stcox i.drug##i.age_60 i.sex , strata(region_nhs)
stcox i.drug i.sex if age_60==0, strata(region_nhs)
stcox i.drug i.sex if age_60==1, strata(region_nhs)

stcox i.drug##i.White age i.sex , strata(region_nhs)
stcox i.drug age i.sex if White==1, strata(region_nhs)
stcox i.drug age i.sex if White==0, strata(region_nhs)

stcox i.drug##i.solid_cancer_new age i.sex , strata(region_nhs)
stcox i.drug age i.sex if solid_cancer_new==1, strata(region_nhs)
stcox i.drug age i.sex if solid_cancer_new==0, strata(region_nhs)
stcox i.drug##i.haema_disease age i.sex , strata(region_nhs)
stcox i.drug age i.sex if haema_disease==1, strata(region_nhs)
stcox i.drug age i.sex if haema_disease==0, strata(region_nhs)
stcox i.drug##i.imid age i.sex , strata(region_nhs)
stcox i.drug age i.sex  if imid==1, strata(region_nhs)
stcox i.drug age i.sex  if imid==0, strata(region_nhs)
stcox i.drug##i.immunosupression_new age i.sex , strata(region_nhs)
stcox i.drug##i.rare_neuro age i.sex , strata(region_nhs)

stcox i.drug##i.bmi_g3 age i.sex , strata(region_nhs)
stcox i.drug age i.sex  if bmi_g3==1, strata(region_nhs)
stcox i.drug age i.sex if bmi_g3==2, strata(region_nhs)
stcox i.drug age i.sex if bmi_g3==3, strata(region_nhs)

stcox i.drug##i.bmi_25 age i.sex , strata(region_nhs)
stcox i.drug age i.sex if bmi_25==0, strata(region_nhs)
stcox i.drug age i.sex if bmi_25==1, strata(region_nhs)

stcox i.drug##i.bmi_30 age i.sex , strata(region_nhs)
stcox i.drug age i.sex if bmi_30==0, strata(region_nhs)
stcox i.drug age i.sex if bmi_30==1, strata(region_nhs)

stcox i.drug##i.diabetes age i.sex , strata(region_nhs)
stcox i.drug age i.sex if diabetes ==0, strata(region_nhs)
stcox i.drug age i.sex if diabetes ==1, strata(region_nhs)

stcox i.drug##i.chronic_cardiac_disease age i.sex , strata(region_nhs)
stcox i.drug age i.sex  if chronic_cardiac_disease==0, strata(region_nhs)
stcox i.drug age i.sex  if chronic_cardiac_disease==1, strata(region_nhs)

stcox i.drug##i.hypertension age i.sex , strata(region_nhs)
stcox i.drug age i.sex if hypertension==0, strata(region_nhs)
stcox i.drug age i.sex  if hypertension==1, strata(region_nhs)

stcox i.drug##i.chronic_respiratory_disease age i.sex  , strata(region_nhs)
stcox i.drug age i.sex  if chronic_respiratory_disease==0, strata(region_nhs)
stcox i.drug age i.sex  if chronic_respiratory_disease==1, strata(region_nhs)

stcox i.drug##i.vaccination_g3 age i.sex , strata(region_nhs)
stcox i.drug age i.sex  if vaccination_g3==1, strata(region_nhs)
stcox i.drug age i.sex  if vaccination_g3==2, strata(region_nhs)
*stcox i.drug age i.sex  if vaccination_3==0, strata(region_nhs)

stcox i.drug##i.d_postest_treat_g2 age i.sex , strata(region_nhs)
stcox i.drug age i.sex if d_postest_treat_g2==0, strata(region_nhs)
stcox i.drug age i.sex if d_postest_treat_g2==1, strata(region_nhs)

stcox i.drug##i.may_31_after i.sex age , strata(region_nhs)
stcox i.drug age i.sex if may_31_after==0, strata(region_nhs)
stcox i.drug age i.sex if may_31_after==1, strata(region_nhs)





*sensitivity analysis*
stset end_date ,  origin(start_date) failure(failure==1)
*additionally adjusting for days between test positive and treatment initiation, and days/months between last vaccination date and treatment initiation; *
stcox i.drug age i.sex i.d_postest_treat_missing i.month_after_vaccinate_missing, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra i.d_postest_treat_missing i.month_after_vaccinate_missing, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* i.d_postest_treat_missing i.month_after_vaccinate_missing, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease i.d_postest_treat_missing i.month_after_vaccinate_missing, strata(region_nhs)
*exclude patients with HIV/AIDS due to small number*
stcox i.drug age i.sex if hiv_aids==0, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if hiv_aids==0, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if hiv_aids==0, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if hiv_aids==0, strata(region_nhs)
*exclude 1-year drug interactions*
stcox i.drug age i.sex if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25))&(drugs_consider_risk>start_date|drugs_consider_risk<(start_date-365.25)), strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25))&(drugs_consider_risk>start_date|drugs_consider_risk<(start_date-365.25)), strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25))&(drugs_consider_risk>start_date|drugs_consider_risk<(start_date-365.25)), strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25))&(drugs_consider_risk>start_date|drugs_consider_risk<(start_date-365.25)), strata(region_nhs)
*excluding patients with treatment records of both sotrovimab and Pax, or with treatment records of any other therapies*
stcox i.drug age i.sex if (sotrovimab_covid_therapeutics==.|paxlovid_covid_therapeutics==.|sotrovimab_covid_therapeutics>start_date_29|paxlovid_covid_therapeutics>start_date_29)&molnupiravir_covid_therapeutics>start_date_29&remdesivir_covid_therapeutics>start_date_29&casirivimab_covid_therapeutics>start_date_29, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if (sotrovimab_covid_therapeutics==.|paxlovid_covid_therapeutics==.|sotrovimab_covid_therapeutics>start_date_29|paxlovid_covid_therapeutics>start_date_29)&molnupiravir_covid_therapeutics>start_date_29&remdesivir_covid_therapeutics>start_date_29&casirivimab_covid_therapeutics>start_date_29, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if (sotrovimab_covid_therapeutics==.|paxlovid_covid_therapeutics==.|sotrovimab_covid_therapeutics>start_date_29|paxlovid_covid_therapeutics>start_date_29)&molnupiravir_covid_therapeutics>start_date_29&remdesivir_covid_therapeutics>start_date_29&casirivimab_covid_therapeutics>start_date_29, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if (sotrovimab_covid_therapeutics==.|paxlovid_covid_therapeutics==.|sotrovimab_covid_therapeutics>start_date_29|paxlovid_covid_therapeutics>start_date_29)&molnupiravir_covid_therapeutics>start_date_29&remdesivir_covid_therapeutics>start_date_29&casirivimab_covid_therapeutics>start_date_29, strata(region_nhs)
*excluding patients who were identified to be pregnant at treatment initiation*
*stcox i.drug age i.sex if pregnancy!=1, strata(region_nhs)
*stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro if pregnancy!=1, strata(region_nhs)
*stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b6.ethnicity_with_missing b5.imd_with_missing i.vaccination_status i.week_after_campaign if pregnancy!=1, strata(region_nhs)
*stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b6.ethnicity_with_missing b5.imd_with_missing i.vaccination_status i.week_after_campaign b1.bmi_g4_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if pregnancy!=1, strata(region_nhs)
*additionally adjusting for rural-urban classification, other comorbidities (dementia, autism, learning disabilities, severe mental illness), and care home residency and housebound status *
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease i.rural_urban_with_missing autism_nhsd care_home_primis dementia_nhsd housebound_opensafely learning_disability_primis serious_mental_illness_nhsd, strata(region_nhs)
*excluding patients who did not have a positive SARS-CoV-2 test record before treatment or initiated treatment after 5 days since positive SARS-CoV-2 test*
stcox i.drug age i.sex if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
*create a 1-day lag in the follow-up start date *
stcox i.drug age i.sex if _t>=2, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if _t>=2, strata(region_nhs) 
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if _t>=2, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if _t>=2, strata(region_nhs)
*create a 2-day lag in the follow-up start date *
stcox i.drug age i.sex if _t>=3, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if _t>=3, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if _t>=3, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if _t>=3, strata(region_nhs)
*use older version of codelist for solid_cancer and immunosupression*
stcox i.drug age i.sex downs_syndrome solid_cancer haema_disease   imid immunosupression   rare_neuro  drugs_consider_risk_contra, strata(region_nhs) 
stcox i.drug age i.sex downs_syndrome solid_cancer haema_disease   imid immunosupression   rare_neuro  drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* , strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer haema_disease   imid immunosupression   rare_neuro  drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease , strata(region_nhs)
*using Cox models with calendar date as the underlying time scale*
gen campaign_start_date=mdy(12,16,2021)
stset end_date ,  origin(campaign_start_date) enter(start_date) failure(failure==1)
stcox i.drug age i.sex, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
*not restricting to primary diagnosis when defining COVID-19 related hospitalisation*
stset end_date_not_primary ,  origin(start_date) failure(failure_not_primary==1)
tab _t drug,m col
by drug, sort: sum _t ,de
tab _t drug if failure_not_primary==1,m col
tab _t drug if failure_not_primary==1&end_date_not_primary==death_with_covid_on_the_death_ce,m col
tab failure_not_primary drug if _st==1,m col
stcox i.drug age i.sex, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
*ATT weight*
psmatch2 drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, logit
drop psweight
gen psweight=cond( drug ==1,1,_pscore/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug
estat phtest,de
*exclude all patients in the non-overlapping parts of the PS distribution*
psmatch2 drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
sum _pscore if drug==0,de
gen _pscore_sot_min=r(min)
gen _pscore_sot_max=r(max)
sum _pscore if drug==1,de
gen _pscore_pax_min=r(min)
gen _pscore_pax_max=r(max)
stset end_date if (drug==0&_pscore>=_pscore_pax_min&_pscore<=_pscore_pax_max)|(drug==1&_pscore>=_pscore_sot_min&_pscore<=_pscore_sot_max) [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug
drop _pscore_sot_min _pscore_sot_max _pscore_pax_min _pscore_pax_max

*cause-specific analysis*
gen failure_covid=(failure_allcause==1&(covid_hospitalisation_outcome_da==end_date_allcause|death_with_covid_on_the_death_ce==end_date_allcause))
gen failure_other=(failure_allcause==1&(covid_hospitalisation_outcome_da!=end_date_allcause&death_with_covid_on_the_death_ce!=end_date_allcause))
tab failure_covid failure_other,m
stset end_date_allcause ,  origin(start_date) failure(failure_covid==1)
stcox i.drug age i.sex, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
psmatch2 drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure_covid) (drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease) if _pscore!=.
tebalance summarize
stset end_date_allcause [pwei=psweight],  origin(start_date) failure(failure_covid==1)
stcox i.drug

stset end_date_allcause ,  origin(start_date) failure(failure_other==1)
stcox i.drug age i.sex, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
psmatch2 drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure_other) (drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease) if _pscore!=.
tebalance summarize
stset end_date_allcause [pwei=psweight],  origin(start_date) failure(failure_other==1)
stcox i.drug

*sensi-all_cause*
stset end_date_allcause ,  origin(start_date) failure(failure_allcause==1)
*additionally adjusting for days between test positive and treatment initiation, and days/months between last vaccination date and treatment initiation; *
stcox i.drug age i.sex i.d_postest_treat_missing i.month_after_vaccinate_missing, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra i.d_postest_treat_missing i.month_after_vaccinate_missing, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* i.d_postest_treat_missing i.month_after_vaccinate_missing, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease i.d_postest_treat_missing i.month_after_vaccinate_missing, strata(region_nhs)
*exclude patients with HIV/AIDS due to small number*
stcox i.drug age i.sex if hiv_aids==0, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if hiv_aids==0, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if hiv_aids==0, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if hiv_aids==0, strata(region_nhs)
*exclude 1-year drug interactions*
stcox i.drug age i.sex if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25))&(drugs_consider_risk>start_date|drugs_consider_risk<(start_date-365.25)), strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25))&(drugs_consider_risk>start_date|drugs_consider_risk<(start_date-365.25)), strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25))&(drugs_consider_risk>start_date|drugs_consider_risk<(start_date-365.25)), strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25))&(drugs_consider_risk>start_date|drugs_consider_risk<(start_date-365.25)), strata(region_nhs)
*excluding patients with treatment records of both sotrovimab and Pax, or with treatment records of any other therapies*
stcox i.drug age i.sex if (sotrovimab_covid_therapeutics==.|paxlovid_covid_therapeutics==.|sotrovimab_covid_therapeutics>start_date_29|paxlovid_covid_therapeutics>start_date_29)&molnupiravir_covid_therapeutics>start_date_29&remdesivir_covid_therapeutics>start_date_29&casirivimab_covid_therapeutics>start_date_29, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if (sotrovimab_covid_therapeutics==.|paxlovid_covid_therapeutics==.|sotrovimab_covid_therapeutics>start_date_29|paxlovid_covid_therapeutics>start_date_29)&molnupiravir_covid_therapeutics>start_date_29&remdesivir_covid_therapeutics>start_date_29&casirivimab_covid_therapeutics>start_date_29, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if (sotrovimab_covid_therapeutics==.|paxlovid_covid_therapeutics==.|sotrovimab_covid_therapeutics>start_date_29|paxlovid_covid_therapeutics>start_date_29)&molnupiravir_covid_therapeutics>start_date_29&remdesivir_covid_therapeutics>start_date_29&casirivimab_covid_therapeutics>start_date_29, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if (sotrovimab_covid_therapeutics==.|paxlovid_covid_therapeutics==.|sotrovimab_covid_therapeutics>start_date_29|paxlovid_covid_therapeutics>start_date_29)&molnupiravir_covid_therapeutics>start_date_29&remdesivir_covid_therapeutics>start_date_29&casirivimab_covid_therapeutics>start_date_29, strata(region_nhs)
*excluding patients who were identified to be pregnant at treatment initiation*
*stcox i.drug age i.sex if pregnancy!=1, strata(region_nhs)
*stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro if pregnancy!=1, strata(region_nhs)
*stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b6.ethnicity_with_missing b5.imd_with_missing i.vaccination_status i.week_after_campaign if pregnancy!=1, strata(region_nhs)
*stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b6.ethnicity_with_missing b5.imd_with_missing i.vaccination_status i.week_after_campaign b1.bmi_g4_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if pregnancy!=1, strata(region_nhs)
*additionally adjusting for rural-urban classification, other comorbidities (dementia, autism, learning disabilities, severe mental illness), and care home residency and housebound status *
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease i.rural_urban_with_missing autism_nhsd care_home_primis dementia_nhsd housebound_opensafely learning_disability_primis serious_mental_illness_nhsd, strata(region_nhs)
*excluding patients who did not have a positive SARS-CoV-2 test record before treatment or initiated treatment after 5 days since positive SARS-CoV-2 test*
stcox i.drug age i.sex if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
*create a 1-day lag in the follow-up start date *
stcox i.drug age i.sex if _t>=2, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if _t>=2, strata(region_nhs) 
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if _t>=2, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if _t>=2, strata(region_nhs)
*create a 2-day lag in the follow-up start date *
stcox i.drug age i.sex if _t>=3, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if _t>=3, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if _t>=3, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if _t>=3, strata(region_nhs)
*use older version of codelist for solid_cancer and immunosupression*
stcox i.drug age i.sex downs_syndrome solid_cancer haema_disease   imid immunosupression   rare_neuro drugs_consider_risk_contra , strata(region_nhs) 
stcox i.drug age i.sex downs_syndrome solid_cancer haema_disease   imid immunosupression   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* , strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer haema_disease   imid immunosupression   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease , strata(region_nhs)
*using Cox models with calendar date as the underlying time scale*
stset end_date_allcause ,  origin(campaign_start_date) enter(start_date) failure(failure_allcause==1)
stcox i.drug age i.sex, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)

*sensi-2m_covid*
stset end_date_2m if start_date<=mdy(9,1,2022),  origin(start_date) failure(failure_2m==1)
*additionally adjusting for days between test positive and treatment initiation, and days/months between last vaccination date and treatment initiation; *
stcox i.drug age i.sex i.d_postest_treat_missing i.month_after_vaccinate_missing, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra i.d_postest_treat_missing i.month_after_vaccinate_missing, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* i.d_postest_treat_missing i.month_after_vaccinate_missing, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease i.d_postest_treat_missing i.month_after_vaccinate_missing, strata(region_nhs)
*exclude patients with HIV/AIDS due to small number*
stcox i.drug age i.sex if hiv_aids==0, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if hiv_aids==0, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if hiv_aids==0, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if hiv_aids==0, strata(region_nhs)
*exclude 1-year drug interactions*
stcox i.drug age i.sex if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25))&(drugs_consider_risk>start_date|drugs_consider_risk<(start_date-365.25)), strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25))&(drugs_consider_risk>start_date|drugs_consider_risk<(start_date-365.25)), strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25))&(drugs_consider_risk>start_date|drugs_consider_risk<(start_date-365.25)), strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25))&(drugs_consider_risk>start_date|drugs_consider_risk<(start_date-365.25)), strata(region_nhs)
*excluding patients with treatment records of both sotrovimab and Pax, or with treatment records of any other therapies*
stcox i.drug age i.sex if (sotrovimab_covid_therapeutics==.|paxlovid_covid_therapeutics==.|sotrovimab_covid_therapeutics>start_date_29|paxlovid_covid_therapeutics>start_date_29)&molnupiravir_covid_therapeutics>start_date_29&remdesivir_covid_therapeutics>start_date_29&casirivimab_covid_therapeutics>start_date_29, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if (sotrovimab_covid_therapeutics==.|paxlovid_covid_therapeutics==.|sotrovimab_covid_therapeutics>start_date_29|paxlovid_covid_therapeutics>start_date_29)&molnupiravir_covid_therapeutics>start_date_29&remdesivir_covid_therapeutics>start_date_29&casirivimab_covid_therapeutics>start_date_29, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if (sotrovimab_covid_therapeutics==.|paxlovid_covid_therapeutics==.|sotrovimab_covid_therapeutics>start_date_29|paxlovid_covid_therapeutics>start_date_29)&molnupiravir_covid_therapeutics>start_date_29&remdesivir_covid_therapeutics>start_date_29&casirivimab_covid_therapeutics>start_date_29, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if (sotrovimab_covid_therapeutics==.|paxlovid_covid_therapeutics==.|sotrovimab_covid_therapeutics>start_date_29|paxlovid_covid_therapeutics>start_date_29)&molnupiravir_covid_therapeutics>start_date_29&remdesivir_covid_therapeutics>start_date_29&casirivimab_covid_therapeutics>start_date_29, strata(region_nhs)
*excluding patients who were identified to be pregnant at treatment initiation*
*stcox i.drug age i.sex if pregnancy!=1, strata(region_nhs)
*stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro if pregnancy!=1, strata(region_nhs)
*stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b6.ethnicity_with_missing b5.imd_with_missing i.vaccination_status i.week_after_campaign if pregnancy!=1, strata(region_nhs)
*stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b6.ethnicity_with_missing b5.imd_with_missing i.vaccination_status i.week_after_campaign b1.bmi_g4_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if pregnancy!=1, strata(region_nhs)
*additionally adjusting for rural-urban classification, other comorbidities (dementia, autism, learning disabilities, severe mental illness), and care home residency and housebound status *
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease i.rural_urban_with_missing autism_nhsd care_home_primis dementia_nhsd housebound_opensafely learning_disability_primis serious_mental_illness_nhsd, strata(region_nhs)
*excluding patients who did not have a positive SARS-CoV-2 test record before treatment or initiated treatment after 5 days since positive SARS-CoV-2 test*
stcox i.drug age i.sex if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
*create a 1-day lag in the follow-up start date *
stcox i.drug age i.sex if _t>=2, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if _t>=2, strata(region_nhs) 
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if _t>=2, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if _t>=2, strata(region_nhs)
*create a 2-day lag in the follow-up start date *
stcox i.drug age i.sex if _t>=3, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if _t>=3, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if _t>=3, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if _t>=3, strata(region_nhs)
*use older version of codelist for solid_cancer and immunosupression*
stcox i.drug age i.sex downs_syndrome solid_cancer haema_disease   imid immunosupression   rare_neuro drugs_consider_risk_contra , strata(region_nhs) 
stcox i.drug age i.sex downs_syndrome solid_cancer haema_disease   imid immunosupression   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* , strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer haema_disease   imid immunosupression   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease , strata(region_nhs)
*using Cox models with calendar date as the underlying time scale*
stset end_date_2m if start_date<=mdy(9,1,2022),  origin(campaign_start_date) enter(start_date) failure(failure_2m==1)
stcox i.drug age i.sex, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)


*safety outcome*
*death not due to covid*
by drug, sort: count if death_date!=.
by drug, sort: count if death_with_covid_on_the_death_ce!=.
by drug, sort: count if death_with_covid_on_the_death_ce==.&death_date!=.
gen death_without_covid=death_date if death_with_covid_on_the_death_ce==.&death_date!=.


*further sensitivity analyses*
clear
use ./output/sensitivity_update.dta
drop if high_risk_group_new==0
*recode Paxlovid as 1*
replace drug=1-drug
label define drug_Paxlovid2 0 "sotrovimab" 1 "Paxlovid"
label values drug drug_Paxlovid2
*not exclude contra*
stset end_date ,  origin(start_date) failure(failure==1)
stcox drug
tab failure drug,m col
*need to adjust for contra*
gen solid_organ_contra=(solid_organ_new==1|solid_organ_therapeutics==1|solid_organ_transplant_snomed<=start_date)
gen liver_contra=(advanced_decompensated_cirrhosis<=start_date|decompensated_cirrhosis_icd10<=start_date|ascitic_drainage_snomed<=start_date|liver_disease_nhsd_icd10<=start_date)
gen renal_contra=1 if renal_disease==1|renal_therapeutics==1|ckd_stages_3_5<=start_date|ckd_primis_stage==3|ckd_primis_stage==4|ckd_primis_stage==5|ckd3_icd10<=start_date|ckd4_icd10<=start_date|ckd5_icd10<=start_date
replace renal_contra=1 if kidney_transplant<=start_date|kidney_transplant_icd10<=start_date|kidney_transplant_procedure<=start_date|dialysis<=start_date|dialysis_icd10<=start_date|dialysis_procedure<=start_date
replace renal_contra=1 if (egfr_creatinine_ctv3<60&creatinine_operator_ctv3!="<")|(egfr_creatinine_snomed<60&creatinine_operator_snomed!="<")|(eGFR_record<60&eGFR_record>0&eGFR_operator!=">"&eGFR_operator!=">=")|(eGFR_short_record<60&eGFR_short_record>0&eGFR_short_operator!=">"&eGFR_short_operator!=">=")
replace renal_contra=0 if renal_contra==.
mkspline calendar_day_spline = day_after_campaign, cubic nknots(4)
stcox i.drug age i.sex, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro solid_organ_contra liver_contra renal_contra drugs_do_not_use_contra drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*  solid_organ_contra liver_contra renal_contra drugs_do_not_use_contra drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease  solid_organ_contra liver_contra renal_contra drugs_do_not_use_contra drugs_consider_risk_contra, strata(region_nhs)
*only use the drugs_do_not_use list*
drop if solid_organ_new==1|solid_organ_therapeutics==1|solid_organ_transplant_snomed<=start_date
drop if advanced_decompensated_cirrhosis<=start_date|decompensated_cirrhosis_icd10<=start_date|ascitic_drainage_snomed<=start_date|liver_disease_nhsd_icd10<=start_date
drop if renal_disease==1|renal_therapeutics==1|ckd_stages_3_5<=start_date|ckd_primis_stage==3|ckd_primis_stage==4|ckd_primis_stage==5|ckd3_icd10<=start_date|ckd4_icd10<=start_date|ckd5_icd10<=start_date
drop if kidney_transplant<=start_date|kidney_transplant_icd10<=start_date|kidney_transplant_procedure<=start_date
drop if dialysis<=start_date|dialysis_icd10<=start_date|dialysis_procedure<=start_date
drop if (egfr_creatinine_ctv3<60&creatinine_operator_ctv3!="<")|(egfr_creatinine_snomed<60&creatinine_operator_snomed!="<")|(eGFR_record<60&eGFR_record>0&eGFR_operator!=">"&eGFR_operator!=">=")|(eGFR_short_record<60&eGFR_short_record>0&eGFR_short_operator!=">"&eGFR_short_operator!=">=")
drop calendar_day_spline*
mkspline calendar_day_spline = day_after_campaign if drugs_do_not_use>start_date|drugs_do_not_use<(start_date-90), cubic nknots(4)
stset end_date ,  origin(start_date) failure(failure==1)
stcox i.drug age i.sex drugs_consider_risk_contra if drugs_do_not_use>start_date|drugs_do_not_use<(start_date-90), strata(region_nhs)
stcox i.drug age i.sex drugs_consider_risk_contra downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro if drugs_do_not_use>start_date|drugs_do_not_use<(start_date-90), strata(region_nhs)
stcox i.drug age i.sex drugs_consider_risk_contra downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* if drugs_do_not_use>start_date|drugs_do_not_use<(start_date-90), strata(region_nhs)
stcox i.drug age i.sex drugs_consider_risk_contra downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if drugs_do_not_use>start_date|drugs_do_not_use<(start_date-90), strata(region_nhs)
*set 3month cut-off for the two drug lists*
drop if drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-90)
drop if drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-90)
drop calendar_day_spline*
mkspline calendar_day_spline = day_after_campaign, cubic nknots(4)
stset end_date ,  origin(start_date) failure(failure==1)
stcox i.drug age i.sex, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline*, strata(region_nhs)
stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro b1.White_with_missing b5.imd_with_missing i.vaccination_status calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)


log close
