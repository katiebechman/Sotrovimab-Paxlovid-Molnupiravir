********************************************************************************
*
*	Do-file:		cox.do
*
*	Project:		sotrovimab-and-Paxlovid
*
*	Programmed by:	Bang Zheng
*
*	Data used:		output/main_mol.dta
*
*	Output:	        logs/cox_mol.log  output/mol/phtest.svg  output/mol/phtest_psw.svg
*
********************************************************************************
*
*	Purpose: This do-file implements stratified Cox regression, propensity score
*   weighted Cox, and subgroup analyses.
*  
********************************************************************************

* Open a log file
cap log close
log using ./logs/cox_mol_update, replace t
clear

use ./output/main_mol_update.dta

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
tab _t drug if failure==0&_t<28&end_date==covid_hosp_date_day_cases,m col
tab _t drug if failure==0&_t<28&end_date==min(sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)&drug==1,m col
tab _t drug if failure==0&_t<28&end_date==min(sotrovimab_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)&drug==0,m col


*un-stratified Cox, with covariate adjustment, complete case*
stcox i.drug
stcox i.drug age i.sex
stcox i.drug age i.sex i.region_nhs
*region_nhs or region_covid_therapeutics? *
stcox i.drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra
stcox i.drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline* 
stcox i.drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline* b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease
*age: 5-year band*
stcox i.drug b7.age_5y_band i.sex i.region_nhs
stcox i.drug b7.age_5y_band i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra
stcox i.drug b7.age_5y_band i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline*
stcox i.drug b7.age_5y_band i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline* b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease
*age: Restricted cubic spline*
stcox i.drug age_spline* i.sex i.region_nhs
stcox i.drug age_spline* i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra
stcox i.drug age_spline* i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline*
stcox i.drug age_spline* i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline* b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease
*PH test*
estat phtest,de

*un-stratified Cox, missing values as a separate category*
stcox i.drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*
stcox i.drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease
*age: 5-year band*
stcox i.drug b7.age_5y_band i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*
stcox i.drug b7.age_5y_band i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease
*age: Restricted cubic spline*
stcox i.drug age_spline* i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* 
stcox i.drug age_spline* i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease
estat phtest,de

*stratified Cox, complete case*
stcox i.drug age i.sex, strata(region_nhs)
*stcox i.drug age i.sex, strata(region_nhs week_after_campaign)
*too few events to allow two-level stratification*
stcox i.drug age i.sex i.region_nhs, strata(week_after_campaign)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline*, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline* b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(week_after_campaign)
stcox i.drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3, strata(week_after_campaign)
stcox i.drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(week_after_campaign)
*age: 5-year band*
stcox i.drug b7.age_5y_band i.sex, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex i.region_nhs, strata(week_after_campaign)
stcox i.drug b7.age_5y_band i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline*, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline* b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(week_after_campaign)
stcox i.drug b7.age_5y_band i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3, strata(week_after_campaign)
stcox i.drug b7.age_5y_band i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(week_after_campaign)
*age: Restricted cubic spline*
stcox i.drug age_spline* i.sex, strata(region_nhs)
stcox i.drug age_spline* i.sex i.region_nhs, strata(week_after_campaign)
stcox i.drug age_spline* i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug age_spline* i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline*, strata(region_nhs)
stcox i.drug age_spline* i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline* b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
estat phtest,de
stcox i.drug age_spline* i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(week_after_campaign)
stcox i.drug age_spline* i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3, strata(week_after_campaign)
stcox i.drug age_spline* i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(week_after_campaign)
estat phtest,de

*stratified Cox, missing values as a separate category*
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3, strata(week_after_campaign)
stcox i.drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(week_after_campaign)
*age: 5-year band*
stcox i.drug b7.age_5y_band i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug b7.age_5y_band i.region_nhs i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3, strata(week_after_campaign)
stcox i.drug b7.age_5y_band i.region_nhs i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(week_after_campaign)
*age: Restricted cubic spline*
stcox i.drug age_spline* i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*, strata(region_nhs)
stcox i.drug age_spline* i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* i.b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
estat phtest,de
stcox i.drug age_spline* i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3, strata(week_after_campaign)
stcox i.drug age_spline* i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(week_after_campaign)
estat phtest,de



*propensity score weighted Cox*
do "analysis/ado/psmatch2.ado"
*age continuous, complete case*
psmatch2 drug age i.sex i.region_nhs, logit
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug age i.sex i.region_nhs ) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug

psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug

psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline*, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline*) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug

psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline* b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline* b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug
estat phtest,de

*age continuous, missing values as a separate categorye*
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug

psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug
estat phtest,de

*age: 5-year band, complete case*
psmatch2 drug b7.age_5y_band i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline*, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de 
teffects ipw (failure) (drug b7.age_5y_band i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline*) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug

psmatch2 drug b7.age_5y_band i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline* b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug b7.age_5y_band i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White b5.imd i.vaccination_g3 calendar_day_spline* b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug
estat phtest,de

*age: 5-year band, missing values as a separate categorye*
psmatch2 drug b7.age_5y_band i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug b7.age_5y_band i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug

psmatch2 drug b7.age_5y_band i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug b7.age_5y_band i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug
estat phtest,de

*age spline, missing values as a separate categorye*
psmatch2 drug age_spline* i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug age_spline* i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug

psmatch2 drug age_spline* i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure) (drug age_spline* i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease) if _pscore!=.
tebalance summarize
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug
estat phtest,de


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
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age_spline* i.sex, strata(region_nhs)
stcox i.drug age_spline* i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug age_spline* i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*, strata(region_nhs)
stcox i.drug age_spline* i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
*PSW age continuous, missing values as a separate categorye*
psmatch2 drug age i.sex i.region_nhs  if _st==1, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure_allcause) (drug age i.sex i.region_nhs ) if _pscore!=.
tebalance summarize
stset end_date_allcause [pwei=psweight],  origin(start_date) failure(failure_allcause==1)
stcox i.drug
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra  if _st==1, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure_allcause) (drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra ) if _pscore!=.
tebalance summarize
stset end_date_allcause [pwei=psweight],  origin(start_date) failure(failure_allcause==1)
stcox i.drug
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*  if _st==1, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure_allcause) (drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* ) if _pscore!=.
tebalance summarize
stset end_date_allcause [pwei=psweight],  origin(start_date) failure(failure_allcause==1)
stcox i.drug
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if _st==1, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure_allcause) (drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease) if _pscore!=.
tebalance summarize
stset end_date_allcause [pwei=psweight],  origin(start_date) failure(failure_allcause==1)
stcox i.drug
*PSW age spline, missing values as a separate categorye*
psmatch2 drug age_spline* i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if _st==1, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure_allcause) (drug age_spline* i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease) if _pscore!=.
tebalance summarize
stset end_date_allcause [pwei=psweight],  origin(start_date) failure(failure_allcause==1)
stcox i.drug


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
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*, strata(region_nhs)
stcox i.drug b7.age_5y_band i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age_spline* i.sex, strata(region_nhs)
stcox i.drug age_spline* i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra, strata(region_nhs)
stcox i.drug age_spline* i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*, strata(region_nhs)
stcox i.drug age_spline* i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
*PSW age continuous, missing values as a separate categorye*
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if _st==1, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure_2m) (drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease) if _pscore!=.
tebalance summarize
stset end_date_2m [pwei=psweight],  origin(start_date) failure(failure_2m==1)
stcox i.drug
*PSW age spline, missing values as a separate categorye*
psmatch2 drug age_spline* i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if _st==1, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
teffects ipw (failure_2m) (drug age_spline* i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease) if _pscore!=.
tebalance summarize
stset end_date_2m [pwei=psweight],  origin(start_date) failure(failure_2m==1)
stcox i.drug




*sensitivity analysis*
stset end_date ,  origin(start_date) failure(failure==1)
*additionally adjusting for days between test positive and treatment initiation, and days/months between last vaccination date and treatment initiation; *
stcox i.drug age i.sex i.d_postest_treat_missing i.month_after_vaccinate_missing, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra i.d_postest_treat_missing i.month_after_vaccinate_missing, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* i.d_postest_treat_missing i.month_after_vaccinate_missing, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease i.d_postest_treat_missing i.month_after_vaccinate_missing, strata(region_nhs)
*excluding patients with treatment records of both mol and Pax, or with treatment records of any other therapies*
stcox i.drug age i.sex if (molnupiravir_covid_therapeutics==.|paxlovid_covid_therapeutics==.|molnupiravir_covid_therapeutics>start_date_29|paxlovid_covid_therapeutics>start_date_29)&sotrovimab_covid_therapeutics>start_date_29&remdesivir_covid_therapeutics>start_date_29&casirivimab_covid_therapeutics>start_date_29, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if (molnupiravir_covid_therapeutics==.|paxlovid_covid_therapeutics==.|molnupiravir_covid_therapeutics>start_date_29|paxlovid_covid_therapeutics>start_date_29)&sotrovimab_covid_therapeutics>start_date_29&remdesivir_covid_therapeutics>start_date_29&casirivimab_covid_therapeutics>start_date_29, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* if (molnupiravir_covid_therapeutics==.|paxlovid_covid_therapeutics==.|molnupiravir_covid_therapeutics>start_date_29|paxlovid_covid_therapeutics>start_date_29)&sotrovimab_covid_therapeutics>start_date_29&remdesivir_covid_therapeutics>start_date_29&casirivimab_covid_therapeutics>start_date_29, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if (molnupiravir_covid_therapeutics==.|paxlovid_covid_therapeutics==.|molnupiravir_covid_therapeutics>start_date_29|paxlovid_covid_therapeutics>start_date_29)&sotrovimab_covid_therapeutics>start_date_29&remdesivir_covid_therapeutics>start_date_29&casirivimab_covid_therapeutics>start_date_29, strata(region_nhs)
*excluding patients who were identified to be pregnant at treatment initiation*
*stcox i.drug age i.sex if pregnancy!=1, strata(region_nhs)
*stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if pregnancy!=1, strata(region_nhs)
*stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* if pregnancy!=1, strata(region_nhs)
*stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if pregnancy!=1, strata(region_nhs)
*additionally adjusting for rural-urban classification, other comorbidities (dementia, autism, learning disabilities, severe mental illness), and care home residency and housebound status *
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease i.rural_urban_with_missing autism_nhsd care_home_primis dementia_nhsd housebound_opensafely learning_disability_primis serious_mental_illness_nhsd, strata(region_nhs)
*excluding patients who did not have a positive SARS-CoV-2 test record before treatment or initiated treatment after 5 days since positive SARS-CoV-2 test*
tab failure drug if d_postest_treat>=0&d_postest_treat<=5,m col
stcox i.drug age i.sex if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if d_postest_treat>=0&d_postest_treat<=5, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug
*create a 1-day lag in the follow-up start date *
stset end_date ,  origin(start_date) failure(failure==1)
tab failure drug if _t>=2,m col
stcox i.drug age i.sex if _t>=2, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if _t>=2, strata(region_nhs) 
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* if _t>=2, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if _t>=2, strata(region_nhs)
*create a 2-day lag in the follow-up start date *
tab failure drug if _t>=3,m col
stcox i.drug age i.sex if _t>=3, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if _t>=3, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* if _t>=3, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if _t>=3, strata(region_nhs)
*use older version of codelist for solid_cancer and immunosupression*
stcox i.drug age i.sex , strata(region_nhs)
stcox i.drug age i.sex  solid_cancer haema_disease   imid immunosupression   rare_neuro drugs_consider_risk_contra , strata(region_nhs) 
stcox i.drug age i.sex  solid_cancer haema_disease   imid immunosupression   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* , strata(region_nhs)
stcox i.drug age i.sex  solid_cancer haema_disease   imid immunosupression   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease , strata(region_nhs)
*stratify by STP*
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*, strata(stp)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(stp)
psmatch2 drug age i.sex i.stp  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug
*exclude missing high_risk_group_new*
stset end_date ,  origin(start_date) failure(failure==1)
tab failure drug if high_risk_group_new==1,m col
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if high_risk_group_new==1, strata(region_nhs)
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if high_risk_group_new==1, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug
*exclude 1-year drug interactions*
stset end_date ,  origin(start_date) failure(failure==1)
tab failure drug if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25)),m col
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25)), strata(region_nhs)
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25)), logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug

stset end_date ,  origin(start_date) failure(failure==1)
tab failure drug if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25))&(drugs_consider_risk>start_date|drugs_consider_risk<(start_date-365.25)),m col
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25))&(drugs_consider_risk>start_date|drugs_consider_risk<(start_date-365.25)), strata(region_nhs)
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25))&(drugs_consider_risk>start_date|drugs_consider_risk<(start_date-365.25)), logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug
*ATT weight*
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, logit
drop psweight
gen psweight=cond( drug ==1,1,_pscore/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug
estat phtest,de
*exclude all patients in the non-overlapping parts of the PS distribution*
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease , logit common
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
sum _pscore if drug==0,de
gen _pscore_mol_min=r(min)
gen _pscore_mol_max=r(max)
sum _pscore if drug==1,de
gen _pscore_pax_min=r(min)
gen _pscore_pax_max=r(max)
stset end_date if (drug==0&_pscore>=_pscore_pax_min&_pscore<=_pscore_pax_max)|(drug==1&_pscore>=_pscore_mol_min&_pscore<=_pscore_mol_max) [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug
drop _pscore_mol_min _pscore_mol_max _pscore_pax_min _pscore_pax_max
*ATE additionally adjust for region*
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease , logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
stset end_date [pwei=psweight],  origin(start_date) failure(failure==1)
stcox i.drug i.region_nhs


*sensitivity analysis-all cause*
stset end_date_allcause ,  origin(start_date) failure(failure_allcause==1)
*additionally adjusting for days between test positive and treatment initiation, and days/months between last vaccination date and treatment initiation; *
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease i.d_postest_treat_missing i.month_after_vaccinate_missing, strata(region_nhs)
*excluding patients with treatment records of both mol and Pax, or with treatment records of any other therapies*
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if (molnupiravir_covid_therapeutics==.|paxlovid_covid_therapeutics==.|molnupiravir_covid_therapeutics>start_date_29|paxlovid_covid_therapeutics>start_date_29)&sotrovimab_covid_therapeutics>start_date_29&remdesivir_covid_therapeutics>start_date_29&casirivimab_covid_therapeutics>start_date_29, strata(region_nhs)
*excluding patients who were identified to be pregnant at treatment initiation*
*stcox i.drug age i.sex if pregnancy!=1, strata(region_nhs)
*stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if pregnancy!=1, strata(region_nhs)
*stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* if pregnancy!=1, strata(region_nhs)
*stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if pregnancy!=1, strata(region_nhs)
*additionally adjusting for rural-urban classification, other comorbidities (dementia, autism, learning disabilities, severe mental illness), and care home residency and housebound status *
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease i.rural_urban_with_missing autism_nhsd care_home_primis dementia_nhsd housebound_opensafely learning_disability_primis serious_mental_illness_nhsd, strata(region_nhs)
*excluding patients who did not have a positive SARS-CoV-2 test record before treatment or initiated treatment after 5 days since positive SARS-CoV-2 test*
tab failure_allcause drug if d_postest_treat>=0&d_postest_treat<=5,m col
stcox i.drug age i.sex if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if d_postest_treat>=0&d_postest_treat<=5, strata(region_nhs)
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if d_postest_treat>=0&d_postest_treat<=5, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
stset end_date_allcause [pwei=psweight],  origin(start_date) failure(failure_allcause==1)
stcox i.drug
*create a 1-day lag in the follow-up start date *
stset end_date_allcause ,  origin(start_date) failure(failure_allcause==1)
tab failure_allcause drug if _t>=2,m col
stcox i.drug age i.sex if _t>=2, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if _t>=2, strata(region_nhs) 
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* if _t>=2, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if _t>=2, strata(region_nhs)
*create a 2-day lag in the follow-up start date *
tab failure_allcause drug if _t>=3,m col
stcox i.drug age i.sex if _t>=3, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra if _t>=3, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* if _t>=3, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if _t>=3, strata(region_nhs)
*use older version of codelist for solid_cancer and immunosupression*
stcox i.drug age i.sex , strata(region_nhs)
stcox i.drug age i.sex  solid_cancer haema_disease   imid immunosupression   rare_neuro drugs_consider_risk_contra , strata(region_nhs) 
stcox i.drug age i.sex  solid_cancer haema_disease   imid immunosupression   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* , strata(region_nhs)
stcox i.drug age i.sex  solid_cancer haema_disease   imid immunosupression   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease , strata(region_nhs)
*stratify by STP*
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*, strata(stp)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(stp)
psmatch2 drug age i.sex i.stp  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
stset end_date_allcause [pwei=psweight],  origin(start_date) failure(failure_allcause==1)
stcox i.drug
*exclude missing high_risk_group_new*
stset end_date_allcause ,  origin(start_date) failure(failure_allcause==1)
tab failure_allcause drug if high_risk_group_new==1,m col
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if high_risk_group_new==1, strata(region_nhs)
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if high_risk_group_new==1, logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
stset end_date_allcause [pwei=psweight],  origin(start_date) failure(failure_allcause==1)
stcox i.drug
*exclude 1-year drug interactions*
stset end_date_allcause ,  origin(start_date) failure(failure_allcause==1)
tab failure_allcause drug if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25)),m col
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25)), strata(region_nhs)
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25)), logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
stset end_date_allcause [pwei=psweight],  origin(start_date) failure(failure_allcause==1)
stcox i.drug

stset end_date_allcause ,  origin(start_date) failure(failure_allcause==1)
tab failure_allcause drug if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25))&(drugs_consider_risk>start_date|drugs_consider_risk<(start_date-365.25)),m col
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25))&(drugs_consider_risk>start_date|drugs_consider_risk<(start_date-365.25)), strata(region_nhs)
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if (drugs_do_not_use>start_date|drugs_do_not_use<(start_date-365.25))&(drugs_consider_risk>start_date|drugs_consider_risk<(start_date-365.25)), logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
stset end_date_allcause [pwei=psweight],  origin(start_date) failure(failure_allcause==1)
stcox i.drug
*ATT weight*
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, logit
drop psweight
gen psweight=cond( drug ==1,1,_pscore/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
stset end_date_allcause [pwei=psweight],  origin(start_date) failure(failure_allcause==1)
stcox i.drug
estat phtest,de
*exclude all patients in the non-overlapping parts of the PS distribution*
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease , logit common
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
sum _pscore if drug==0,de
gen _pscore_mol_min=r(min)
gen _pscore_mol_max=r(max)
sum _pscore if drug==1,de
gen _pscore_pax_min=r(min)
gen _pscore_pax_max=r(max)
stset end_date_allcause if (drug==0&_pscore>=_pscore_pax_min&_pscore<=_pscore_pax_max)|(drug==1&_pscore>=_pscore_mol_min&_pscore<=_pscore_mol_max) [pwei=psweight],  origin(start_date) failure(failure_allcause==1)
stcox i.drug
drop _pscore_mol_min _pscore_mol_max _pscore_pax_min _pscore_pax_max
*ATE additionally adjust for region*
psmatch2 drug age i.sex i.region_nhs  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease , logit
drop psweight
gen psweight=cond( drug ==1,1/_pscore,1/(1-_pscore)) if _pscore!=.
sum psweight,de
by drug, sort: sum _pscore ,de
stset end_date_allcause [pwei=psweight],  origin(start_date) failure(failure_allcause==1)
stcox i.drug i.region_nhs



*subgroup analysis*
stset end_date ,  origin(start_date) failure(failure==1)
stcox i.drug##i.sex age  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if sex==0, strata(region_nhs)
stcox i.drug age  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if sex==1, strata(region_nhs)

stcox i.drug##i.age_group3 i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)

stcox i.drug##i.age_50 i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)

stcox i.drug##i.age_55 i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if age_55==0, strata(region_nhs)
stcox i.drug i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if age_55==1, strata(region_nhs)

stcox i.drug##i.age_60 i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if age_60==0, strata(region_nhs)
stcox i.drug i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if age_60==1, strata(region_nhs)

stcox i.drug##i.White age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if White==1, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if White==0, strata(region_nhs)

stcox i.drug##i.solid_cancer_new age i.sex   haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease , strata(region_nhs)
stcox i.drug age i.sex   haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if solid_cancer_new==1, strata(region_nhs)
stcox i.drug age i.sex   haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if solid_cancer_new==0, strata(region_nhs)
stcox i.drug##i.haema_disease age i.sex  solid_cancer_new    imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new    imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if haema_disease==1, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new    imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if haema_disease==0, strata(region_nhs)
stcox i.drug##i.imid age i.sex  solid_cancer_new haema_disease    immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease , strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease    immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if imid==1, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease    immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if imid==0, strata(region_nhs)
stcox i.drug##i.immunosupression_new age i.sex  solid_cancer_new haema_disease   imid   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease , strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if immunosupression_new==1, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if immunosupression_new==0, strata(region_nhs)
stcox i.drug##i.rare_neuro age i.sex  solid_cancer_new haema_disease   imid immunosupression_new drugs_consider_risk_contra   b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease , strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new  drugs_consider_risk_contra  b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if rare_neuro==0, strata(region_nhs)

stcox i.drug##i.bmi_g3 age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*  diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*  diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if bmi_g3==1, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if bmi_g3==2, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if bmi_g3==3, strata(region_nhs)

stcox i.drug##i.bmi_25 age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*  diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*  diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if bmi_25==0, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if bmi_25==1, strata(region_nhs)

stcox i.drug##i.bmi_30 age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*  diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline*  diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if bmi_30==0, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if bmi_30==1, strata(region_nhs)

stcox i.drug##i.diabetes age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing  chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing chronic_cardiac_disease hypertension chronic_respiratory_disease if diabetes ==0, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing chronic_cardiac_disease hypertension chronic_respiratory_disease if diabetes ==1, strata(region_nhs)

stcox i.drug##i.chronic_cardiac_disease age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes  hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes  hypertension chronic_respiratory_disease if chronic_cardiac_disease==0, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes  hypertension chronic_respiratory_disease if chronic_cardiac_disease==1, strata(region_nhs)

stcox i.drug##i.hypertension age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease  chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease  chronic_respiratory_disease if hypertension==0, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease  chronic_respiratory_disease if hypertension==1, strata(region_nhs)

stcox i.drug##i.chronic_respiratory_disease age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension , strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension if chronic_respiratory_disease==0, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension if chronic_respiratory_disease==1, strata(region_nhs)

stcox i.drug##i.vaccination_g3 age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing  calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing  calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if vaccination_g3==1, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing  calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if vaccination_g3==2, strata(region_nhs)
*stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing  calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if vaccination_3==0, strata(region_nhs)
*stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing  calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if vaccination_status==0, strata(region_nhs)

stcox i.drug##i.d_postest_treat_g2 age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if d_postest_treat_g2==0, strata(region_nhs)
stcox i.drug age i.sex  solid_cancer_new haema_disease   imid immunosupression_new   rare_neuro drugs_consider_risk_contra b1.White_with_missing b5.imd_with_missing i.vaccination_g3 calendar_day_spline* b1.bmi_g3_with_missing diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease if d_postest_treat_g2==1, strata(region_nhs)

*use minimal-adjusted model*
stcox i.drug##i.sex age , strata(region_nhs)
stcox i.drug age if sex==0, strata(region_nhs)
stcox i.drug age if sex==1, strata(region_nhs)

stcox i.drug##i.age_group3 i.sex , strata(region_nhs)
stcox i.drug i.sex  if age_group3==0, strata(region_nhs)
stcox i.drug i.sex if age_group3==1, strata(region_nhs)
stcox i.drug i.sex if age_group3==2, strata(region_nhs)

stcox i.drug##i.age_50 i.sex , strata(region_nhs)
stcox i.drug i.sex  if age_50==0, strata(region_nhs)
stcox i.drug i.sex  if age_50==1, strata(region_nhs)

stcox i.drug##i.age_55 i.sex , strata(region_nhs)
stcox i.drug i.sex  if age_55==0, strata(region_nhs)
stcox i.drug i.sex if age_55==1, strata(region_nhs)

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
*stcox i.drug age i.sex if vaccination_status==0, strata(region_nhs)

stcox i.drug##i.d_postest_treat_g2 age i.sex , strata(region_nhs)
stcox i.drug age i.sex if d_postest_treat_g2==0, strata(region_nhs)
stcox i.drug age i.sex if d_postest_treat_g2==1, strata(region_nhs)



*safety outcome*
*death not due to covid*
by drug, sort: count if death_date!=.
by drug, sort: count if death_with_covid_on_the_death_ce!=.
by drug, sort: count if death_with_covid_on_the_death_ce==.&death_date!=.
gen death_without_covid=death_date if death_with_covid_on_the_death_ce==.&death_date!=.

log close
