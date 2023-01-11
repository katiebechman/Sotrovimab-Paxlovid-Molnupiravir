********************************************************************************
*
*	Do-file:		multiple_imputation.do
*
*	Project:		sotrovimab-and-Paxlovid
*
*	Programmed by:	Bang Zheng
*
*	Data used:		output/main.dta
*
*	Output:	        logs/MI.log  
*
********************************************************************************
*
*	Purpose: This do-file implements stratified Cox regression after multiple imputation for covariates.
*  
********************************************************************************

* Open a log file
cap log close
log using ./logs/MI, replace t
clear

use ./output/main.dta

stset end_date ,  origin(start_date) failure(failure==1)
keep if _st==1

*MI*
*install ice package by changing ado filepath*
sysdir
sysdir set PLUS "analysis/ado"
sysdir set PERSONAL "analysis/ado"

set seed 1000

ice m.White m.bmi_g3  m.imd   drug age i.sex i.region_nhs downs_syndrome solid_cancer_new haema_disease imid immunosupression_new  rare_neuro i.vaccination_status calendar_day_spline* diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease failure, m(5) saving(imputed,replace)  
clear
use imputed
mi import ice, imputed(White bmi_g3  imd)
mi stset end_date ,  origin(start_date) failure(failure==1)
mi estimate, hr: stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease imid immunosupression_new  rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline*,strata(region_nhs)
mi estimate, hr: stcox i.drug age i.sex downs_syndrome solid_cancer_new haema_disease imid immunosupression_new  rare_neuro b1.White b5.imd i.vaccination_status calendar_day_spline* b1.bmi_g3 diabetes chronic_cardiac_disease hypertension chronic_respiratory_disease, strata(region_nhs)


log close
