********************************************************************************
*
*	Do-file:			define covariates.do
*	Project:			Sotrovimab-Paxlovid-Molnupiravir
*   Date:  				15/2/23
*	Programmed by:		Katie Bechman
* 	Description:		data management, reformat variables, categorise variables, label variables 
*	Data used:			data in memory (from output/input.csv) 
*	Data created:		analysis files/main.dta  (main analysis dataset)
*	Other output:		logfiles, printed to folder $Logdir
*	User installed ado: (place .ado file(s) in analysis folder)

****************************************************************************************************************
**Set filepaths
global projectdir "C:\Users\k1635179\OneDrive - King's College London\Katie\OpenSAFELY\Sotrovimab-Paxlovid-Molnupiravir"
//global projectdir `c(pwd)'
di "$projectdir"
capture mkdir "$projectdir/output/data"
capture mkdir "$projectdir/output/figures"
capture mkdir "$projectdir/output/tables"
global logdir "$projectdir/logs"
di "$logdir"
* Open a log file
cap log close
log using "$logdir/cleaning_dataset.log", replace

* import dataset
import delimited "$projectdir/output/input.csv", clear

*Set Ado file path
adopath + "$projectdir/analysis/extra_ados"

//describe
//codebook

*  Convert strings to dates  *
foreach var of varlist 	 covid_test_positive_date				///
						 covid_test_positive_date2 				///
						 covid_symptoms_snomed 					///	
						 prior_covid_date				    	///
						 sotrovimab								///
						 molnupiravir							/// 
						 paxlovid								///
						 remdesivir								///
						 casirivimab							///
						 sotrovimab_not_start					///
						 molnupiravir_not_start					///
						 paxlovid_not_start						///
						 sotrovimab_stopped						///
						 molnupiravir_stopped					///
						 paxlovid_stopped						///
						 date_treated							///
						 last_vaccination_date 					///
						 death_date								///
						 dereg_date       						///
						 covid_hosp_discharge					///
						 any_covid_hosp_discharge				///			
						 ae_diverticulitis_icd					///
						 ae_diverticulitis_snomed				///
					     ae_diarrhoea_snomed					///
						 ae_taste_snomed						///
						 ae_taste_icd							///
						 ae_rheumatoid_arthritis				///
						 ae_sle									///
						 ae_psoriasis_snomed					///
						 ae_psoriatic_arthritis_snomed			///
						 ae_ankylosing_spondylitis_ctv			///
						 ae_ibd_snomed							///
						 rheumatoid_arthritis_nhsd_snomed		///
						 rheumatoid_arthritis_nhsd_icd10		///
						 sle_nhsd_ctv							///
						 sle_nhsd_icd10							///
						 psoriasis_nhsd							///
						 psoriatic_arthritis_nhsd				///
						 ankylosing_spondylitis_nhsd			///
						 ibd_ctv								///
						 all_hosp_date 							///
						 covid_hosp_date 						///
						 emerg_covid_hosp_date 					///
						 covid_hosp_date_mabs_procedure 		///
						 died_date_ons							///
						 died_ons_covid							///
						 {					 
	capture confirm string variable `var'
	if _rc==0 {
	rename `var' a
	gen `var' = date(a, "YMD")
	drop a
	format %td `var'
 }
}

*check hosp/death event date range*
codebook covid_test_positive_date all_hosp_date died_date_ons

***Exposure*
*check drug given within covid test + 5 days 
foreach var of varlist sotrovimab molnupiravir paxlovid {
    gen `var'_check = 1 if `var'>=covid_test_positive_date & `var'<=covid_test_positive_date+5 & `var'!=.
	replace `var'_check = 0 if (`var' < covid_test_positive_date | `var' > covid_test_positive_date + 5) & `var'!=.
	codebook `var'
	tab `var'_check, m  // should be no 0s
	sum `var'_not_start // number prescribed but not started
	gen `var'_started = 1 if `var'_not_start==. & `var'_check==1
	tab `var'_started, m
}
gen drug=1 if sotrovimab==date_treated & sotrovimab_start==1
replace drug=2 if paxlovid==date_treated & paxlovid_start==1
replace drug=3 if molnupiravir==date_treated & molnupiravir_start==1
replace drug=0 if drug==.
label define drug 0 "control" 1 "sotrovimab" 2 "paxlovid" 3"molnupiravir", replace
label values drug drug
tab drug, m
gen start_date=date_treated if drug>0
format start_date %td
egen median_time = median(date_treated - covid_test_positive_date) if drug>0
egen median_max = max(median_time)
replace start_date = covid_test_positive_date + median_max if drug==0

***Inclusion criteria*
keep if age>=18 & age<110
keep if sex=="F"|sex=="M"
keep if has_died==0
*check IMID on drug should only include those with imid AND (imid_drug or drug_HCD)* 
tab imid_on_drug if imid_on_drug==1 & imid_nhsd==1 & imid_drug==0 & imid_drug_hcd==0 // should be no 0
keep if eligible==1	// should be no 0s		
*check covid positive, and not repeat covid test after an infection within 30 days prior
tab covid_test_positive covid_positive_previous_30_days, m
keep if covid_test_positive==1 & covid_positive_previous_30_days==0

***Exclusion criteria*
*capture and exclude COVID-hosp admission/death or deregistration is on the start date [not to exclude other causes of admission on start date - drug reactions etc]
bys drug: count if start_date>=dereg_date & start_date !=.
bys drug: count if start_date>death_date
drop if start_date>=death_date | start_date>=dereg_date

*** Primary outcome - AESI
gen imae = 0
gen imae_serious = 0
gen aesi_spc = 0
gen aesi_spc_serious = 0
global aesi  	ae_diverticulitis_icd					///
				ae_diverticulitis_snomed				///
				ae_diarrhoea_snomed						///
				ae_taste_snomed							///
				ae_taste_icd	
global imae  	ae_rheumatoid_arthritis_snomed			///
				ae_rheumatoid_arthritis_icd				///
				ae_sle_ctv								///
				ae_sle_icd								///
				ae_psoriasis_snomed						///
				ae_psoriatic_arthritis_snomed			///
				ae_ankylosing_spondylitis_ctv			///
				ae_ibd_snomed
foreach x in $aesi {
				display "`x'"
				bys drug: count if (`x' < covid_test_positive_date | `x' > covid_test_positive_date + 28) & `x'!=.
				replace `x'=. if (`x' < covid_test_positive_date | `x' > covid_test_positive_date + 28) & `x'!=.
				replace aesi_spc = 1 if `x' !=.
				replace aesi_spc_serious = 1 if `x' !=. & (`x'==all_hosp_date)
}
foreach y in $imae {
				display "`y'"
				bys drug: count if (`y' < covid_test_positive_date | `y' > covid_test_positive_date + 28) & `x'!=.
				replace `y'=. if (`y' < covid_test_positive_date | `y' > covid_test_positive_date + 28) & `y'!=.
				replace imae  = 1 if `y' !=.
				replace imae_serious = 1 if `y' !=. & (`y'==all_hosp_date)
}
gen new_ae_ra = ae_rheumatoid_arthritis if rheumatoid_arthritis_nhsd_snomed==0 & rheumatoid_arthritis_nhsd_icd10==0  // need to do icd for psa, psorasis, ibd, ankspon
gen new_ae_sle = ae_sle if sle_nhsd_ctv==0 & sle_nhsd_icd10==0
gen new_ae_psoriasis = ae_psoriasis_snomed if psoriasis_nhsd==0
gen new_ae_psoriatic_arthritis = ae_psoriatic_arthritis_snomed if psoriatic_arthritis_nhsd==0
gen new_ae_ankylosing_spondylitis = ae_ankylosing_spondylitis_ctv if ankylosing_spondylitis_nhsd==0
gen new_ae_ibd = ae_ae_ibd_snomed if ibd_ctv==0

gen aesi = 0
gen aesi_serious = 0




*** Secondary outcome - SAEs hospitalisation or death including COVID-19
*correcting COVID hosp events:  ignore any day cases or sotro initiators who had COVID hosp record with mab procedure codes [nb have only excluded if same day, not day+1]
by drug, sort: count if covid_hosp_date == covid_hosp_discharge & covid_hosp_date!=. 
by drug, sort: count if covid_hosp_date == covid_hosp_date_mabs_procedure & covid_hosp_date_mabs_procedure!=. 
by drug, sort: count if covid_hosp_date == covid_hosp_date_mabs_procedure & covid_hosp_date!=. & covid_hosp_date==covid_hosp_discharge
replace covid_hosp_date=. if covid_hosp_date==covid_hosp_date_mabs_procedure & covid_hosp_date_mabs_procedure!=. & drug==1



*** Secondary outcome - severe drug reactions (including DRESS, SJS, TEN, anaphylaxis) 



*capture and exclude COVID-hospital admission/death on the start date
by drug, sort: count if start_date != covid_hosp_date & covid_hosp_date !=.
by drug, sort: count if start_date == covid_hosp_date 
drop if start_date>=covid_hospitalisation_outcome_da| start_date>=death_with_covid_on_the_death_ce|start_date>=death_date|start_date>=dereg_date


*time to admission
gen days_to_covid_admission=covid_hosp_date-start_date if covid_hospitalisation_outcome_da!=.
by drug days_to_covid_admission, sort: count if covid_hospitalisation_outcome_da!=.

by drug, sort: count if covid_hosp_outcome_date2==covid_hosp_discharge_date2&covid_hosp_outcome_date2<=(start_date+28)&days_to_covid_admission>=2
by drug days_to_covid_admission, sort: count if covid_hosp_outcome_date2==covid_hosp_discharge_date2&covid_hosp_outcome_date2<=(start_date+28)&days_to_covid_admission>=2
by drug, sort: count if covid_hosp_outcome_date2==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure<=(start_date+28)&days_to_covid_admission>=2
by drug days_to_covid_admission, sort: count if covid_hosp_outcome_date2==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure<=(start_date+28)&days_to_covid_admission>=2
by drug, sort: count if covid_hosp_outcome_date2==covid_hosp_date_mabs_procedure&covid_hosp_outcome_date2<=(start_date+28)&covid_hosp_outcome_date2==covid_hosp_discharge_date2&days_to_covid_admission>=2
by drug, sort: count if covid_hosp_outcome_date2==covid_hosp_date_mabs_procedure&covid_hosp_outcome_date2<=(start_date+28)&(covid_hosp_discharge_date2 - covid_hosp_outcome_date2)==1&days_to_covid_admission>=2
by drug, sort: count if covid_hosp_outcome_date2==covid_hosp_date_mabs_procedure&covid_hosp_outcome_date2<=(start_date+28)&(covid_hosp_discharge_date2 - covid_hosp_outcome_date2)==2&days_to_covid_admission>=2
count if covid_hosp_outcome_date2==covid_hosp_discharge_date2&covid_hosp_outcome_date2<=(start_date+28)&days_to_covid_admission>=2&sotrovimab_covid_therapeutics==.&casirivimab_covid_therapeutics==.
count if covid_hosp_outcome_date2==covid_hosp_discharge_date2&covid_hosp_outcome_date2<=(start_date+28)&days_to_covid_admission>=2&sotrovimab_covid_therapeutics==.&casirivimab_covid_therapeutics==.
by drug, sort: count if covid_hosp_outcome_date2==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure<=(start_date+28)&days_to_covid_admission>=2&sotrovimab_covid_therapeutics==.&casirivimab_covid_therapeutics==.
by drug days_to_covid_admission, sort: count if covid_hosp_outcome_date2==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure<=(start_date+28)&days_to_covid_admission>=2&sotrovimab_covid_therapeutics==.&casirivimab_covid_therapeutics==.
*ignore and censor day cases on or after day 2 from this analysis*
*ignore and censor admissions for mab procedure >= day 2 and with same-day or 1-day discharge*
gen covid_hosp_date_day_cases_mab=covid_hospitalisation_outcome_da if covid_hosp_outcome_date2==covid_hosp_discharge_date2&covid_hosp_outcome_date2!=.&days_to_covid_admission>=2
replace covid_hosp_date_day_cases_mab=covid_hospitalisation_outcome_da if covid_hosp_outcome_date2==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure!=.&days_to_covid_admission>=2&(covid_hosp_discharge_date2-covid_hosp_outcome_date2)<=1&drug==1
replace covid_hospitalisation_outcome_da=. if covid_hosp_outcome_date2==covid_hosp_discharge_date2&covid_hosp_outcome_date2!=.&days_to_covid_admission>=2
replace covid_hospitalisation_outcome_da=. if covid_hosp_outcome_date2==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure!=.&days_to_covid_admission>=2&(covid_hosp_discharge_date2-covid_hosp_outcome_date2)<=1&drug==1
*check hosp_admission_method*
tab covid_hosp_admission_method,m
tab drug covid_hosp_admission_method, row chi
*by drug days_to_covid_admission, sort: count if covid_hospitalisation_outcome_da!=covid_hosp_date_emergency&covid_hospitalisation_outcome_da!=.


*count recorded day cases or regulars*
by drug, sort: count if covid_hosp_outcome_date0==.&covid_hosp_outcome_day_date0==start_date
by drug, sort: count if covid_hosp_outcome_date1==.&covid_hosp_outcome_day_date1==start_date+1
by drug, sort: count if covid_hosp_outcome_date2>covid_hosp_outcome_day_date2
by drug, sort: count if covid_hosp_outcome_date0==.&covid_hosp_outcome_day_date0==start_date&(covid_hosp_discharge_day_date0-covid_hosp_outcome_day_date0)==0
by drug, sort: count if covid_hosp_outcome_date1==.&covid_hosp_outcome_day_date1==start_date+1&(covid_hosp_discharge_day_date1-covid_hosp_outcome_day_date1)==0
by drug, sort: count if covid_hosp_outcome_date2>covid_hosp_outcome_day_date2&(covid_hosp_discharge_day_date2-covid_hosp_outcome_day_date2)==0
by drug, sort: count if covid_hosp_outcome_date0==.&covid_hosp_outcome_day_date0==start_date&(covid_hosp_discharge_day_date0-covid_hosp_outcome_day_date0)==1
by drug, sort: count if covid_hosp_outcome_date1==.&covid_hosp_outcome_day_date1==start_date+1&(covid_hosp_discharge_day_date1-covid_hosp_outcome_day_date1)==1
by drug, sort: count if covid_hosp_outcome_date2>covid_hosp_outcome_day_date2&(covid_hosp_discharge_day_date2-covid_hosp_outcome_day_date2)==1
by drug, sort: count if covid_hosp_outcome_date0==.&covid_hosp_outcome_day_date0==start_date&covid_hosp_outcome_day_date0==covid_hosp_date_mabs_day
by drug, sort: count if covid_hosp_outcome_date1==.&covid_hosp_outcome_day_date1==start_date+1&covid_hosp_outcome_day_date1==covid_hosp_date_mabs_day
by drug, sort: count if covid_hosp_outcome_date2>covid_hosp_outcome_day_date2&covid_hosp_outcome_day_date2==covid_hosp_date_mabs_day

count if covid_hosp_outcome_date0==.&covid_hosp_outcome_day_date0==start_date&sotrovimab_covid_therapeutics==.&casirivimab_covid_therapeutics==.
count if covid_hosp_outcome_date1==.&covid_hosp_outcome_day_date1==start_date+1&sotrovimab_covid_therapeutics==.&casirivimab_covid_therapeutics==.
count if covid_hosp_outcome_date2>covid_hosp_outcome_day_date2&sotrovimab_covid_therapeutics==.&casirivimab_covid_therapeutics==.
count if covid_hosp_outcome_date0==.&covid_hosp_outcome_day_date0==start_date&covid_hosp_outcome_day_date0==covid_hosp_date_mabs_day&sotrovimab_covid_therapeutics==.&casirivimab_covid_therapeutics==.
count if covid_hosp_outcome_date1==.&covid_hosp_outcome_day_date1==start_date+1&covid_hosp_outcome_day_date1==covid_hosp_date_mabs_day&sotrovimab_covid_therapeutics==.&casirivimab_covid_therapeutics==.
count if covid_hosp_outcome_date2>covid_hosp_outcome_day_date2&covid_hosp_outcome_day_date2==covid_hosp_date_mabs_day&sotrovimab_covid_therapeutics==.&casirivimab_covid_therapeutics==.







*define outcome and follow-up time*
gen study_end_date=mdy(11,28,2022)
gen start_date_29=start_date+28
by drug, sort: count if covid_hospitalisation_outcome_da!=.
by drug, sort: count if death_with_covid_on_the_death_ce!=.
by drug, sort: count if covid_hospitalisation_outcome_da==.&death_with_covid_on_the_death_ce!=.
by drug, sort: count if death_with_covid_on_the_death_ce==.&covid_hospitalisation_outcome_da!=.
*primary outcome*
gen event_date=min( covid_hospitalisation_outcome_da, death_with_covid_on_the_death_ce )
gen failure=(event_date!=.&event_date<=min(study_end_date,start_date_29,molnupiravir_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)) if drug==1
replace failure=(event_date!=.&event_date<=min(study_end_date,start_date_29,sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)) if drug==0
tab failure,m
gen end_date=event_date if failure==1
replace end_date=min(death_date, dereg_date, study_end_date, start_date_29,molnupiravir_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics,covid_hosp_date_day_cases_mab) if failure==0&drug==1
replace end_date=min(death_date, dereg_date, study_end_date, start_date_29,sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics,covid_hosp_date_day_cases_mab) if failure==0&drug==0
format %td event_date end_date study_end_date start_date_29

stset end_date ,  origin(start_date) failure(failure==1)
stcox drug

*secondary outcome: within 2 months*
gen start_date_2m=start_date+60
gen failure_2m=(event_date!=.&event_date<=min(study_end_date,start_date_2m,molnupiravir_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)) if drug==1
replace failure_2m=(event_date!=.&event_date<=min(study_end_date,start_date_2m,sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)) if drug==0
tab failure_2m,m
gen end_date_2m=event_date if failure_2m==1
replace end_date_2m=min(death_date, dereg_date, study_end_date, start_date_2m,molnupiravir_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics,covid_hosp_date_day_cases_mab) if failure_2m==0&drug==1
replace end_date_2m=min(death_date, dereg_date, study_end_date, start_date_2m,sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics,covid_hosp_date_day_cases_mab) if failure_2m==0&drug==0
format %td end_date_2m start_date_2m

stset end_date_2m ,  origin(start_date) failure(failure_2m==1)
stcox drug

*secondary outcome: within 3 months*
gen start_date_3m=start_date+90
gen failure_3m=(event_date!=.&event_date<=min(study_end_date,start_date_3m,molnupiravir_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)) if drug==1
replace failure_3m=(event_date!=.&event_date<=min(study_end_date,start_date_3m,sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)) if drug==0
tab failure_3m,m
gen end_date_3m=event_date if failure_3m==1
replace end_date_3m=min(death_date, dereg_date, study_end_date, start_date_3m,molnupiravir_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics,covid_hosp_date_day_cases_mab) if failure_3m==0&drug==1
replace end_date_3m=min(death_date, dereg_date, study_end_date, start_date_3m,sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics,covid_hosp_date_day_cases_mab) if failure_3m==0&drug==0
format %td end_date_3m start_date_3m

stset end_date_3m ,  origin(start_date) failure(failure_3m==1)
stcox drug

*secondary outcome: all-cause hosp/death within 29 days*
*correct all cause hosp date *
count if hospitalisation_outcome_date0!=start_date&hospitalisation_outcome_date0!=.
count if hospitalisation_outcome_date1!=(start_date+1)&hospitalisation_outcome_date1!=.
by drug, sort: count if hospitalisation_outcome_date0==hosp_discharge_date0&hospitalisation_outcome_date0!=.
by drug, sort: count if hospitalisation_outcome_date1==hosp_discharge_date1&hospitalisation_outcome_date1!=.
by drug, sort: count if hospitalisation_outcome_date0==covid_hosp_date_mabs_all_cause&covid_hosp_date_mabs_all_cause!=.
by drug, sort: count if hospitalisation_outcome_date1==covid_hosp_date_mabs_all_cause&covid_hosp_date_mabs_all_cause!=.
by drug, sort: count if hospitalisation_outcome_date0==covid_hosp_date_mabs_all_cause&hospitalisation_outcome_date0==hosp_discharge_date0&hospitalisation_outcome_date0!=.
by drug, sort: count if hospitalisation_outcome_date1==covid_hosp_date_mabs_all_cause&hospitalisation_outcome_date1==hosp_discharge_date1&hospitalisation_outcome_date1!=.
replace hospitalisation_outcome_date0=. if hospitalisation_outcome_date0==hosp_discharge_date0&hospitalisation_outcome_date0!=.
replace hospitalisation_outcome_date1=. if hospitalisation_outcome_date1==hosp_discharge_date1&hospitalisation_outcome_date1!=.
replace hospitalisation_outcome_date0=. if hospitalisation_outcome_date0==covid_hosp_date_mabs_all_cause&covid_hosp_date_mabs_all_cause!=.&drug==1
replace hospitalisation_outcome_date1=. if hospitalisation_outcome_date1==covid_hosp_date_mabs_all_cause&covid_hosp_date_mabs_all_cause!=.&drug==1

gen hospitalisation_outcome_date=hospitalisation_outcome_date2
replace hospitalisation_outcome_date=hospitalisation_outcome_date1 if hospitalisation_outcome_date1!=.
replace hospitalisation_outcome_date=hospitalisation_outcome_date0 if hospitalisation_outcome_date0!=.&drug==0
replace hospitalisation_outcome_date=hospitalisation_outcome_date0 if hospitalisation_outcome_date0!=.&drug==1

gen days_to_any_hosp_admission=hospitalisation_outcome_date-start_date if hospitalisation_outcome_date!=.
by drug days_to_any_hosp_admission, sort: count if hospitalisation_outcome_date!=.
by drug, sort: count if hospitalisation_outcome_date2==hosp_discharge_date2&hospitalisation_outcome_date2!=.&days_to_any_hosp_admission>=2
by drug days_to_any_hosp_admission, sort: count if hospitalisation_outcome_date2==hosp_discharge_date2&hospitalisation_outcome_date2!=.&days_to_any_hosp_admission>=2
by drug, sort: count if hospitalisation_outcome_date2==covid_hosp_date_mabs_all_cause&covid_hosp_date_mabs_all_cause!=.&days_to_any_hosp_admission>=2
by drug days_to_any_hosp_admission, sort: count if hospitalisation_outcome_date2==covid_hosp_date_mabs_all_cause&covid_hosp_date_mabs_all_cause!=.&days_to_any_hosp_admission>=2
by drug, sort: count if hospitalisation_outcome_date2==hosp_discharge_date2&hospitalisation_outcome_date2!=.&hospitalisation_outcome_date2==covid_hosp_date_mabs_all_cause&days_to_any_hosp_admission>=2
by drug, sort: count if hospitalisation_outcome_date2==covid_hosp_date_mabs_all_cause&hospitalisation_outcome_date2!=.&(hosp_discharge_date2 - hospitalisation_outcome_date2)==1&days_to_any_hosp_admission>=2
by drug, sort: count if hospitalisation_outcome_date2==covid_hosp_date_mabs_all_cause&hospitalisation_outcome_date2!=.&(hosp_discharge_date2 - hospitalisation_outcome_date2)==2&days_to_any_hosp_admission>=2
*ignore and censor day cases on or after day 2 from this analysis*
*ignore and censor admissions for mab procedure >= day 2 and with same-day or 1-day discharge*
gen hosp_date_day_cases_mab=hospitalisation_outcome_date if hospitalisation_outcome_date2==hosp_discharge_date2&hospitalisation_outcome_date2!=.&hospitalisation_outcome_date0==.&hospitalisation_outcome_date1==.
replace hosp_date_day_cases_mab=hospitalisation_outcome_date if hospitalisation_outcome_date2==covid_hosp_date_mabs_all_cause&covid_hosp_date_mabs_all_cause!=.&drug==1&(hosp_discharge_date2-hospitalisation_outcome_date2)<=1&hospitalisation_outcome_date0==.&hospitalisation_outcome_date1==.
replace hospitalisation_outcome_date=. if hospitalisation_outcome_date2==hosp_discharge_date2&hospitalisation_outcome_date2!=.&hospitalisation_outcome_date0==.&hospitalisation_outcome_date1==.
replace hospitalisation_outcome_date=. if hospitalisation_outcome_date2==covid_hosp_date_mabs_all_cause&covid_hosp_date_mabs_all_cause!=.&drug==1&(hosp_discharge_date2-hospitalisation_outcome_date2)<=1&hospitalisation_outcome_date0==.&hospitalisation_outcome_date1==.

gen event_date_allcause=min( death_date, hospitalisation_outcome_date,covid_hospitalisation_outcome_da )
gen failure_allcause=(event_date_allcause!=.&event_date_allcause<=min(study_end_date,start_date_29,molnupiravir_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)) if drug==1
replace failure_allcause=(event_date_allcause!=.&event_date_allcause<=min(study_end_date,start_date_29,sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)) if drug==0
tab failure_allcause,m
gen end_date_allcause=event_date_allcause if failure_allcause==1
replace end_date_allcause=min(death_date, dereg_date, study_end_date, start_date_29,molnupiravir_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics,hosp_date_day_cases_mab) if failure_allcause==0&drug==1
replace end_date_allcause=min(death_date, dereg_date, study_end_date, start_date_29,sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics,hosp_date_day_cases_mab) if failure_allcause==0&drug==0
format %td event_date_allcause end_date_allcause  

stset end_date_allcause ,  origin(start_date) failure(failure_allcause==1)
stcox drug

*sensitivity analysis for primary outcome: only emergency admissions, ignore non-emergency admissions*
*correct hosp date*
count if covid_hosp_date_emergency0!=start_date&covid_hosp_date_emergency0!=.
count if covid_hosp_date_emergency1!=(start_date+1)&covid_hosp_date_emergency1!=.
by drug, sort: count if covid_hosp_date_emergency0==covid_emerg_discharge_date0&covid_hosp_date_emergency0!=.
by drug, sort: count if covid_hosp_date_emergency1==covid_emerg_discharge_date1&covid_hosp_date_emergency1!=.
by drug, sort: count if covid_hosp_date_emergency0==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure!=.
by drug, sort: count if covid_hosp_date_emergency1==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure!=.
by drug, sort: count if covid_hosp_date_emergency0==covid_hosp_date_mabs_procedure&covid_hosp_date_emergency0==covid_emerg_discharge_date0&covid_hosp_date_emergency0!=.
by drug, sort: count if covid_hosp_date_emergency1==covid_hosp_date_mabs_procedure&covid_hosp_date_emergency1==covid_emerg_discharge_date1&covid_hosp_date_emergency1!=.
replace covid_hosp_date_emergency0=. if covid_hosp_date_emergency0==covid_emerg_discharge_date0&covid_hosp_date_emergency0!=.
replace covid_hosp_date_emergency1=. if covid_hosp_date_emergency1==covid_emerg_discharge_date1&covid_hosp_date_emergency1!=.
replace covid_hosp_date_emergency0=. if covid_hosp_date_emergency0==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure!=.&drug==1
replace covid_hosp_date_emergency1=. if covid_hosp_date_emergency1==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure!=.&drug==1

gen covid_hosp_date_emergency=covid_hosp_date_emergency2
replace covid_hosp_date_emergency=covid_hosp_date_emergency1 if covid_hosp_date_emergency1!=.
replace covid_hosp_date_emergency=covid_hosp_date_emergency0 if covid_hosp_date_emergency0!=.
*ignore and censor day cases on or after day 2 from this analysis*
*ignore and censor admissions for mab procedure >= day 2 and with same-day or 1-day discharge*
gen hosp_emergency_day_cases_mab=covid_hosp_date_emergency if covid_hosp_date_emergency2==covid_emerg_discharge_date2&covid_hosp_date_emergency2!=.&covid_hosp_date_emergency0==.&covid_hosp_date_emergency1==.
replace hosp_emergency_day_cases_mab=covid_hosp_date_emergency if covid_hosp_date_emergency2==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure!=.&drug==1&(covid_emerg_discharge_date2-covid_hosp_date_emergency2)<=1&covid_hosp_date_emergency0==.&covid_hosp_date_emergency1==.
replace covid_hosp_date_emergency=. if covid_hosp_date_emergency2==covid_emerg_discharge_date2&covid_hosp_date_emergency2!=.&covid_hosp_date_emergency0==.&covid_hosp_date_emergency1==.
replace covid_hosp_date_emergency=. if covid_hosp_date_emergency2==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure!=.&drug==1&(covid_emerg_discharge_date2-covid_hosp_date_emergency2)<=1&covid_hosp_date_emergency0==.&covid_hosp_date_emergency1==.

gen event_date_emergency=min( covid_hosp_date_emergency, death_with_covid_on_the_death_ce )
gen failure_emergency=(event_date_emergency!=.&event_date_emergency<=min(study_end_date,start_date_29,molnupiravir_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)) if drug==1
replace failure_emergency=(event_date_emergency!=.&event_date_emergency<=min(study_end_date,start_date_29,sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)) if drug==0
tab failure_emergency,m
gen end_date_emergency=event_date_emergency if failure_emergency==1
replace end_date_emergency=min(death_date, dereg_date, study_end_date, start_date_29,molnupiravir_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics,hosp_emergency_day_cases_mab) if failure_emergency==0&drug==1
replace end_date_emergency=min(death_date, dereg_date, study_end_date, start_date_29,sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics,hosp_emergency_day_cases_mab) if failure_emergency==0&drug==0
format %td event_date_emergency end_date_emergency  

stset end_date_emergency ,  origin(start_date) failure(failure_emergency==1)
stcox drug

*sensitivity analysis for primary outcome: not require covid as primary diagnosis*
*correct hosp date*
count if covid_hosp_date0_not_primary!=start_date&covid_hosp_date0_not_primary!=.
count if covid_hosp_date1_not_primary!=(start_date+1)&covid_hosp_date1_not_primary!=.
by drug, sort: count if covid_hosp_date0_not_primary==covid_discharge_date0_not_pri&covid_hosp_date0_not_primary!=.
by drug, sort: count if covid_hosp_date1_not_primary==covid_discharge_date1_not_pri&covid_hosp_date1_not_primary!=.
by drug, sort: count if covid_hosp_date0_not_primary==covid_hosp_date_mabs_not_pri&covid_hosp_date_mabs_not_pri!=.
by drug, sort: count if covid_hosp_date1_not_primary==covid_hosp_date_mabs_not_pri&covid_hosp_date_mabs_not_pri!=.
by drug, sort: count if covid_hosp_date0_not_primary==covid_hosp_date_mabs_not_pri&covid_hosp_date0_not_primary==covid_discharge_date0_not_pri&covid_hosp_date0_not_primary!=.
by drug, sort: count if covid_hosp_date1_not_primary==covid_hosp_date_mabs_not_pri&covid_hosp_date1_not_primary==covid_discharge_date1_not_pri&covid_hosp_date1_not_primary!=.
replace covid_hosp_date0_not_primary=. if covid_hosp_date0_not_primary==covid_discharge_date0_not_pri&covid_hosp_date0_not_primary!=.
replace covid_hosp_date1_not_primary=. if covid_hosp_date1_not_primary==covid_discharge_date1_not_pri&covid_hosp_date1_not_primary!=.
replace covid_hosp_date0_not_primary=. if covid_hosp_date0_not_primary==covid_hosp_date_mabs_not_pri&covid_hosp_date_mabs_not_pri!=.&drug==1
replace covid_hosp_date1_not_primary=. if covid_hosp_date1_not_primary==covid_hosp_date_mabs_not_pri&covid_hosp_date_mabs_not_pri!=.&drug==1

gen covid_hosp_date_not_primary=covid_hosp_date2_not_primary
replace covid_hosp_date_not_primary=covid_hosp_date1_not_primary if covid_hosp_date1_not_primary!=.
replace covid_hosp_date_not_primary=covid_hosp_date0_not_primary if covid_hosp_date0_not_primary!=.&drug==0
replace covid_hosp_date_not_primary=covid_hosp_date0_not_primary if covid_hosp_date0_not_primary!=.&drug==1

gen days_to_covid_adm_not_pri=covid_hosp_date_not_primary-start_date if covid_hosp_date2_not_primary!=.
by drug days_to_covid_adm_not_pri, sort: count if covid_hosp_date_not_primary!=.
by drug, sort: count if covid_hosp_date2_not_primary==covid_discharge_date2_not_pri&covid_hosp_date2_not_primary!=.&days_to_covid_adm_not_pri>=2
by drug days_to_covid_adm_not_pri, sort: count if covid_hosp_date2_not_primary==covid_discharge_date2_not_pri&covid_hosp_date2_not_primary!=.&days_to_covid_adm_not_pri>=2
by drug, sort: count if covid_hosp_date2_not_primary==covid_hosp_date_mabs_not_pri&covid_hosp_date_mabs_not_pri!=.&days_to_covid_adm_not_pri>=2
by drug days_to_covid_adm_not_pri, sort: count if covid_hosp_date2_not_primary==covid_hosp_date_mabs_not_pri&covid_hosp_date_mabs_not_pri!=.&days_to_covid_adm_not_pri>=2
by drug, sort: count if covid_hosp_date2_not_primary==covid_hosp_date_mabs_not_pri&covid_hosp_date2_not_primary!=.&covid_hosp_date2_not_primary==covid_discharge_date2_not_pri&days_to_covid_adm_not_pri>=2
by drug, sort: count if covid_hosp_date2_not_primary==covid_hosp_date_mabs_not_pri&covid_hosp_date2_not_primary!=.&(covid_discharge_date2_not_pri - covid_hosp_date2_not_primary)==1&days_to_covid_adm_not_pri>=2
by drug, sort: count if covid_hosp_date2_not_primary==covid_hosp_date_mabs_not_pri&covid_hosp_date2_not_primary!=.&(covid_discharge_date2_not_pri - covid_hosp_date2_not_primary)==2&days_to_covid_adm_not_pri>=2
*ignore and censor day cases on or after day 2 from this analysis*
*ignore and censor admissions for mab procedure >= day 2 and with same-day or 1-day discharge*
gen covid_day_cases_mab_not_pri=covid_hosp_date_not_primary if covid_hosp_date2_not_primary==covid_discharge_date2_not_pri&covid_hosp_date2_not_primary!=.&covid_hosp_date0_not_primary==.&covid_hosp_date1_not_primary==.
replace covid_day_cases_mab_not_pri=covid_hosp_date_not_primary if covid_hosp_date2_not_primary==covid_hosp_date_mabs_not_pri&covid_hosp_date_mabs_not_pri!=.&drug==1&covid_hosp_date0_not_primary==.&covid_hosp_date1_not_primary==.&(covid_discharge_date2_not_pri-covid_hosp_date2_not_primary)<=1
replace covid_hosp_date_not_primary=. if covid_hosp_date2_not_primary==covid_discharge_date2_not_pri&covid_hosp_date2_not_primary!=.&covid_hosp_date0_not_primary==.&covid_hosp_date1_not_primary==.
replace covid_hosp_date_not_primary=. if covid_hosp_date2_not_primary==covid_hosp_date_mabs_not_pri&covid_hosp_date_mabs_not_pri!=.&drug==1&covid_hosp_date0_not_primary==.&covid_hosp_date1_not_primary==.&(covid_discharge_date2_not_pri-covid_hosp_date2_not_primary)<=1

gen event_date_not_primary=min( covid_hosp_date_not_primary, death_with_covid_on_the_death_ce,covid_hospitalisation_outcome_da )
gen failure_not_primary=(event_date_not_primary!=.&event_date_not_primary<=min(study_end_date,start_date_29,molnupiravir_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)) if drug==1
replace failure_not_primary=(event_date_not_primary!=.&event_date_not_primary<=min(study_end_date,start_date_29,sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)) if drug==0
tab failure_not_primary,m
gen end_date_not_primary=event_date_not_primary if failure_not_primary==1
replace end_date_not_primary=min(death_date, dereg_date, study_end_date, start_date_29,molnupiravir_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics,covid_day_cases_mab_not_pri) if failure_not_primary==0&drug==1
replace end_date_not_primary=min(death_date, dereg_date, study_end_date, start_date_29,sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics,covid_day_cases_mab_not_pri) if failure_not_primary==0&drug==0
format %td event_date_not_primary end_date_not_primary  

stset end_date_not_primary ,  origin(start_date) failure(failure_not_primary==1)
stcox drug

*count censored due to second therapy*
count if failure==0&drug==1&min(molnupiravir_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)==end_date
count if failure==0&drug==0&min(molnupiravir_covid_therapeutics,sotrovimab_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)==end_date
*count covid death during day1-28 and before censor*
count if failure==1&drug==1&death_with_covid_on_the_death_ce<=min(study_end_date,start_date_29,molnupiravir_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)
count if failure==1&drug==0&death_with_covid_on_the_death_ce<=min(study_end_date,start_date_29,sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)
count if failure==1&drug==1&death_with_covid_underlying_date<=min(study_end_date,start_date_29,molnupiravir_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)
count if failure==1&drug==0&death_with_covid_underlying_date<=min(study_end_date,start_date_29,sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)
*count covid hosp during day1-28 and before censor*
by drug, sort: count if failure==1&covid_hospitalisation_outcome_da==end_date
*count covid death after covid hosp during day1-28 and before censor*
count if failure==1&covid_hospitalisation_outcome_da==end_date&drug==1&death_with_covid_on_the_death_ce>=covid_hospitalisation_outcome_da&death_with_covid_on_the_death_ce<=min(study_end_date,start_date_29,molnupiravir_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)
count if failure==1&covid_hospitalisation_outcome_da==end_date&drug==0&death_with_covid_on_the_death_ce>=covid_hospitalisation_outcome_da&death_with_covid_on_the_death_ce<=min(study_end_date,start_date_29,sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)
count if failure==1&covid_hospitalisation_outcome_da==end_date&drug==1&death_with_covid_underlying_date>=covid_hospitalisation_outcome_da&death_with_covid_underlying_date<=min(study_end_date,start_date_29,molnupiravir_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)
count if failure==1&covid_hospitalisation_outcome_da==end_date&drug==0&death_with_covid_underlying_date>=covid_hospitalisation_outcome_da&death_with_covid_underlying_date<=min(study_end_date,start_date_29,sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)
*count critical care within day1-28*
tab drug covid_hosp_critical_care,m row
tab drug covid_hosp_critical_care if failure==1&covid_hospitalisation_outcome_da==end_date,m row
*average hospital stay*
gen covid_hosp_days=covid_hosp_discharge_date2-covid_hospitalisation_outcome_da if covid_hosp_discharge_date2!=.&covid_hosp_discharge_date2>covid_hospitalisation_outcome_da
replace covid_hosp_days=covid_hosp_discharge_date1-covid_hospitalisation_outcome_da if covid_hosp_discharge_date1!=.&covid_hosp_discharge_date1>covid_hospitalisation_outcome_da
replace covid_hosp_days=covid_hosp_discharge_date0-covid_hospitalisation_outcome_da if covid_hosp_discharge_date0!=.&covid_hosp_discharge_date0>covid_hospitalisation_outcome_da
by drug, sort: sum covid_hosp_days if failure==1&covid_hospitalisation_outcome_da==end_date, de




*covariates* 
*10 high risk groups: downs_syndrome, solid_cancer, haematological_disease, renal_disease, liver_disease, imid, 
*immunosupression, hiv_aids, solid_organ_transplant, rare_neurological_conditions, high_risk_group_combined	
tab high_risk_cohort_covid_therapeut,m
by drug,sort: tab high_risk_cohort_covid_therapeut,m
gen downs_therapeutics= 1 if strpos(high_risk_cohort_covid_therapeut, "Downs syndrome")
gen solid_cancer_therapeutics=1 if strpos(high_risk_cohort_covid_therapeut, "solid cancer")
gen haema_disease_therapeutics=1 if strpos(high_risk_cohort_covid_therapeut, "haematological malignancies")
replace haema_disease_therapeutics=1 if strpos(high_risk_cohort_covid_therapeut, "sickle cell disease")
replace haema_disease_therapeutics=1 if strpos(high_risk_cohort_covid_therapeut, "haematological diseases")
replace haema_disease_therapeutics=1 if strpos(high_risk_cohort_covid_therapeut, "stem cell transplant")
gen renal_therapeutics= 1 if strpos(high_risk_cohort_covid_therapeut, "renal disease")
gen liver_therapeutics= 1 if strpos(high_risk_cohort_covid_therapeut, "liver disease")
gen imid_therapeutics= 1 if strpos(high_risk_cohort_covid_therapeut, "IMID")
gen immunosup_therapeutics= 1 if strpos(high_risk_cohort_covid_therapeut, "primary immune deficiencies")
gen hiv_aids_therapeutics= 1 if strpos(high_risk_cohort_covid_therapeut, "HIV or AIDS")
gen solid_organ_therapeutics= 1 if strpos(high_risk_cohort_covid_therapeut, "solid organ recipients")
replace solid_organ_therapeutics= 1 if strpos(high_risk_cohort_covid_therapeut, "solid organ transplant")
gen rare_neuro_therapeutics= 1 if strpos(high_risk_cohort_covid_therapeut, "rare neurological conditions")
*check if all diseases have been captured*
by drug,sort: count if high_risk_cohort_covid_therapeut!=""&high_risk_cohort_covid_therapeut!="other"& ///
   min(downs_therapeutics,solid_cancer_therapeutics,haema_disease_therapeutics,renal_therapeutics,liver_therapeutics,imid_therapeutics,immunosup_therapeutics,hiv_aids_therapeutics,solid_organ_therapeutics,rare_neuro_therapeutics)==.

replace oral_steroid_drugs_nhsd=. if oral_steroid_drug_nhsd_3m_count < 2 & oral_steroid_drug_nhsd_12m_count < 4
gen imid_nhsd=min(oral_steroid_drugs_nhsd, immunosuppresant_drugs_nhsd)
gen rare_neuro_nhsd = min(multiple_sclerosis_nhsd, motor_neurone_disease_nhsd, myasthenia_gravis_nhsd, huntingtons_disease_nhsd)

*gen downs_syndrome=(downs_syndrome_nhsd<=start_date|downs_therapeutics==1)
*gen solid_cancer=(cancer_opensafely_snomed<=start_date|solid_cancer_therapeutics==1)
*gen solid_cancer_new=(cancer_opensafely_snomed_new<=start_date|solid_cancer_therapeutics==1)
*gen haema_disease=( haematological_disease_nhsd <=start_date|haema_disease_therapeutics==1)
*gen renal_disease=( ckd_stage_5_nhsd <=start_date|renal_therapeutics==1)
*gen liver_disease=( liver_disease_nhsd <=start_date|liver_therapeutics==1)
*gen imid=( imid_nhsd <=start_date|imid_therapeutics==1)
*gen immunosupression=( immunosupression_nhsd <=start_date|immunosup_therapeutics==1)
*gen immunosupression_new=( immunosupression_nhsd_new <=start_date|immunosup_therapeutics==1)
*gen hiv_aids=( hiv_aids_nhsd <=start_date|hiv_aids_therapeutics==1)
*gen solid_organ=( solid_organ_transplant_nhsd<=start_date|solid_organ_therapeutics==1)
*gen solid_organ_new=( solid_organ_transplant_nhsd_new<=start_date|solid_organ_therapeutics==1)
*gen rare_neuro=( rare_neuro_nhsd <=start_date|rare_neuro_therapeutics==1)
*gen high_risk_group=(( downs_syndrome + solid_cancer + haema_disease + renal_disease + liver_disease + imid + immunosupression + hiv_aids + solid_organ + rare_neuro )>0)
*tab high_risk_group,m
*gen high_risk_group_new=(( downs_syndrome + solid_cancer_new + haema_disease + renal_disease + liver_disease + imid + immunosupression_new + hiv_aids + solid_organ_new + rare_neuro )>0)
*tab high_risk_group_new,m
*high risk group only based on codelists*
gen downs_syndrome=(downs_syndrome_nhsd<=start_date)
gen solid_cancer=(cancer_opensafely_snomed<=start_date)
gen solid_cancer_new=(cancer_opensafely_snomed_new<=start_date)
gen haema_disease=( haematological_disease_nhsd <=start_date)
gen renal_disease=( ckd_stage_5_nhsd <=start_date)
gen liver_disease=( liver_disease_nhsd <=start_date)
gen imid=( imid_nhsd <=start_date)
gen immunosupression=( immunosupression_nhsd <=start_date)
gen immunosupression_new=( immunosupression_nhsd_new <=start_date)
gen hiv_aids=( hiv_aids_nhsd <=start_date)
gen solid_organ=( solid_organ_transplant_nhsd<=start_date)
gen solid_organ_new=( solid_organ_transplant_nhsd_new<=start_date)
gen rare_neuro=( rare_neuro_nhsd <=start_date)
gen high_risk_group=(( downs_syndrome + solid_cancer + haema_disease + renal_disease + liver_disease + imid + immunosupression + hiv_aids + solid_organ + rare_neuro )>0)
tab high_risk_group,m
gen high_risk_group_new=(( downs_syndrome + solid_cancer_new + haema_disease + renal_disease + liver_disease + imid + immunosupression_new + hiv_aids + solid_organ_new + rare_neuro )>0)
tab high_risk_group_new,m

*Time between positive test and treatment*
gen d_postest_treat=start_date - covid_test_positive_date
tab d_postest_treat,m
replace d_postest_treat=. if d_postest_treat<0|d_postest_treat>7
gen d_postest_treat_g2=(d_postest_treat>=3) if d_postest_treat<=5
label define d_postest_treat_g2_Pax 0 "<3 days" 1 "3-5 days" 
label values d_postest_treat_g2 d_postest_treat_g2_Pax
gen d_postest_treat_missing=d_postest_treat_g2
replace d_postest_treat_missing=9 if d_postest_treat_g2==.
label define d_postest_treat_missing_Pax 0 "<3 days" 1 "3-5 days" 9 "missing" 
label values d_postest_treat_missing d_postest_treat_missing_Pax
*demo*
gen age_group3=(age>=40)+(age>=60)
label define age_group3_Paxlovid 0 "18-39" 1 "40-59" 2 ">=60" 
label values age_group3 age_group3_Paxlovid
tab age_group3,m
egen age_5y_band=cut(age), at(18,25,30,35,40,45,50,55,60,65,70,75,80,85,110) label
tab age_5y_band,m
gen age_50=(age>=50)
gen age_55=(age>=55)
gen age_60=(age>=60)

tab sex,m
rename sex sex_str
gen sex=0 if sex_str=="M"
replace sex=1 if sex_str=="F"
label define sex_Paxlovid 0 "Male" 1 "Female"
label values sex sex_Paxlovid

tab ethnicity,m
rename ethnicity ethnicity_with_missing_str
encode  ethnicity_with_missing_str ,gen(ethnicity_with_missing)
label list ethnicity_with_missing
gen ethnicity=ethnicity_with_missing
replace ethnicity=. if ethnicity_with_missing_str=="Missing"
label values ethnicity ethnicity_with_missing
gen White=1 if ethnicity==6
replace White=0 if ethnicity!=6&ethnicity!=.
gen White_with_missing=White
replace White_with_missing=9 if White==.


tab imd,m
replace imd=. if imd==0
label define imd_Paxlovid 1 "most deprived" 5 "least deprived"
label values imd imd_Paxlovid
gen imd_with_missing=imd
replace imd_with_missing=9 if imd==.

tab region_nhs,m
rename region_nhs region_nhs_str 
encode  region_nhs_str ,gen(region_nhs)
label list region_nhs

tab region_covid_therapeutics ,m
rename region_covid_therapeutics region_covid_therapeutics_str
encode  region_covid_therapeutics_str ,gen( region_covid_therapeutics )
label list region_covid_therapeutics

tab stp ,m
rename stp stp_str
encode  stp_str ,gen(stp)
label list stp
*combine stps with low N (<100) as "Other"*
by stp, sort: gen stp_N=_N if stp!=.
replace stp=99 if stp_N<100
tab stp ,m

tab rural_urban,m
replace rural_urban=. if rural_urban<1
replace rural_urban=3 if rural_urban==4
replace rural_urban=5 if rural_urban==6
replace rural_urban=7 if rural_urban==8
tab rural_urban,m
gen rural_urban_with_missing=rural_urban
replace rural_urban_with_missing=99 if rural_urban==.
*comor*
tab autism_nhsd,m
tab care_home_primis,m
tab dementia_nhsd,m
tab housebound_opensafely,m
tab learning_disability_primis,m
tab serious_mental_illness_nhsd,m
sum bmi,de
replace bmi=. if bmi<10|bmi>60
rename bmi bmi_all
*latest BMI within recent 10 years*
gen bmi=bmi_all if bmi_date_measured!=.&bmi_date_measured>=start_date-365*10&(age+((bmi_date_measured-start_date)/365)>=18)
gen bmi_5y=bmi_all if bmi_date_measured!=.&bmi_date_measured>=start_date-365*5&(age+((bmi_date_measured-start_date)/365)>=18)
gen bmi_2y=bmi_all if bmi_date_measured!=.&bmi_date_measured>=start_date-365*2&(age+((bmi_date_measured-start_date)/365)>=18)
gen bmi_group4=(bmi>=18.5)+(bmi>=25.0)+(bmi>=30.0) if bmi!=.
label define bmi_Paxlovid 0 "underweight" 1 "normal" 2 "overweight" 3 "obese"
label values bmi_group4 bmi_Paxlovid
gen bmi_g4_with_missing=bmi_group4
replace bmi_g4_with_missing=9 if bmi_group4==.
gen bmi_g3=bmi_group4
replace bmi_g3=1 if bmi_g3==0
label values bmi_g3 bmi_Paxlovid
gen bmi_g3_with_missing=bmi_g3
replace bmi_g3_with_missing=9 if bmi_g3==.
gen bmi_25=(bmi>=25) if bmi!=.
gen bmi_30=(bmi>=30) if bmi!=.

tab diabetes,m
tab chronic_cardiac_disease,m
tab hypertension,m
tab chronic_respiratory_disease,m
*vac and variant*
tab vaccination_status,m
rename vaccination_status vaccination_status_g6
gen vaccination_status=0 if vaccination_status_g6=="Un-vaccinated"|vaccination_status_g6=="Un-vaccinated (declined)"
replace vaccination_status=1 if vaccination_status_g6=="One vaccination"
replace vaccination_status=2 if vaccination_status_g6=="Two vaccinations"
replace vaccination_status=3 if vaccination_status_g6=="Three vaccinations"
replace vaccination_status=4 if vaccination_status_g6=="Four or more vaccinations"
label define vac_Paxlovid 0 "Un-vaccinated" 1 "One vaccination" 2 "Two vaccinations" 3 "Three vaccinations" 4 "Four or more vaccinations"
label values vaccination_status vac_Paxlovid
gen vaccination_3=1 if vaccination_status==3|vaccination_status==4
replace vaccination_3=0 if vaccination_status<3
gen vaccination_g3=vaccination_3 
replace vaccination_g3=2 if vaccination_status==4
gen pre_infection=(covid_test_positive_pre_date<=(covid_test_positive_date - 30)&covid_test_positive_pre_date>mdy(1,1,2020)&covid_test_positive_pre_date!=.)
tab pre_infection,m
tab sgtf,m
tab sgtf_new, m
label define sgtf_new_Paxlovid 0 "S gene detected" 1 "confirmed SGTF" 9 "NA"
label values sgtf_new sgtf_new_Paxlovid
tab variant_recorded ,m
*tab sgtf variant_recorded ,m
*Time between last vaccination and treatment*
gen d_vaccinate_treat=start_date - last_vaccination_date
sum d_vaccinate_treat,de
gen month_after_vaccinate=ceil(d_vaccinate_treat/30)
tab month_after_vaccinate,m
gen week_after_vaccinate=ceil(d_vaccinate_treat/7)
tab week_after_vaccinate,m
*calendar time*
gen month_after_campaign=ceil((start_date-mdy(12,15,2021))/30)
tab month_after_campaign,m
gen week_after_campaign=ceil((start_date-mdy(12,15,2021))/7)
tab week_after_campaign,m
gen day_after_campaign=start_date-mdy(12,15,2021)
sum day_after_campaign,de


*exclude those with contraindications for Pax*
*solid organ transplant*
tab drug if solid_organ==1|solid_organ_therapeutics==1|solid_organ_transplant_snomed<=start_date
tab drug if solid_organ_new==1|solid_organ_therapeutics==1|solid_organ_transplant_snomed<=start_date
*liver*
tab drug if advanced_decompensated_cirrhosis<=start_date
tab drug if decompensated_cirrhosis_icd10<=start_date
tab drug if ascitic_drainage_snomed<=start_date
tab drug if ascitic_drainage_snomed<=start_date&ascitic_drainage_snomed>=(start_date-3*365.25)
tab drug if liver_disease_nhsd_icd10<=start_date
*renal*
tab drug if renal_disease==1|renal_therapeutics==1
tab drug if ckd_stages_3_5<=start_date
tab drug ckd_primis_stage,row
replace ckd_primis_stage=. if ckd_primis_stage_date>start_date
tab drug if ckd3_icd10<=start_date|ckd4_icd10<=start_date|ckd5_icd10<=start_date
tab drug if kidney_transplant<=start_date|kidney_transplant_icd10<=start_date|kidney_transplant_procedure<=start_date
tab drug if dialysis<=start_date|dialysis_icd10<=start_date|dialysis_procedure<=start_date
*egfr: adapted from https://github.com/opensafely/COVID-19-vaccine-breakthrough/blob/updates-feb/analysis/data_process.R*
tab creatinine_operator_ctv3,m
replace creatinine_ctv3 = . if !inrange(creatinine_ctv3, 20, 3000)| creatinine_ctv3_date>start_date
tab creatinine_operator_ctv3 if creatinine_ctv3!=.,m
replace creatinine_ctv3 = creatinine_ctv3/88.4
gen min_creatinine_ctv3=.
replace min_creatinine_ctv3 = (creatinine_ctv3/0.7)^-0.329 if sex==1
replace min_creatinine_ctv3 = (creatinine_ctv3/0.9)^-0.411 if sex==0
replace min_creatinine_ctv3 = 1 if min_creatinine_ctv3<1
gen max_creatinine_ctv3=.
replace max_creatinine_ctv3 = (creatinine_ctv3/0.7)^-1.209 if sex==1
replace max_creatinine_ctv3 = (creatinine_ctv3/0.9)^-1.209 if sex==0
replace max_creatinine_ctv3 = 1 if max_creatinine_ctv3>1
gen egfr_creatinine_ctv3 = min_creatinine_ctv3*max_creatinine_ctv3*141*(0.993^age_creatinine_ctv3) if age_creatinine_ctv3>0&age_creatinine_ctv3<=120
replace egfr_creatinine_ctv3 = egfr_creatinine_ctv3*1.018 if sex==1

tab creatinine_operator_snomed,m
tab creatinine_operator_snomed if creatinine_snomed!=.,m
replace creatinine_snomed = . if !inrange(creatinine_snomed, 20, 3000)| creatinine_snomed_date>start_date
replace creatinine_snomed_date = creatinine_short_snomed_date if missing(creatinine_snomed)
replace creatinine_operator_snomed = creatinine_operator_short_snomed if missing(creatinine_snomed)
replace age_creatinine_snomed = age_creatinine_short_snomed if missing(creatinine_snomed)
replace creatinine_snomed = creatinine_short_snomed if missing(creatinine_snomed)
replace creatinine_snomed = . if !inrange(creatinine_snomed, 20, 3000)| creatinine_snomed_date>start_date
replace creatinine_snomed = creatinine_snomed/88.4
gen min_creatinine_snomed=.
replace min_creatinine_snomed = (creatinine_snomed/0.7)^-0.329 if sex==1
replace min_creatinine_snomed = (creatinine_snomed/0.9)^-0.411 if sex==0
replace min_creatinine_snomed = 1 if min_creatinine_snomed<1
gen max_creatinine_snomed=.
replace max_creatinine_snomed = (creatinine_snomed/0.7)^-1.209 if sex==1
replace max_creatinine_snomed = (creatinine_snomed/0.9)^-1.209 if sex==0
replace max_creatinine_snomed = 1 if max_creatinine_snomed>1
gen egfr_creatinine_snomed = min_creatinine_snomed*max_creatinine_snomed*141*(0.993^age_creatinine_snomed) if age_creatinine_snomed>0&age_creatinine_snomed<=120
replace egfr_creatinine_snomed = egfr_creatinine_snomed*1.018 if sex==1

tab eGFR_operator if eGFR_record!=.,m
tab eGFR_short_operator if eGFR_short_record!=.,m
tab drug if (egfr_creatinine_ctv3<60&creatinine_operator_ctv3!="<")|(egfr_creatinine_snomed<60&creatinine_operator_snomed!="<")|(eGFR_record<60&eGFR_record>0&eGFR_operator!=">"&eGFR_operator!=">=")|(eGFR_short_record<60&eGFR_short_record>0&eGFR_short_operator!=">"&eGFR_short_operator!=">=")

*drug interactions*
tab drug if drugs_do_not_use<=start_date
tab drug if drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-3*365.25)
tab drug if drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-365.25)
tab drug if drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-180)
tab drug if drugs_consider_risk<=start_date
tab drug if drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-3*365.25)
tab drug if drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-365.25)
tab drug if drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-180)
gen drugs_do_not_use_contra=(drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-180))
gen drugs_consider_risk_contra=(drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-180))

save ./output/sensitivity.dta, replace

drop if solid_organ_new==1|solid_organ_therapeutics==1|solid_organ_transplant_snomed<=start_date
drop if advanced_decompensated_cirrhosis<=start_date|decompensated_cirrhosis_icd10<=start_date|ascitic_drainage_snomed<=start_date|liver_disease_nhsd_icd10<=start_date
drop if renal_disease==1|renal_therapeutics==1|ckd_stages_3_5<=start_date|ckd_primis_stage==3|ckd_primis_stage==4|ckd_primis_stage==5|ckd3_icd10<=start_date|ckd4_icd10<=start_date|ckd5_icd10<=start_date
drop if kidney_transplant<=start_date|kidney_transplant_icd10<=start_date|kidney_transplant_procedure<=start_date
drop if dialysis<=start_date|dialysis_icd10<=start_date|dialysis_procedure<=start_date
drop if (egfr_creatinine_ctv3<60&creatinine_operator_ctv3!="<")|(egfr_creatinine_snomed<60&creatinine_operator_snomed!="<")|(eGFR_record<60&eGFR_record>0&eGFR_operator!=">"&eGFR_operator!=">=")|(eGFR_short_record<60&eGFR_short_record>0&eGFR_short_operator!=">"&eGFR_short_operator!=">=")
*drop if drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-365.25)
*drop if drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-365.25)
drop if drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-180)
*drop if drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-180)

tab drug if liver_disease_nhsd_snomed<=start_date
tab drug if liver_disease==1
tab drug if drugs_do_not_use<=start_date
tab drug if drugs_consider_risk<=start_date

*clean covariates*
tab month_after_vaccinate,m
*combine month7 and over due to small N*
replace month_after_vaccinate=7 if month_after_vaccinate>=7&month_after_vaccinate!=.
gen month_after_vaccinate_missing=month_after_vaccinate
replace month_after_vaccinate_missing=99 if month_after_vaccinate_missing==.
*calendar time*
tab week_after_campaign,m
*combine 9/10 and 26/27 due to small N*
*replace week_after_campaign=10 if week_after_campaign==9
*replace week_after_campaign=26 if week_after_campaign==27
*combine stps with low N (<100) as "Other"*
drop stp_N
by stp, sort: gen stp_N=_N if stp!=.
replace stp=99 if stp_N<100
tab stp ,m


*descriptives by drug groups*
by drug,sort: sum age,de
ttest age , by( drug )
by drug,sort: sum bmi,de
ttest bmi, by( drug )
sum d_postest_treat ,de
by drug,sort: sum d_postest_treat ,de
ttest d_postest_treat , by( drug )
ranksum d_postest_treat,by(drug)
sum week_after_campaign,de
by drug,sort: sum week_after_campaign,de
ttest week_after_campaign , by( drug )
ranksum week_after_campaign,by(drug)
sum week_after_vaccinate,de
by drug,sort: sum week_after_vaccinate,de
ttest week_after_vaccinate , by( drug )
ranksum week_after_vaccinate,by(drug)
sum d_vaccinate_treat,de
by drug,sort: sum d_vaccinate_treat,de
ttest d_vaccinate_treat , by( drug )
ranksum d_vaccinate_treat,by(drug)

tab drug sex,row chi
tab drug ethnicity,row chi
tab drug White,row chi
tab drug imd,row chi
ranksum imd,by(drug)
tab drug rural_urban,row chi
ranksum rural_urban,by(drug)
tab drug region_nhs,row chi
tab drug region_covid_therapeutics,row chi
*need to address the error of "too many values"*
tab stp if drug==0
tab stp if drug==1
tab drug age_group3 ,row chi
tab drug d_postest_treat_g2 ,row chi
tab drug d_postest_treat ,row
tab drug downs_syndrome ,row chi
tab drug solid_cancer ,row chi
tab drug solid_cancer_new ,row chi
tab drug haema_disease ,row chi
tab drug renal_disease ,row chi
tab drug liver_disease ,row chi
tab drug imid ,row chi
tab drug immunosupression ,row chi
tab drug immunosupression_new ,row chi
tab drug hiv_aids ,row chi
tab drug solid_organ ,row chi
tab drug solid_organ_new ,row chi
tab drug rare_neuro ,row chi
tab drug high_risk_group ,row chi
tab drug high_risk_group_new ,row chi
tab drug autism_nhsd ,row chi
tab drug care_home_primis ,row chi
tab drug dementia_nhsd ,row chi
tab drug housebound_opensafely ,row chi
tab drug learning_disability_primis ,row chi
tab drug serious_mental_illness_nhsd ,row chi
tab drug bmi_group4 ,row chi
tab drug bmi_g3 ,row chi
tab drug diabetes ,row chi
tab drug chronic_cardiac_disease ,row chi
tab drug hypertension ,row chi
tab drug chronic_respiratory_disease ,row chi
tab drug vaccination_status ,row chi
tab drug month_after_vaccinate,row chi
tab drug sgtf ,row chi
tab drug sgtf_new ,row chi
tab drug drugs_consider_risk_contra,row chi
*tab drug variant_recorded ,row chi
tab drug if covid_test_positive_pre_date!=.
stset end_date ,  origin(start_date) failure(failure==1)
stcox drug


*check treatment status*
count if drug==0&paxlovid_covid_therapeutics==paxlovid_covid_approved
count if drug==0&paxlovid_covid_therapeutics==paxlovid_covid_complete
count if drug==0&paxlovid_covid_therapeutics==paxlovid_covid_not_start
count if drug==0&paxlovid_covid_therapeutics==paxlovid_covid_stopped
count if drug==0&paxlovid_covid_approved!=.
count if drug==0&paxlovid_covid_complete!=.
count if drug==0&paxlovid_covid_not_start!=.
count if drug==0&paxlovid_covid_stopped!=.
count if drug==1&sotrovimab_covid_therapeutics==sotrovimab_covid_approved
count if drug==1&sotrovimab_covid_therapeutics==sotrovimab_covid_complete
count if drug==1&sotrovimab_covid_therapeutics==sotrovimab_covid_not_start
count if drug==1&sotrovimab_covid_therapeutics==sotrovimab_covid_stopped
count if drug==1&sotrovimab_covid_approved!=.
count if drug==1&sotrovimab_covid_complete!=.
count if drug==1&sotrovimab_covid_not_start!=.
count if drug==1&sotrovimab_covid_stopped!=.


*compare characteristics between those with detected high-risk group category and those without*
by drug,sort: tab high_risk_group_new,m

by drug,sort: sum age if high_risk_group_new==0,de
by drug,sort: sum bmi if high_risk_group_new==0,de
by drug,sort: sum d_postest_treat if high_risk_group_new==0,de
by drug,sort: sum week_after_campaign if high_risk_group_new==0,de
by drug,sort: sum week_after_vaccinate if high_risk_group_new==0,de
by drug,sort: sum d_vaccinate_treat if high_risk_group_new==0,de

tab drug sex if high_risk_group_new==0,row chi
tab drug ethnicity if high_risk_group_new==0,row chi
tab drug White if high_risk_group_new==0,row chi
tab drug imd if high_risk_group_new==0,row chi
tab drug rural_urban if high_risk_group_new==0,row chi
tab drug region_nhs if high_risk_group_new==0,row chi
tab drug region_covid_therapeutics if high_risk_group_new==0,row chi
tab drug age_group3  if high_risk_group_new==0,row chi
tab drug d_postest_treat_g2  if high_risk_group_new==0,row chi
tab drug d_postest_treat  if high_risk_group_new==0,row
tab drug downs_therapeutics  if high_risk_group_new==0,row
tab drug solid_cancer_therapeutics  if high_risk_group_new==0,row
tab drug haema_disease_therapeutics  if high_risk_group_new==0,row
tab drug renal_therapeutics  if high_risk_group_new==0,row
tab drug liver_therapeutics  if high_risk_group_new==0,row
tab drug imid_therapeutics  if high_risk_group_new==0,row
tab drug immunosup_therapeutics  if high_risk_group_new==0,row
tab drug hiv_aids_therapeutics  if high_risk_group_new==0,row
tab drug solid_organ_therapeutics  if high_risk_group_new==0,row
tab drug rare_neuro_therapeutics  if high_risk_group_new==0,row
tab drug autism_nhsd  if high_risk_group_new==0,row chi
tab drug care_home_primis  if high_risk_group_new==0,row chi
tab drug dementia_nhsd  if high_risk_group_new==0,row chi
tab drug housebound_opensafely  if high_risk_group_new==0,row chi
tab drug learning_disability_primis  if high_risk_group_new==0,row chi
tab drug serious_mental_illness_nhsd  if high_risk_group_new==0,row chi
tab drug bmi_group4  if high_risk_group_new==0,row chi
tab drug bmi_g3  if high_risk_group_new==0,row chi
tab drug diabetes  if high_risk_group_new==0,row chi
tab drug chronic_cardiac_disease  if high_risk_group_new==0,row chi
tab drug hypertension  if high_risk_group_new==0,row chi
tab drug chronic_respiratory_disease  if high_risk_group_new==0,row chi
tab drug vaccination_status  if high_risk_group_new==0,row chi
tab drug month_after_vaccinate if high_risk_group_new==0,row chi
tab drug drugs_consider_risk_contra if high_risk_group_new==0,row chi
tab failure drug if high_risk_group_new==0,m col



drop if high_risk_group_new==0
*descriptives by drug groups*
by drug,sort: sum age,de
ttest age , by( drug )
by drug,sort: sum bmi,de
ttest bmi, by( drug )
sum d_postest_treat ,de
by drug,sort: sum d_postest_treat ,de
ttest d_postest_treat , by( drug )
ranksum d_postest_treat,by(drug)
sum week_after_campaign,de
by drug,sort: sum week_after_campaign,de
ttest week_after_campaign , by( drug )
ranksum week_after_campaign,by(drug)
sum week_after_vaccinate,de
by drug,sort: sum week_after_vaccinate,de
ttest week_after_vaccinate , by( drug )
ranksum week_after_vaccinate,by(drug)
sum d_vaccinate_treat,de
by drug,sort: sum d_vaccinate_treat,de
ttest d_vaccinate_treat , by( drug )
ranksum d_vaccinate_treat,by(drug)

tab drug sex,row chi
tab drug ethnicity,row chi
tab drug White,row chi
tab drug imd,row chi
ranksum imd,by(drug)
tab drug rural_urban,row chi
ranksum rural_urban,by(drug)
tab drug region_nhs,row chi
tab drug region_covid_therapeutics,row chi
*need to address the error of "too many values"*
tab stp if drug==0
tab stp if drug==1
tab drug age_group3 ,row chi
tab drug d_postest_treat_g2 ,row chi
tab drug d_postest_treat ,row
tab drug downs_syndrome ,row chi
tab drug solid_cancer ,row chi
tab drug solid_cancer_new ,row chi
tab drug haema_disease ,row chi
tab drug renal_disease ,row chi
tab drug liver_disease ,row chi
tab drug imid ,row chi
tab drug immunosupression ,row chi
tab drug immunosupression_new ,row chi
tab drug hiv_aids ,row chi
tab drug solid_organ ,row chi
tab drug solid_organ_new ,row chi
tab drug rare_neuro ,row chi
tab drug high_risk_group ,row chi
tab drug autism_nhsd ,row chi
tab drug care_home_primis ,row chi
tab drug dementia_nhsd ,row chi
tab drug housebound_opensafely ,row chi
tab drug learning_disability_primis ,row chi
tab drug serious_mental_illness_nhsd ,row chi
tab drug bmi_group4 ,row chi
tab drug bmi_g3 ,row chi
tab drug diabetes ,row chi
tab drug chronic_cardiac_disease ,row chi
tab drug hypertension ,row chi
tab drug chronic_respiratory_disease ,row chi
tab drug vaccination_status ,row chi
tab drug month_after_vaccinate,row chi
tab drug month_after_campaign,row chi
tab drug sgtf ,row chi
tab drug sgtf_new ,row chi
tab drug pre_infection,row chi
tab drug drugs_consider_risk_contra,row chi
*tab drug variant_recorded ,row chi
tab drug if covid_test_positive_pre_date!=.
stset end_date ,  origin(start_date) failure(failure==1)
stcox drug

*recode Paxlovid as 1*
replace drug=1-drug
label define drug_Paxlovid2 0 "sotrovimab" 1 "Paxlovid"
label values drug drug_Paxlovid2
*gen splines*
mkspline age_spline = age, cubic nknots(4)
mkspline calendar_day_spline = day_after_campaign, cubic nknots(4)


save ./output/main.dta, replace

log close




