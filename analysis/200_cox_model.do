********************************************************************************
*
*	Do-file:			cox model
*	Project:			Sotrovimab-Paxlovid-Molnupiravir
*   Date:  				06/4/23
*	Programmed by:		Katie Bechman
* 	Description:		run cox models
*	Data used:			main.dta
*	Data created:		coxoutput
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


/* SET Index date ===========================================================*/
global indexdate 			= "01/03/2020"

* Event
foreach var of varlist ae_spc* ae_imae* ae_all* covid_hosp_date all_hosp_date{
				display "`var'"
				by drug, sort: count if `var'!=. 
}
foreach var of varlist ae_spc_all ae_imae_all ae_all covid_hosp_date all_hosp_date{
				gen fail_`var' = `var'
} 
 	
 
	*generate censor date
	gen diecensor = mdy(10,01,2020)
	format diecensor %td
	
	egen stopdied = rmin(died_ons_date diecensor)
	egen stopicu = rmin(icu_or_death_covid_date diecensor)
	egen stopswab = rmin(first_pos_test_sgss_date diecensor)


	
	gen exitdied = died_ons_covid_flag_any
	gen exiticu = icu_or_death_covid
	gen exitswab = first_pos_test_sgss 



 

gen event_date=ae_all
by drug, sort: count if event_date!=.
by drug, sort:  count if event_date<start_date & event_date!=.
drop if event_date<start_date & event_date!=. & drug==0 // REMOVE control from control group if event occurs before 'treatment start'
* Failure - if event occurs before end date or second drug initiation 
foreach var of varlist
gen failure = (event_date!=. & event_date<= min(study_end_date, start_date_29, paxlovid_date_started, molnupiravir_date_started, remdesivir_date_started, casirivimab_date_started)) if drug==1
replace failure = (event_date!=. & event_date<= min(study_end_date, start_date_29, sotrovimab_date_started, molnupiravir_date_started, remdesivir_date_started, casirivimab_date_started)) if drug==2
replace failure = (event_date!=. & event_date<= min(study_end_date, start_date_29, sotrovimab_date_started, paxlovid_date_started, remdesivir_date_started, casirivimab_date_started)) if drug==3
replace failure = (event_date!=. & event_date<= min(study_end_date, start_date_29, sotrovimab_date_started, paxlovid_date_started, molnupiravir_date_started, remdesivir_date_started, casirivimab_date_started)) if drug==0
tab drug failure, m
// Censoring
gen end_date=event_date if failure==1
replace end_date=min(death_date, dereg_date, study_end_date, start_date_29, paxlovid_date_started, molnupiravir_date_started, remdesivir_date_started, casirivimab_date_started) if failure==0&drug==1
replace end_date=min(death_date, dereg_date, study_end_date, start_date_29, sotrovimab_date_started, molnupiravir_date_started, remdesivir_date_started, casirivimab_date_started) if failure==0&drug==2
replace end_date=min(death_date, dereg_date, study_end_date, start_date_29, sotrovimab_date_started, paxlovid_date_started, remdesivir_date_started, casirivimab_date_started) if failure==0&drug==3
replace end_date=min(death_date, dereg_date, study_end_date, start_date_29, sotrovimab_date_started, paxlovid_date_started, molnupiravir_date_started, remdesivir_date_started, casirivimab_date_started) if failure==0&drug==0
format %td event_date end_date study_end_date start_date_29


/* SET Index date ===========================================================*/
global crude
global agesex i.agegroup male
global adjusted_main i.agegroup male i.ethnicity i.imd bmicat bowel skin joint chronic_cardiac_disease cancer stroke i.diabcat steroidcat 
global adjusted_extra i.agegroup male i.ethnicity i.imd bmicat bowel skin joint chronic_cardiac_disease cancer stroke i.diabcat steroidcat i.ckd chronic_liver_disease chronic_respiratory_disease
tempname coxoutput
	postfile `coxoutput' str20(cohort) str20(model) str20(failure) ///
		ptime_exposed events_exposed rate_exposed /// 
		ptime_comparator events_comparator rate_comparator hr lc uc ///
		using $projectdir/output/data/cox_model_summary, replace						

use $projectdir/output/data/file_`f', replace


	 

	foreach fail in died icu swab {

		stset stop`fail', id(patient_id) failure(fail`fail'==1) origin(time enter_date)  enter(time enter_date) scale(365.25) 
						
		foreach model in crude agesex adjusted_main adjusted_extra {
				
			stcox `f' $`model', vce(robust)
						matrix b = r(table)
						local hr = b[1,1]
						local lc = b[5,1]
						local uc = b[6,1]

			stptime if `f' == 1
						local rate_exposed = `r(rate)'
						local ptime_exposed = `r(ptime)'
						local events_exposed .
						if `r(failures)' == 0 | `r(failures)' > 5 local events_exposed `r(failures)'
						
			stptime if `f' == 0
						local rate_comparator = `r(rate)'
						local ptime_comparator = `r(ptime)'
						local events_comparator .
						if `r(failures)' == 0 | `r(failures)' > 5 local events_comparator `r(failures)'

			post `coxoutput' ("`f'") ("`model'") ("`fail'") (`ptime_exposed') (`events_exposed') (`rate_exposed') ///
						(`ptime_comparator') (`events_comparator') (`rate_comparator') ///
						(`hr') (`lc') (`uc')
		}
	}
}

postclose `coxoutput'




log close



* Follow-up time
gen study_end_date=mdy(01,08,2023)
gen start_date_29=start_date+28
stset end_date, origin(start_date) failure(failure==1)
stcox i.drug

gen event_date_allcause=min( death_date, hospitalisation_outcome_date,covid_hospitalisation_outcome_da )
gen failure_allcause=(event_date_allcause!=.&event_date_allcause<=min(study_end_date,start_date_29,sotrovimab_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)) if drug==1
replace failure_allcause=(event_date_allcause!=.&event_date_allcause<=min(study_end_date,start_date_29,sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics)) if drug==0
tab failure_allcause,m
gen end_date_allcause=event_date_allcause if failure_allcause==1
replace end_date_allcause=min(death_date, dereg_date, study_end_date, start_date_29,sotrovimab_covid_therapeutics,paxlovid_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics,hosp_date_day_cases) if failure_allcause==0&drug==1
replace end_date_allcause=min(death_date, dereg_date, study_end_date, start_date_29,sotrovimab_covid_therapeutics,molnupiravir_covid_therapeutics,remdesivir_covid_therapeutics,casirivimab_covid_therapeutics,hosp_date_day_cases) if failure_allcause==0&drug==0
format %td event_date_allcause end_date_allcause  

stset end_date_allcause ,  origin(start_date) failure(failure_allcause==1)
stcox drug



















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