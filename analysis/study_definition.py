from cohortextractor import StudyDefinition, patients, codelist_from_csv, codelist, filter_codes_by_category, combine_codelists, Measure
import datetime
from codelists import *

## Define study time variables
from datetime import timedelta, date, datetime 
campaign_start = "2021-12-16"
end_date = date.today().isoformat()

def comorbidity_snomed(dx_codelist):
  return patients.with_these_clinical_events(
      dx_codelist,
      returning="date",
      on_or_before ="covid_test_positive_date",
      find_last_match_in_period=True,
      date_format="YYYY-MM-DD",
      return_expectations={
          "incidence": 0.2,
          "date": {"earliest": "1950-01-01", "latest": end_date},
      },
  )
def comorbidity_snomed_6m(dx_codelist):
  return patients.with_these_clinical_events(
      dx_codelist,
      returning="date",
      between = ["covid_test_positive_date - 6 months", "covid_test_positive_date"],
      find_last_match_in_period=True,
      date_format="YYYY-MM-DD",
      return_expectations={
          "incidence": 0.2,
          "date": {"earliest": "2020-07-16", "latest": end_date},
      },
  )
def comorbidity_snomed_12m(dx_codelist):
  return patients.with_these_clinical_events(
      dx_codelist,
      returning="date",
      between = ["covid_test_positive_date - 12 months", "covid_test_positive_date"],
      find_last_match_in_period=True,
      date_format="YYYY-MM-DD",
      return_expectations={
          "incidence": 0.2,
          "date": {"earliest": "2020-12-16", "latest": end_date},
      },
  )
def comorbidity_icd(dx_codelist):
  return patients.admitted_to_hospital(
      with_these_diagnoses = dx_codelist,
      returning="date_admitted",
      on_or_before = "covid_test_positive_date",
      find_last_match_in_period=True,
      date_format="YYYY-MM-DD",
      return_expectations={
          "incidence": 0.2,
          "date": {"earliest": "1950-01-01", "latest": end_date},
      },
  )
def comorbidity_icd_12m(dx_codelist):
  return patients.admitted_to_hospital(
      with_these_diagnoses = dx_codelist,
      returning="date_admitted",
      between = ["covid_test_positive_date - 12 months", "covid_test_positive_date"],
      find_last_match_in_period=True,
      date_format="YYYY-MM-DD",
      return_expectations={
          "incidence": 0.2,
          "date": {"earliest": "2020-12-16", "latest": end_date},
      },
  )
def comorbidity_ops(dx_codelist):
  return patients.admitted_to_hospital(
      with_these_procedures=dx_codelist,
      returning="date_admitted",
      on_or_before = "covid_test_positive_date",
      find_last_match_in_period=True,
      date_format="YYYY-MM-DD",
      return_expectations={
          "incidence": 0.01,
          "date": {"earliest": "2020-12-16", "latest": end_date},
          "rate": "exponential_increase",
      },
  )
def comorbidity_ops_12m(dx_codelist):
  return patients.admitted_to_hospital(
      with_these_procedures=dx_codelist,
      returning="date_admitted",
      between = ["covid_test_positive_date - 12 months", "covid_test_positive_date"],
      find_last_match_in_period=True,
      date_format="YYYY-MM-DD",
      return_expectations={
          "incidence": 0.01,
          "date": {"earliest": "2020-02-01", "latest": end_date},
          "rate": "exponential_increase",
      },
  )
def drug_6m(dx_codelist):
  return patients.with_these_medications(
      dx_codelist,
      returning="date",
      between = ["covid_test_positive_date - 6 months", "covid_test_positive_date"],
      find_last_match_in_period=True,
      date_format="YYYY-MM-DD",
      return_expectations={
          "incidence": 0.01,
          "date": {"earliest": "2021-06-16", "latest": end_date},
      },
  )
def drug_12m(dx_codelist):
  return patients.with_these_medications(
      dx_codelist,
      returning="date",
      between = ["covid_test_positive_date - 12 months", "covid_test_positive_date"],
      find_last_match_in_period=True,
      date_format="YYYY-MM-DD",
      return_expectations={
          "incidence": 0.01,
          "date": {"earliest": "2020-12-16", "latest": end_date},
      },
  )
def high_cost_drugs_3m(dx_codelist):
  return patients.with_high_cost_drugs(
      returning="date",
      between = ["covid_test_positive_date - 3 months", "covid_test_positive_date"],
      find_last_match_in_period=True,
      date_format="YYYY-MM",
      return_expectations={
          "incidence": 0.01,
          "date": {"earliest": "2021-09-16", "latest": end_date},
      },
  ) 
def covid_therapeutics(dx_codelist):
  return patients.with_covid_therapeutics(
      dx_codelist,
      with_these_indications = "non_hospitalised",
      between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
      find_last_match_in_period=True,
      returning="date",
      date_format="YYYY-MM-DD",
      return_expectations={
          "incidence": 0.2,
          "date": {"earliest": "2021-12-16"},
      },
  )
def covid_therapeutics_approved(dx_codelist):
  return patients.with_covid_therapeutics(
      dx_codelist,
      with_these_statuses = ["Approved"],
      with_these_indications = "non_hospitalised",
      between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
      find_last_match_in_period=True,
      returning="date",
      date_format="YYYY-MM-DD",
      return_expectations={
          "incidence": 0.2,
          "date": {"earliest": "2021-12-16"},
      },
  )
def covid_therapeutics_complete(dx_codelist):
  return patients.with_covid_therapeutics(
      dx_codelist,
      with_these_statuses = ["Treatment Complete"], 
      with_these_indications = "non_hospitalised",
      between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
      find_last_match_in_period=True,
      returning="date",
      date_format="YYYY-MM-DD",
      return_expectations={
          "incidence": 0.2,
          "date": {"earliest": "2021-12-16"},
      },
  )
def covid_therapeutics_could_not_start(dx_codelist):
  return patients.with_covid_therapeutics(
      dx_codelist,
      with_these_indications = "non_hospitalised",
      with_these_statuses = ["Treatment Not Started"],
      between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
      find_last_match_in_period=True,
      returning="date",
      date_format="YYYY-MM-DD",
      return_expectations={
          "incidence": 0.2,
          "date": {"earliest": "2021-12-16"},
      },
  )
def covid_therapeutics_stopped(dx_codelist):
  return patients.with_covid_therapeutics(
      dx_codelist,
      with_these_indications = "non_hospitalised",
      with_these_statuses = ["Treatment Stopped"],
      between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
      find_last_match_in_period=True,
      returning="date",
      date_format="YYYY-MM-DD",
      return_expectations={
          "incidence": 0.2,
          "date": {"earliest": "2021-12-16"},
      },
  )
def adverse_outcome_icd(dx_codelist):
  return patients.admitted_to_hospital(
      with_these_diagnoses = dx_codelist,
      returning="date_admitted",
      between = ["covid_test_positive_date", "covid_test_positive_date + 28 days"],
      find_first_match_in_period=True,
      date_format="YYYY-MM-DD",
      return_expectations={
          "incidence": 0.2,
          "date": {"earliest": "2021-12-16"},
      },
  )
def adverse_outcome_snomed(dx_codelist):
  return patients.with_these_clinical_events(
      dx_codelist,
      returning="date",
      between = ["covid_test_positive_date", "covid_test_positive_date + 28 days"],
      find_first_match_in_period=True,
      date_format="YYYY-MM-DD",
      return_expectations={
          "incidence": 0.2,
          "date": {"earliest": "2021-12-16"},
      },
  )
def imae_snomed(dx_codelist):  ## to idenitfy NEW conditions, will include all first diagnosis but then on stata remove those diagnosed before covid test positive date
  return patients.with_these_clinical_events(
      dx_codelist,
      returning="date",
      on_or_before ="covid_test_positive_date + 28 days",
      find_first_match_in_period=True,
      date_format="YYYY-MM-DD",
      return_expectations={
          "incidence": 0.2,
          "date": {"earliest": "2021-12-16"},
      },
  )
def imae_icd(dx_codelist):  ## to idenitfy NEW conditions, will include all first diagnosis but then on stata remove those diagnosed before covid test positive date
  return patients.admitted_to_hospital(
      with_these_diagnoses = dx_codelist,
      returning="date_admitted",
      on_or_before ="covid_test_positive_date + 28 days",
      find_first_match_in_period=True,
      date_format="YYYY-MM-DD",
      return_expectations={
          "incidence": 0.2,
          "date": {"earliest": "2021-12-16"},
      },
  )
study = StudyDefinition(
  # Configure the expectations framework
  default_expectations = {
    "date": {"earliest": "2021-11-01", "latest": "today"},
    "rate": "uniform",
    "incidence": 0.05,
  },
      
  ## Define index date
  index_date = campaign_start,

  population = patients.satisfying(
    """
    age >= 18 AND age < 110
    AND NOT has_died
    AND covid_test_positive
    AND eligble
    """,
  ),
        
  ## COVID TEST POSITIVE
  ## First positive SARS-CoV-2 test. [patients are eligible for treatment if diagnosed <=5d ago. Restricted to first positive covid test after index date]
  covid_test_positive = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2", test_result = "positive", returning = "binary_flag", on_or_after = "index_date - 5 days", find_first_match_in_period = True, 
    restrict_to_earliest_specimen_date = True,  return_expectations = {"incidence": 0.9 },
  ),
  covid_test_positive_date = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2", test_result = "positive", returning = "date", date_format = "YYYY-MM-DD", on_or_after = "index_date - 5 days", find_first_match_in_period = True, 
    restrict_to_earliest_specimen_date = True, return_expectations = {"date": {"earliest": "2021-12-11", "latest": "today"}, "incidence": 0.9},
  ),
  covid_positive_test_type = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2", test_result = "positive", returning = "case_category", on_or_after = "index_date - 5 days", find_first_match_in_period = True, 
    restrict_to_earliest_specimen_date = True, return_expectations = {"category": {"ratios": {"LFT_Only": 0.4, "PCR_Only": 0.4, "LFT_WithPCR": 0.2}},"incidence": 0.2, },
  ),
  ## Second positive SARS-CoV-2 test
  covid_test_positive_date2 = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2", test_result = "positive", returning = "date", date_format = "YYYY-MM-DD", on_or_after = "covid_test_positive_date + 30 days", find_first_match_in_period = True,
    restrict_to_earliest_specimen_date = True, return_expectations = {"date": {"earliest": "2021-12-20", "latest": "today"}, "incidence": 0.1 },
  ),
  ## Positive covid test 30 days prior to positive test [will only apply to patients who first tested positive towards the beginning of the study period]
  covid_positive_previous_30_days = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2", test_result = "positive", returning = "binary_flag",  between = ["covid_test_positive_date - 31 days", "covid_test_positive_date - 1 day"],
    find_last_match_in_period = True, restrict_to_earliest_specimen_date = False, return_expectations = {"incidence": 0.05 },
  ),
  ## Onset of symptoms of COVID-19
  symptomatic_covid_test = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2", test_result = "any", returning = "symptomatic", on_or_after = "index_date - 5 days", find_first_match_in_period = True,
    restrict_to_earliest_specimen_date = False, return_expectations={"incidence": 0.1, "category": {"ratios": {"": 0.2, "N": 0.2,"Y": 0.6, } }, },
  ),
  covid_symptoms_snomed = patients.with_these_clinical_events(
    covid_symptoms_snomed_codes, returning = "date", date_format = "YYYY-MM-DD", find_first_match_in_period = True, on_or_after = "index_date - 5 days",
  ),
  prior_covid = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2", test_result = "positive", returning = "binary_flag", on_or_before = "covid_test_positive_date - 1 day", find_last_match_in_period = True, 
    restrict_to_earliest_specimen_date = True,  return_expectations = {"incidence": 0.9 },
  ),
  prior_covid_date = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2", test_result = "positive", returning = "date", date_format = "YYYY-MM-DD", on_or_before = "covid_test_positive_date - 1 day", find_last_match_in_period = True, 
    restrict_to_earliest_specimen_date = True, return_expectations = {"date": {"earliest": "2021-12-11", "latest": "today"}, "incidence": 0.9},
  ),
  prior_covid_test_type = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2", test_result = "positive", returning = "case_category", on_or_before = "covid_test_positive_date - 1 day", find_last_match_in_period = True, 
    restrict_to_earliest_specimen_date = True, return_expectations = {"category": {"ratios": {"LFT_Only": 0.4, "PCR_Only": 0.4, "LFT_WithPCR": 0.2}},"incidence": 0.2, },
  ),
  ## SGTF indicator and Variant
  sgtf=patients.with_test_result_in_sgss(
      pathogen="SARS-CoV-2", test_result="positive", find_first_match_in_period=True, between=["covid_test_positive_date","covid_test_positive_date + 30 days"],
      returning="s_gene_target_failure", return_expectations={"rate": "universal", "category": {"ratios": {"0": 0.7, "1": 0.1, "9": 0.1, "": 0.1}}, },
  ), 
  ## New sgtf data in "all tests dataset" - not resitrcted to earliest speciman 
  sgtf_new=patients.with_test_result_in_sgss(
       pathogen="SARS-CoV-2", test_result="positive",find_first_match_in_period=True, between=["covid_test_positive_date","covid_test_positive_date + 30 days"],
       restrict_to_earliest_specimen_date=False, returning="s_gene_target_failure", return_expectations={"rate": "universal", "category": {"ratios": {"0": 0.7, "1": 0.1, "9": 0.1, "": 0.1}},},
  ), 
  variant_recorded=patients.with_test_result_in_sgss(
      pathogen="SARS-CoV-2", test_result="positive", find_first_match_in_period=True, between=["covid_test_positive_date","covid_test_positive_date + 30 days"],
      restrict_to_earliest_specimen_date=False, returning="variant", return_expectations={"rate": "universal", "category": {"ratios": {"B.1.617.2": 0.7, "VOC-21JAN-02": 0.2, "": 0.1}},},
  ), 
  

  # Eligable based on comorbidities 
  
  ## Down's syndrome
  downs_syndrome_snomed=comorbidity_snomed(downs_syndrome_nhsd_snomed_codes),
  downs_syndrome=comorbidity_icd(downs_syndrome_nhsd_icd10_codes),

  ## Solid cancer
  cancer_opensafely_snomed=comorbidity_snomed_6m(combine_codelists(non_haematological_cancer_opensafely_snomed_codes, lung_cancer_opensafely_snomed_codes, chemotherapy_radiotherapy_opensafely_snomed_codes)),
  
  ## Haematological diseases (malignancy taken as 12m instead of 24m)
  haematopoietic_stem_cell_snomed=comorbidity_snomed_12m(haematopoietic_stem_cell_transplant_nhsd_snomed_codes),
  haematopoietic_stem_cell_icd10=comorbidity_icd_12m(haematopoietic_stem_cell_transplant_nhsd_icd10_codes),
  haematological_malignancies_snomed=comorbidity_snomed_12m(haematological_malignancies_nhsd_snomed_codes),
  haematological_malignancies_icd10=comorbidity_icd_12m(haematological_malignancies_nhsd_icd10_codes),
  haematopoietic_stem_cell_opcs4=comorbidity_ops_12m(haematopoietic_stem_cell_transplant_nhsd_opcs4_codes),
  sickle_cell_disease_nhsd_snomed=comorbidity_snomed(sickle_cell_disease_nhsd_snomed_codes),
  sickle_cell_disease_nhsd_icd10=comorbidity_icd(sickle_cell_disease_nhsd_icd10_codes),
  haematological_disease_nhsd = patients.minimum_of("haematopoietic_stem_cell_snomed", 
                                                    "haematopoietic_stem_cell_icd10", 
                                                    "haematopoietic_stem_cell_opcs4", 
                                                    "haematological_malignancies_snomed", 
                                                    "haematological_malignancies_icd10",
                                                    "sickle_cell_disease_nhsd_snomed", 
                                                    "sickle_cell_disease_nhsd_icd10"), 
                                            
  ## Renal disease  (to get access to UK Renal Registry - advanced CKD stages 4 and 5, dialysis,  kidney transplantation in 2y care- need aproval)
  ckd_stage_5_nhsd_snomed=comorbidity_snomed(ckd_stage_5_nhsd_snomed_codes),
  ckd_stage_5_nhsd_icd10 =comorbidity_icd(ckd_stage_5_nhsd_icd10_codes),
  ckd_stage_5_nhsd = patients.minimum_of("ckd_stage_5_nhsd_snomed", "ckd_stage_5_nhsd_icd10"),

  ## Liver disease
  liver_disease_nhsd_snomed=comorbidity_snomed(liver_disease_nhsd_snomed_codes),
  liver_disease_nhsd_icd10=comorbidity_icd(liver_disease_nhsd_icd10_codes),
  # # advanced_decompensated_cirrhosis=comorbidity_snomed(advanced_decompensated_cirrhosis_snomed_codes),
  # # decompensated_cirrhosis_icd10=comorbidity_icd(advanced_decompensated_cirrhosis_icd10_codes),
  # # ascitic_drainage_snomed=comorbidity_snomed(ascitic_drainage_snomed_codes),

  liver_disease_nhsd = patients.minimum_of("liver_disease_nhsd_snomed", "liver_disease_nhsd_icd10"), 

  ## Immune-mediated inflammatory disorders (IMID) 
  rheumatoid_arthritis_nhsd_snomed=comorbidity_snomed(rheumatoid_arthritis_snowmed),
  rheumatoid_arthritis_nhsd_icd10=comorbidity_icd(rheumatoid_arthritis_icd10),
  SLE_nhsd_ctv=comorbidity_snomed(SLE_ctv),
  SLE_nhsd_icd10=comorbidity_icd(SLE_icd10), 
  Psoriasis_nhsd=comorbidity_snomed(Psoriasis_ctv3),
  Psoriatic_arthritis_nhsd=comorbidity_snomed(Psoriatic_arthritis_snomed),
  Ankylosing_Spondylitis_nhsd=comorbidity_snomed(Ankylosing_Spondylitis_ctv3),  
  IBD_nhsd=comorbidity_snomed(IBD_ctv3),
  
  ## Immunosuppression - Treatment steriods (4x prescriptions steroids in 6m or high dose) / CYC / MMF/ TAC / CIC
  immunosuppresant_drugs_nhsd=drug_6m(combine_codelists(immunosuppresant_drugs_dmd_codes, immunosuppresant_drugs_snomed_codes)),
  oral_steroid_drugs_nhsd=drug_12m(combine_codelists(oral_steroid_drugs_dmd_codes, oral_steroid_drugs_snomed_codes)),
  methotrexate_drugs_nhsd=drug_6m(combine_codelists(oral_methotrexate_drugs_snomed_codes, inj_methotrexate_drugs_snomed_codes)),
  ciclosporin_drugs_nhsd=drug_6m(oral_ciclosporin_snomed_codes),
  mycophenolate_drugs_nhsd=drug_6m(oral_mycophenolate_drugs_snomed_codes),
  abatacept_highcostdrugs = high_cost_drugs_3m(abatacept_high_cost_drugs_codes),
  adalimumab_highcostdrugs = high_cost_drugs_3m(adalimumab_high_cost_drugs_codes),
  alemtuzumab_highcostdrugs = high_cost_drugs_3m(alemtuzumab_high_cost_drugs_codes),
  baricitinib_highcostdrugs = high_cost_drugs_3m(baricitinib_high_cost_drugs_codes), 
  brodalumab_highcostdrugs = high_cost_drugs_3m(brodalumab_high_cost_drugs_codes),
  certolizumab_highcostdrugs = high_cost_drugs_3m(certolizumab_high_cost_drugs_codes),
  etanercept_highcostdrugs = high_cost_drugs_3m(etanercept_high_cost_drugs_codes),
  golimumab_highcostdrugs = high_cost_drugs_3m(golimumab_high_cost_drugs_codes),
  guselkumab_highcostdrugs = high_cost_drugs_3m(guselkumab_high_cost_drugs_codes),
  infliximab_highcostdrugs = high_cost_drugs_3m(infliximab_high_cost_drugs_codes),
  mepolizumab_highcostdrugs = high_cost_drugs_3m(mepolizumab_high_cost_drugs_codes),
  risankizumab_highcostdrugs = high_cost_drugs_3m(risankizumab_high_cost_drugs_codes),
  sarilumab_highcostdrugs = high_cost_drugs_3m(sarilumab_high_cost_drugs_codes),
  secukinumab_highcostdrugs = high_cost_drugs_3m(secukinumab_high_cost_drugs_codes),
  tildrakizumab_highcostdrugs = high_cost_drugs_3m(tildrakizumab_high_cost_drugs_codes),
  tocilizumab_highcostdrugs = high_cost_drugs_3m(tocilizumab_high_cost_drugs_codes),
  tofacitinib_highcostdrugs = high_cost_drugs_3m(tofacitinib_high_cost_drugs_codes),
  upadacitinib_highcostdrugs = high_cost_drugs_3m(upadacitinib_high_cost_drugs_codes),
  ustekinumab_highcostdrugs = high_cost_drugs_3m(ustekinumab_high_cost_drugs_codes),
  vedolizumab_highcostdrugs = high_cost_drugs_3m(vedolizumab_high_cost_drugs_codes),

  oral_steroid_drug_nhsd_6m_count = patients.with_these_medications(
    codelist = combine_codelists(oral_steroid_drugs_dmd_codes, oral_steroid_drugs_snomed_codes),
    returning = "number_of_matches_in_period",
    between = ["covid_test_positive_date - 6 months", "covid_test_positive_date"],
    return_expectations = {"incidence": 0.1,
      "int": {"distribution": "normal", "mean": 2, "stddev": 1},
    },
  ),
  
  imid_nhsd = patients.minimum_of("rheumatoid_arthritis_nhsd_snomed", "rheumatoid_arthritis_nhsd_icd10", "SLE_nhsd_ctv", "SLE_nhsd_icd10", "Psoriasis_nhsd", 
                                   "Psoriatic_arthritis_nhsd", "Ankylosing_Spondylitis_nhsd", "IBD_nhsd"), 
  imid_drug = patients.minimum_of("immunosuppresant_drugs_nhsd", "oral_steroid_drugs_nhsd", "methotrexate_drugs_nhsd", "ciclosporin_drugs_nhsd", "mycophenolate_drugs_nhsd"),
  imid_drug_HCD = patients.minimum_of("abatacept_highcostdrugs", "adalimumab_highcostdrugs", "alemtuzumab_highcostdrugs", "baricitinib_highcostdrugs", "brodalumab_highcostdrugs", 
                                   "certolizumab_highcostdrugs","etanercept_highcostdrugs", "golimumab_highcostdrugs",  "guselkumab_highcostdrugs", "infliximab_highcostdrugs",  
                                   "mepolizumab_highcostdrugs","risankizumab_highcostdrugs", "sarilumab_highcostdrugs", "secukinumab_highcostdrugs", "tildrakizumab_highcostdrugs",
                                   "tocilizumab_highcostdrugs", "tofacitinib_highcostdrugs", "upadacitinib_highcostdrugs","ustekinumab_highcostdrugs", "vedolizumab_highcostdrugs"),
  imid_on_drug = patients.satisfying(
    """
    imid_nhsd AND
    imid_drug OR imid_drug_HCD
    """,
  ),

  ## Primary immune deficiencies / HIV/AIDs
  immunosupression_nhsd=comorbidity_snomed(immunosupression_nhsd_codes),
  hiv_aids_nhsd_snomed=comorbidity_snomed(hiv_aids_nhsd_snomed_codes),
  hiv_aids_nhsd_icd10=comorbidity_icd(hiv_aids_nhsd_icd10_codes),
  hiv_aids_nhsd = patients.minimum_of("hiv_aids_nhsd_snomed", "hiv_aids_nhsd_icd10"),

  ## Solid organ transplant
  ## kb put in - solid_organ_transplant_nhsd_snomed_new, replacement_of_organ_transplant_nhsd_opcs4_codes
  ## kb removed - thymus_gland_transplant_nhsd_opcs4_codes, conjunctiva_transplant_nhsd_opcs4_codes 
  ## kb remove - stomach_transplant_nhsd_opcs4_codes transplant_stomach_opcs4 transplant_ileum_1_opcs4 which all used [between = ["transplant_all_y_codes_opcs4","transplant_all_y_codes_opcs4"],]
  solid_organ_transplant_nhsd_snomed=comorbidity_snomed(solid_organ_transplant_nhsd_snomed_codes),
  solid_organ_transplant_nhsd_snomed_new=comorbidity_snomed(solid_organ_transplant_nhsd_snomed_codes_new),
  solid_organ_transplant_nhsd_opcs4=comorbidity_ops(solid_organ_transplant_nhsd_opcs4_codes),
  solid_organ_replacement_nhsd_opcs4=comorbidity_ops(replacement_of_organ_transplant_nhsd_opcs4_codes),
  transplant_conjunctiva_y_code_opcs4=comorbidity_ops(conjunctiva_y_codes_transplant_nhsd_opcs4_codes),
  transplant_ileum_1_Y_codes_opcs4=comorbidity_ops(ileum_1_y_codes_transplant_nhsd_opcs4_codes),
  transplant_ileum_2_Y_codes_opcs4=comorbidity_ops(ileum_2_y_codes_transplant_nhsd_opcs4_codes), 
  solid_organ_transplant_nhsd = patients.minimum_of("solid_organ_transplant_nhsd_snomed", "solid_organ_transplant_nhsd_snomed_new", 
                                                    "solid_organ_transplant_nhsd_opcs4","solid_organ_replacement_nhsd_opcs4", "transplant_conjunctiva_y_code_opcs4", 
                                                    "transplant_ileum_1_Y_codes_opcs4","transplant_ileum_2_Y_codes_opcs4"),
 
  ### Neurological conditions
  multiple_sclerosis_nhsd_snomed=comorbidity_snomed(multiple_sclerosis_nhsd_snomed_codes),
  multiple_sclerosis_nhsd_icd10=comorbidity_icd(multiple_sclerosis_nhsd_icd10_codes),
  multiple_sclerosis_nhsd = patients.minimum_of("multiple_sclerosis_nhsd_snomed", "multiple_sclerosis_nhsd_icd10"), 
  motor_neurone_disease_nhsd_snomed=comorbidity_snomed(motor_neurone_disease_nhsd_snomed_codes),
  motor_neurone_disease_nhsd_icd10=comorbidity_icd(motor_neurone_disease_nhsd_icd10_codes),
  motor_neurone_disease_nhsd = patients.minimum_of("motor_neurone_disease_nhsd_snomed", "motor_neurone_disease_nhsd_icd10"),
  myasthenia_gravis_nhsd_snomed=comorbidity_snomed(myasthenia_gravis_nhsd_snomed_codes),
  myasthenia_gravis_nhsd_icd10=comorbidity_icd(myasthenia_gravis_nhsd_icd10_codes),
  myasthenia_gravis_nhsd = patients.minimum_of("myasthenia_gravis_nhsd_snomed", "myasthenia_gravis_nhsd_icd10"),
  huntingtons_disease_nhsd_snomed=comorbidity_snomed(huntingtons_disease_nhsd_snomed_codes),
  huntingtons_disease_nhsd_icd10=comorbidity_icd(huntingtons_disease_nhsd_icd10_codes),
  huntingtons_disease_nhsd = patients.minimum_of("huntingtons_disease_nhsd_snomed", "huntingtons_disease_nhsd_icd10"),
  
  ## Eligble 
  eligble = patients.minimum_of(
    "cancer_opensafely_snomed",
    "haematological_disease_nhsd",
    "ckd_stage_5_nhsd",
    "liver_disease_nhsd",
    "immunosupression_nhsd",
    "hiv_aids_nhsd",
    "solid_organ_transplant_nhsd",
    "multiple_sclerosis_nhsd",
    "motor_neurone_disease_nhsd",
    "myasthenia_gravis_nhsd",
    "huntingtons_disease_nhsd",
    ),
  
  ## TREATMENT - MAB + Antivirals. Either use on_or_after = "index date" or between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"]
  sotrovimab_covid_therapeutics=covid_therapeutics("Sotrovimab"),
  molnupiravir_covid_therapeutics =covid_therapeutics("Molnupiravir"),
  paxlovid_covid_therapeutics=covid_therapeutics("Paxlovid"),
  remdesivir_covid_therapeutics=covid_therapeutics("Remdesivir"),
  casirivimab_covid_therapeutics=covid_therapeutics("Casirivimab and imdevimab"),
  # # sotrovimab_covid_approved=covid_therapeutics_approved("Sotrovimab"),
  # # molnupiravir_covid_approved=covid_therapeutics_approved("Molnupiravir"),
  # # paxlovid_covid_approved=covid_therapeutics_approved("Paxlovid"),
  # # sotrovimab_covid_complete=covid_therapeutics_complete("Sotrovimab"), 
  # # molnupiravir_covid_complete =covid_therapeutics_complete("Molnupiravir"),
  # # paxlovid_covid_complete=covid_therapeutics_complete("Paxlovid"),
  # # sotrovimab_covid_stopped=covid_therapeutics_stopped("Sotrovimab"), 
  # # molnupiravir_covid_stopped =covid_therapeutics_stopped("Molnupiravir"),
  # # paxlovid_covid_stopped=covid_therapeutics_stopped("Paxlovid"),
  # # sotrovimab_covid_covid_not_start=covid_therapeutics_could_not_start("Sotrovimab"), 
  # # molnupiravir_covid_covid_not_start=covid_therapeutics_could_not_start("Molnupiravir"),
  # # paxlovid_covid_covid_not_start=covid_therapeutics_could_not_start("Paxlovid"),
  
  ## Date treated
  date_treated = patients.minimum_of(
    "sotrovimab_covid_therapeutics","paxlovid_covid_therapeutics", "molnupiravir_covid_therapeutics",
  ),

  # DEMOGRAPHIC COVARIATES & COMORBIDITY
  age = patients.age_as_of(
    "covid_test_positive_date",
    return_expectations = {
      "rate": "universal",
      "int": {"distribution": "population_ages"},
      "incidence" : 0.9
    },
  ),
  sex = patients.sex(
    return_expectations = {
      "rate": "universal",
      "category": {"ratios": {"M": 0.49, "F": 0.51}},
    }
  ),
  ethnicity = patients.categorised_as(
            {"Missing": "DEFAULT",
            "White": "eth='1' OR (NOT eth AND ethnicity_sus='1')", 
            "Mixed": "eth='2' OR (NOT eth AND ethnicity_sus='2')", 
            "South Asian": "eth='3' OR (NOT eth AND ethnicity_sus='3')", 
            "Black": "eth='4' OR (NOT eth AND ethnicity_sus='4')",  
            "Other": "eth='5' OR (NOT eth AND ethnicity_sus='5')",
            }, 
            return_expectations={
            "category": {"ratios": {"White": 0.6, "Mixed": 0.1, "South Asian": 0.1, "Black": 0.1, "Other": 0.1}},
            "incidence": 0.4,
            },
            ethnicity_sus = patients.with_ethnicity_from_sus(
                returning="group_6",  
                use_most_frequent_code=True,
                return_expectations={
                    "category": {"ratios": {"1": 0.6, "2": 0.1, "3": 0.1, "4": 0.1, "5": 0.1}},
                    "incidence": 0.4,
                    },
            ),
            eth=patients.with_these_clinical_events(
                ethnicity_primis_snomed_codes,
                returning="category",
                find_last_match_in_period=True,
                on_or_before="today",
                return_expectations={
                    "category": {"ratios": {"1": 0.6, "2": 0.1, "3": 0.1, "4":0.1,"5": 0.1}},
                    "incidence": 0.75,
                },
            ),
    ),
  imd = patients.categorised_as(
    {     "0": "DEFAULT",
          "1": "index_of_multiple_deprivation >= 0 AND index_of_multiple_deprivation < 32800*1/5",
          "2": "index_of_multiple_deprivation >= 32800*1/5 AND index_of_multiple_deprivation < 32800*2/5",
          "3": "index_of_multiple_deprivation >= 32800*2/5 AND index_of_multiple_deprivation < 32800*3/5",
          "4": "index_of_multiple_deprivation >= 32800*3/5 AND index_of_multiple_deprivation < 32800*4/5",
          "5": "index_of_multiple_deprivation >= 32800*4/5 AND index_of_multiple_deprivation <= 32800",
    },
    index_of_multiple_deprivation = patients.address_as_of(
      "covid_test_positive_date",
      returning = "index_of_multiple_deprivation",
      round_to_nearest = 100,
    ),
    return_expectations = {
      "rate": "universal",
      "category": {
        "ratios": {"0": 0.01,"1": 0.20,"2": 0.20, "3": 0.20, "4": 0.20,"5": 0.19,}
      },
    },
  ),
  region_covid_therapeutics = patients.with_covid_therapeutics(
    with_these_therapeutics = ["Sotrovimab", "Molnupiravir", "Paxlovid"],
    with_these_indications = "non_hospitalised",
    on_or_after = "covid_test_positive_date",
    find_first_match_in_period = True,
    returning = "region",
    return_expectations = {
      "rate": "universal",
      "category": {
        "ratios": {
          "North East": 0.1, "North West": 0.1,"Yorkshire and The Humber": 0.1, "East Midlands": 0.1,"West Midlands": 0.1, "East": 0.1, "London": 0.2, "South West": 0.1, "South East": 0.1,},
      },
    },
  ),

  # STP [NHS administration region based on geography, currently closest match to CMDU] 
  stp = patients.registered_practice_as_of(
    "covid_test_positive_date",
    returning = "stp_code",
    return_expectations = {
      "rate": "universal",
      "category": {
        "ratios": {
          "STP1": 0.1, "STP2": 0.1,"STP3": 0.1,"STP4": 0.1,"STP5": 0.1,"STP6": 0.1,"STP7": 0.1,"STP8": 0.1, "STP9": 0.1, "STP10": 0.1, }
      },
    },
  ),

  vaccination_status = patients.categorised_as(
    { "Un-vaccinated": "DEFAULT",
      "Un-vaccinated (declined)": """ covid_vax_declined AND NOT (covid_vax_1 OR covid_vax_2 OR covid_vax_3)""",
      "One vaccination": """ covid_vax_1 AND NOT covid_vax_2 """,
      "Two vaccinations": """ covid_vax_2 AND NOT covid_vax_3 """,
      "Three or more vaccinations": """ covid_vax_3 """
    },
    covid_vax_1 = patients.with_tpp_vaccination_record(
      target_disease_matches = "SARS-2 CORONAVIRUS",
      between = ["2020-06-08", "covid_test_positive_date"],
      find_first_match_in_period = True,
      returning = "date",
      date_format = "YYYY-MM-DD"
    ),
    covid_vax_2 = patients.with_tpp_vaccination_record(
      target_disease_matches = "SARS-2 CORONAVIRUS",
      between = ["covid_vax_1 + 19 days", "covid_test_positive_date"],
      find_first_match_in_period = True,
      returning = "date",
      date_format = "YYYY-MM-DD"
    ),
    covid_vax_3 = patients.with_tpp_vaccination_record(
      target_disease_matches = "SARS-2 CORONAVIRUS",
      between = ["covid_vax_2 + 56 days", "covid_test_positive_date"],
      find_first_match_in_period = True,
      returning = "date",
      date_format = "YYYY-MM-DD"
    ),
   covid_vax_declined = patients.with_these_clinical_events(
      covid_vaccine_declined_codes,
      returning="binary_flag",
      on_or_before = "covid_test_positive_date",
    ),
    return_expectations = {
      "rate": "universal",
      "category": {
        "ratios": {
          "Un-vaccinated": 0.1,
          "Un-vaccinated (declined)": 0.1,
          "One vaccination": 0.1,
          "Two vaccinations": 0.2,
          "Three or more vaccinations": 0.5,
        }
      },
    },
  ),
  last_vaccination_date = patients.with_tpp_vaccination_record(
      target_disease_matches = "SARS-2 CORONAVIRUS",
      on_or_before = "covid_test_positive_date",
      find_last_match_in_period = True,
      returning = "date",
      date_format = "YYYY-MM-DD",
      return_expectations={
            "date": {"earliest": "2020-06-08", "latest": "today"},
            "incidence": 0.95,
      }
  ),

  bmi=patients.most_recent_bmi(
      on_or_before="covid_test_positive_date",
      minimum_age_at_measurement=18,
      include_measurement_date=True,
      date_format="YYYY-MM-DD",
      return_expectations={
            "date": {"earliest": "2020-01-01", "latest": "today"},
            "float": {"distribution": "normal", "mean": 28, "stddev": 8},
            "incidence": 0.95,
        }
  ),
  ## to return date, so to make binary in stata
  diabetes=comorbidity_snomed(diabetes_codes),
  chronic_cardiac_disease=comorbidity_snomed(chronic_cardiac_dis_codes),
  hypertension=comorbidity_snomed(hypertension_codes),
  chronic_respiratory_disease=comorbidity_snomed(chronic_respiratory_dis_codes),


  ### EXCLUSIONS 
  ## Died
  has_died = patients.died_from_any_cause(
    on_or_before = "covid_test_positive_date - 1 day",
    returning = "binary_flag",
  ),
  death_date = patients.died_from_any_cause(
    returning = "date_of_death",
    date_format = "YYYY-MM-DD",
    on_or_after = "covid_test_positive_date - 1 day",
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date"},
      "incidence": 0.1
    },
  ),  
  
  ## De-registration
  dereg_date = patients.date_deregistered_from_all_supported_practices(
    on_or_after = "covid_test_positive_date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date"},
      "incidence": 0.1
    },
  ), 
  
  ## Hospitalised with covid between positive test and 5 days
  primary_covid_hospital_discharge_date = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], 
    between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
    date_format = "YYYY-MM-DD",
    find_first_match_in_period = False,
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date - 1 day"},
      "rate": "uniform",
      "incidence": 0.05
    },
  ),
  any_covid_hospital_discharge_date = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_these_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], 
    with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], 
    between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
    date_format = "YYYY-MM-DD",
    find_first_match_in_period = False,
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date - 1 day"},
      "rate": "uniform",
      "incidence": 0.05
    },
  ),
  
  #  pregnancy record in last 36 weeks.  pregnancy OR delivery code since latest pregnancy record:
  preg_36wks_date = patients.with_these_clinical_events(
    pregnancy_primis_codes,
    returning = "date",
    find_last_match_in_period = True,
    between = ["covid_test_positive_date - 252 days", "covid_test_positive_date - 1 day"],
    date_format = "YYYY-MM-DD",
  ),
  pregdel = patients.with_these_clinical_events(
    pregdel_primis_codes,
    returning = "binary_flag", find_last_match_in_period = True, between = ["preg_36wks_date + 1 day", "covid_test_positive_date - 1 day"], date_format = "YYYY-MM-DD",
  ),
  pregnancy = patients.satisfying(
    """
    gender = 'F' AND preg_age <= 50
    AND (preg_36wks_date AND NOT pregdel)
    """,
    
    gender = patients.sex(return_expectations = {"rate": "universal", "category": {"ratios": {"M": 0.49, "F": 0.51}},}
    ),
    preg_age = patients.age_as_of("preg_36wks_date", return_expectations = {"rate": "universal","int": {"distribution": "population_ages"}, "incidence" : 0.9 },
    ),
  ),

  ## OUTCOME 
  ## Time to outcome - coded in stata [28 days after +covid or 28day after Tx]
  ## 1) DRUG REACTIONS AND AESI (diarrhoea, diverticulitis, altered taste)
  ae_diverticulitis_icd=adverse_outcome_icd(diverticulitis_icd_codes),
  ae_diverticulitis_snomed=adverse_outcome_snomed(diverticulitis_snomed_codes),
  ae_diarrhoea_snomed=adverse_outcome_snomed(diarrhoea_snomed_codes), 
  ae_taste_snomed=adverse_outcome_snomed(taste_snomed_codes),
  ae_taste_icd=adverse_outcome_icd(taste_icd_codes),
  ae_rheumatoid_arthritis_snomed=imae_snomed(rheumatoid_arthritis_snowmed),
  ae_rheumatoid_arthritis_icd=imae_icd(rheumatoid_arthritis_icd10),
  ae_SLE_ctv=imae_snomed(SLE_ctv),
  ae_SLE_icd=imae_icd(SLE_icd10), 
  ae_Psoriasis_snomed=imae_snomed(Psoriasis_ctv3),
  ae_Psoriatic_arthritis_snomed=imae_snomed(Psoriatic_arthritis_snomed),
  ae_Ankylosing_Spondylitis_ctv=imae_snomed(Ankylosing_Spondylitis_ctv3),  
  ae_IBD_snomed=imae_snomed(IBD_ctv3),

  ## 2) ALL SAE INCLUDING COVID [Note need to consider patients admitted for MAB infusion]
  all_hosp_date = patients.admitted_to_hospital(
    returning = "date_admitted", with_patient_classification = ["1"], with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], 
    between = ["covid_test_positive_date", "covid_test_positive_date + 28 day"], find_first_match_in_period = True, date_format = "YYYY-MM-DD",
    return_expectations = {"date": {"earliest": "2022-02-16"},"rate": "uniform", "incidence": 0.1},
  ), 
  ## all_hosp_cause = patients.admitted_to_hospital(returning = "primary_diagnosis", with_patient_classification = ["1"], between = ["covid_test_positive_date", "covid_test_positive_date + 28 day"], find_first_match_in_period = True, return_expectations=list, ), 

  covid_hosp_date = patients.admitted_to_hospital(
    returning = "date_admitted", with_these_diagnoses = covid_icd10_codes, with_patient_classification = ["1"], with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], 
    between = ["covid_test_positive_date", "covid_test_positive_date + 28 day"], find_first_match_in_period = True, date_format = "YYYY-MM-DD",
    return_expectations = {"date": {"earliest": "2022-02-16"}, "rate": "uniform", "incidence": 0.1 },
  ),  
  emerg_covid_hosp_date = patients.admitted_to_hospital(
    returning = "date_admitted", with_these_diagnoses = covid_icd10_codes, with_patient_classification = ["1"], with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], 
    between = ["covid_test_positive_date", "covid_test_positive_date + 28 day"], find_first_match_in_period = True, date_format = "YYYY-MM-DD",
    return_expectations = {"date": {"earliest": "2022-02-16"}, "rate": "uniform", "incidence": 0.1 },
  ),  
  covid_hosp_date_mabs_procedure = patients.admitted_to_hospital(
    returning = "date_admitted", with_patient_classification = ["1"], with_these_procedures=mabs_procedure_codes, 
    between = ["covid_test_positive_date", "covid_test_positive_date + 28 day"], find_first_match_in_period = True, date_format = "YYYY-MM-DD",
    return_expectations = {"date": {"earliest": "2022-02-16"},"rate": "uniform", "incidence": 0.1},
  ),
  died_date_ons=patients.died_from_any_cause(
    returning="date_of_death", 
    between = ["covid_test_positive_date", "covid_test_positive_date + 28 day"], date_format = "YYYY-MM-DD",
    return_expectations = {"date": {"earliest": "2021-01-01", "latest" : "today"}, "rate": "uniform","incidence": 0.6},
  ),  
  # # death_code = patients.died_from_any_cause(returning = "underlying_cause_of_death", between = ["covid_test_positive_date", "covid_test_positive_date + 28 day"]),

  died_ons_covid=patients.with_these_codes_on_death_certificate(
    covid_icd10_codes, returning="date_of_death",       
    between = ["covid_test_positive_date", "covid_test_positive_date + 28 day"], date_format = "YYYY-MM-DD",  match_only_underlying_cause=True,
    return_expectations = {"date": {"earliest": "2021-01-01", "latest" : "today"}, "rate": "uniform", "incidence": 0.6},
  ),
)