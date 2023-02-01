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
          "date": {"earliest": "2020-02-01", "latest": end_date},
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
          "date": {"earliest": "2020-02-01", "latest": end_date},
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
          "date": {"earliest": "2020-02-01", "latest": end_date},
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
          "date": {"earliest": "2020-02-01", "latest": end_date},
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
          "date": {"earliest": "2020-02-01", "latest": end_date},
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
          "date": {"earliest": "2020-02-01", "latest": end_date},
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
        
  ## NOT died 
  has_died = patients.died_from_any_cause(
    on_or_before = "index_date - 1 day",
    returning = "binary_flag",
  ),

  ### First positive SARS-CoV-2 test
  # Note patients are eligible for treatment if diagnosed <=5d ago. Restricted to first positive covid test after index date
  covid_test_positive = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    returning = "binary_flag",
    on_or_after = "index_date - 5 days",
    find_first_match_in_period = True,
    restrict_to_earliest_specimen_date = True,
    return_expectations = {
      "incidence": 0.9
    },
  ),

  covid_test_positive_date = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    find_first_match_in_period = True,
    restrict_to_earliest_specimen_date = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    on_or_after = "index_date - 5 days",
    return_expectations = {
      "date": {"earliest": "2021-12-11", "latest": "today"},
      "incidence": 0.9
    },
  ),

  ### Second positive SARS-CoV-2 test
  covid_test_positive_date2 = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    find_first_match_in_period = True,
    restrict_to_earliest_specimen_date = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    on_or_after = "covid_test_positive_date + 30 days",
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "today"},
      "incidence": 0.1
    },
  ),

  ## Covid test type - whether PCR or lat flow (possibly not available) 
  covid_positive_test_type = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    returning = "case_category",
    on_or_after = "index_date - 5 days",
    restrict_to_earliest_specimen_date = True,
    return_expectations = {
      "category": {"ratios": {"LFT_Only": 0.4, "PCR_Only": 0.4, "LFT_WithPCR": 0.2}},
      "incidence": 0.2,
    },
  ),

  ### Positive covid test 30 days prior to positive test
  # (note this will only apply to patients who first tested positive towards the beginning of the study period)
  covid_positive_previous_30_days = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    returning = "binary_flag",
    between = ["covid_test_positive_date - 31 days", "covid_test_positive_date - 1 day"],
    find_last_match_in_period = True,
    restrict_to_earliest_specimen_date = False,
    return_expectations = {
      "incidence": 0.05
    },
  ),

  ### Onset of symptoms of COVID-19
  symptomatic_covid_test = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "any",
    returning = "symptomatic",
    on_or_after = "index_date - 5 days",
    find_first_match_in_period = True,
    restrict_to_earliest_specimen_date = False,
    return_expectations={
      "incidence": 0.1,
      "category": {
        "ratios": {
          "": 0.2,
          "N": 0.2,
          "Y": 0.6,
        }
      },
    },
  ),

  covid_symptoms_snomed = patients.with_these_clinical_events(
    covid_symptoms_snomed_codes,
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_first_match_in_period = True,
    on_or_after = "index_date - 5 days",
  ),

  ### SGTF indicator and Variant
  sgtf=patients.with_test_result_in_sgss(
      pathogen="SARS-CoV-2",
      test_result="positive",
      find_first_match_in_period=True,
      between=["covid_test_positive_date","covid_test_positive_date + 30 days"],
      returning="s_gene_target_failure",
      return_expectations={
          "rate": "universal",
          "category": {"ratios": {"0": 0.7, "1": 0.1, "9": 0.1, "": 0.1}},
      },
  ), 
  
  # new sgtf data in "all tests dataset" - not resitrcted to earliest speciman 
  sgtf_new=patients.with_test_result_in_sgss(
       pathogen="SARS-CoV-2",
       test_result="positive",
       find_first_match_in_period=True,
       restrict_to_earliest_specimen_date=False,
       between=["covid_test_positive_date","covid_test_positive_date + 30 days"],
       returning="s_gene_target_failure",
       return_expectations={
            "rate": "universal",
            "category": {"ratios": {"0": 0.7, "1": 0.1, "9": 0.1, "": 0.1}},
       },
  ), 
  variant_recorded=patients.with_test_result_in_sgss(
      pathogen="SARS-CoV-2",
      test_result="positive",
      find_first_match_in_period=True,
      restrict_to_earliest_specimen_date=False,
      between=["covid_test_positive_date","covid_test_positive_date + 30 days"],
      returning="variant",
      return_expectations={
          "rate": "universal",
          "category": {"ratios": {"B.1.617.2": 0.7, "VOC-21JAN-02": 0.2, "": 0.1}},
      },
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
  IMID_nhsd = patients.minimum_of("rheumatoid_arthritis_nhsd_snomed", "rheumatoid_arthritis_nhsd_icd10", "SLE_nhsd_ctv", "SLE_nhsd_icd10", "Psoriasis_nhsd", "Psoriatic_arthritis_nhsd", "Ankylosing_Spondylitis_nhsd", "IBD_nhsd"), 
  
  ## Immunosuppression - Treatment steriods (4x prescriptions steroids in 6m or high dose) / CYC / MMF/ TAC / CIC
  immunosuppresant_drugs_nhsd=drug_6m(combine_codelists(immunosuppresant_drugs_dmd_codes, immunosuppresant_drugs_snomed_codes)),
  oral_steroid_drugs_nhsd=drug_12m(combine_codelists(oral_steroid_drugs_dmd_codes, oral_steroid_drugs_snomed_codes)),
  methotrexate_drugs_nhsd=drug_6m(combine_codelists(oral_methotrexate_drugs_snomed_codes, inj_methotrexate_drugs_snomed_codes)),
  ciclosporin_drugs_nhsd=drug_6m(oral_ciclosporin_snomed_codes),
  mycophenolate_drugs_nhsd=drug_6m(oral_mycophenolate_drugs_snomed_codes),

  oral_steroid_drug_nhsd_6m_count = patients.with_these_medications(
    codelist = combine_codelists(oral_steroid_drugs_dmd_codes, oral_steroid_drugs_snomed_codes),
    returning = "number_of_matches_in_period",
    between = ["covid_test_positive_date - 6 months", "covid_test_positive_date"],
    return_expectations = {"incidence": 0.1,
      "int": {"distribution": "normal", "mean": 2, "stddev": 1},
    },
  ),

  imid_on_drug_nhsd = patients.minimum_of("immunosuppresant_drugs_nhsd", "oral_steroid_drugs_nhsd"), 

  # # Where are the high cost drugs

  # # **medication_dates("methotrexate", "opensafely-methotrexate-oral", False, True),
  # # **medication_dates("methotrexate_inj", "opensafely-methotrexate-injectable", False, True),
  # # **medication_dates("mycophenolate", "opensafely-mycophenolate", False, True),
  # # **medication_dates("ciclosporin", "opensafely-ciclosporin-oral-dmd", False, True),
  # # **medication_dates("abatacept", "opensafely-high-cost-drugs-abatacept", True, False),
  # # **medication_dates("adalimumab", "opensafely-high-cost-drugs-adalimumab", True, False),
  # # **medication_dates("baricitinib", "opensafely-high-cost-drugs-baricitinib", True, False),
  # # **medication_dates("certolizumab", "opensafely-high-cost-drugs-certolizumab", True, False),
  # # **medication_dates("etanercept", "opensafely-high-cost-drugs-etanercept", True, False),
  # # **medication_dates("golimumab", "opensafely-high-cost-drugs-golimumab", True, False),
  # # **medication_dates("guselkumab", "opensafely-high-cost-drugs-guselkumab", True, False),
  # # **medication_dates("infliximab", "opensafely-high-cost-drugs-infliximab", True, False),
  # # **medication_dates("ixekizumab", "opensafely-high-cost-drugs-ixekizumab", True, False),
  # # **medication_dates("methotrexate_hcd", "opensafely-high-cost-drugs-methotrexate", True, False),
  # # **medication_dates("rituximab", "opensafely-high-cost-drugs-rituximab", True, False),
  # # **medication_dates("sarilumab", "opensafely-high-cost-drugs-sarilumab", True, False),
  # # **medication_dates("secukinumab", "opensafely-high-cost-drugs-secukinumab", True, False),
  # # **medication_dates("tocilizumab", "opensafely-high-cost-drugs-tocilizumab", True, False),
  # # **medication_dates("tofacitinib", "opensafely-high-cost-drugs-tofacitinib", True, False),
  # # **medication_dates("upadacitinib", "opensafely-high-cost-drugs-upadacitinib", True, False),
  # # **medication_dates("ustekinumab", "opensafely-high-cost-drugs-ustekinumab", True, False),

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
    "imid_on_drug_nhsd",
    "immunosupression_nhsd",
    "hiv_aids_nhsd",
    "solid_organ_transplant_nhsd",
    "multiple_sclerosis_nhsd",
    "motor_neurone_disease_nhsd",
    "myasthenia_gravis_nhsd",
    "huntingtons_disease_nhsd",
    ),
  

  ## TREATMENT - MAB + Antivirals. Either use on_or_after = "index date" or between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"]
  
  ## Sotrovimab
  sotrovimab_covid_therapeutics = patients.with_covid_therapeutics(
    #with_these_statuses = ["Approved", "Treatment Complete"],
    with_these_therapeutics = "Sotrovimab",
    with_these_indications = "non_hospitalised",
    between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.2
    },
  ),
  # restrict by status
  sotrovimab_covid_approved = patients.with_covid_therapeutics(
    with_these_statuses = ["Approved"],
    with_these_therapeutics = "Sotrovimab",
    with_these_indications = "non_hospitalised",
    between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.2
    },
  ),
  sotrovimab_covid_complete = patients.with_covid_therapeutics(
    with_these_statuses = ["Treatment Complete"],
    with_these_therapeutics = "Sotrovimab",
    with_these_indications = "non_hospitalised",
    between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.2
    },
  ),
  sotrovimab_covid_not_start = patients.with_covid_therapeutics(
    with_these_statuses = ["Treatment Not Started"],
    with_these_therapeutics = "Sotrovimab",
    with_these_indications = "non_hospitalised",
    between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.2
    },
  ),
  sotrovimab_covid_stopped = patients.with_covid_therapeutics(
    with_these_statuses = ["Treatment Stopped"],
    with_these_therapeutics = "Sotrovimab",
    with_these_indications = "non_hospitalised",
    between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.4
    },
  ),
  ### Molnupiravir
  molnupiravir_covid_therapeutics = patients.with_covid_therapeutics(
    #with_these_statuses = ["Approved", "Treatment Complete"],
    with_these_therapeutics = "Molnupiravir",
    with_these_indications = "non_hospitalised",
    between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.2
    },
  ),
  # restrict by status
  molnupiravir_covid_approved = patients.with_covid_therapeutics(
    with_these_statuses = ["Approved"],
    with_these_therapeutics = "Molnupiravir",
    with_these_indications = "non_hospitalised",
    between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.2
    },
  ),
  molnupiravir_covid_complete = patients.with_covid_therapeutics(
    with_these_statuses = ["Treatment Complete"],
    with_these_therapeutics = "Molnupiravir",
    with_these_indications = "non_hospitalised",
    between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.2
    },
  ),
  molnupiravir_covid_not_start = patients.with_covid_therapeutics(
    with_these_statuses = ["Treatment Not Started"],
    with_these_therapeutics = "Molnupiravir",
    with_these_indications = "non_hospitalised",
    between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.2
    },
  ),
  molnupiravir_covid_stopped = patients.with_covid_therapeutics(
    with_these_statuses = ["Treatment Stopped"],
    with_these_therapeutics = "Molnupiravir",
    with_these_indications = "non_hospitalised",
    between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.2
    },
  ),
  ### Paxlovid
  paxlovid_covid_therapeutics = patients.with_covid_therapeutics(
    with_these_therapeutics = "Paxlovid",
    with_these_indications = "non_hospitalised",
    between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.2
    },
  ),
  # restrict by status
  paxlovid_covid_approved = patients.with_covid_therapeutics(
    with_these_statuses = ["Approved"],
    with_these_therapeutics = "Paxlovid",
    with_these_indications = "non_hospitalised",
    between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.2
    },
  ),
  paxlovid_covid_complete = patients.with_covid_therapeutics(
    with_these_statuses = ["Treatment Complete"],
    with_these_therapeutics = "Paxlovid",
    with_these_indications = "non_hospitalised",
    between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.2
    },
  ),
  paxlovid_covid_not_start = patients.with_covid_therapeutics(
    with_these_statuses = ["Treatment Not Started"],
    with_these_therapeutics = "Paxlovid",
    with_these_indications = "non_hospitalised",
    between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.2
    },
  ),
  paxlovid_covid_stopped = patients.with_covid_therapeutics(
    with_these_statuses = ["Treatment Stopped"],
    with_these_therapeutics = "Paxlovid",
    with_these_indications = "non_hospitalised",
    between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.2
    },
  ),
  ## Remdesivir
  remdesivir_covid_therapeutics = patients.with_covid_therapeutics(
    # with_these_statuses = ["Approved", "Treatment Complete"],
    with_these_therapeutics = "Remdesivir",
    with_these_indications = "non_hospitalised",
    on_or_after = "covid_test_positive_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.02
    },
  ),
    
  ### Casirivimab and imdevimab
  casirivimab_covid_therapeutics = patients.with_covid_therapeutics(
    # with_these_statuses = ["Approved", "Treatment Complete"],
    with_these_therapeutics = "Casirivimab and imdevimab",
    with_these_indications = "non_hospitalised",
    on_or_after = "covid_test_positive_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.02
    },
  ), 
    
  ## Date treated
  date_treated = patients.minimum_of(
    "sotrovimab_covid_therapeutics",
    "paxlovid_covid_therapeutics",
    "molnupiravir_covid_therapeutics",
  ),

  ### Exclusion 
  # Require hospitalisation for COVID-19 (note this data lags behind the therapeutics/testing data so may be missing)
  # HELP - how to exclude admitted for covid covid_test_positive_date & covid_test_positive_date + 5 days  
  primary_covid_hospital_discharge_date = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders. see https://docs.opensafely.org/study-def-variables/#sus for more info
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
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    between = ["covid_test_positive_date", "covid_test_positive_date + 5 days"],
    date_format = "YYYY-MM-DD",
    find_first_match_in_period = False,
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date - 1 day"},
      "rate": "uniform",
      "incidence": 0.05
    },
  ),

  # CLINICAL/DEMOGRAPHIC COVARIATES ----

  ## Age 
  age = patients.age_as_of(
    "covid_test_positive_date",
    return_expectations = {
      "rate": "universal",
      "int": {"distribution": "population_ages"},
      "incidence" : 0.9
    },
  ),
  ## Sex
  sex = patients.sex(
    return_expectations = {
      "rate": "universal",
      "category": {"ratios": {"M": 0.49, "F": 0.51}},
    }
  ),
  ## Ethnicity
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
  ## Index of multiple deprivation
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
        "ratios": {
          "0": 0.01,
          "1": 0.20,
          "2": 0.20,
          "3": 0.20,
          "4": 0.20,
          "5": 0.19,
        }
      },
    },
  ),
  
  region_covid_therapeutics = patients.with_covid_therapeutics(
    #with_these_statuses = ["Approved", "Treatment Complete"],
    with_these_therapeutics = ["Sotrovimab", "Molnupiravir", "Paxlovid"],
    with_these_indications = "non_hospitalised",
    on_or_after = "covid_test_positive_date",
    find_first_match_in_period = True,
    returning = "region",
    return_expectations = {
      "rate": "universal",
      "category": {
        "ratios": {
          "North East": 0.1,
          "North West": 0.1,
          "Yorkshire and The Humber": 0.1,
          "East Midlands": 0.1,
          "West Midlands": 0.1,
          "East": 0.1,
          "London": 0.2,
          "South West": 0.1,
          "South East": 0.1,},},
    },
  ),

 # STP (NHS administration region based on geography, currently closest match to CMDU) 
 stp = patients.registered_practice_as_of(
    "covid_test_positive_date",
    returning = "stp_code",
    return_expectations = {
      "rate": "universal",
      "category": {
        "ratios": {
          "STP1": 0.1,
          "STP2": 0.1,
          "STP3": 0.1,
          "STP4": 0.1,
          "STP5": 0.1,
          "STP6": 0.1,
          "STP7": 0.1,
          "STP8": 0.1,
          "STP9": 0.1,
          "STP10": 0.1,
        }
      },
    },
  ),

  ## Vaccination status
  vaccination_status = patients.categorised_as(
    { "Un-vaccinated": "DEFAULT",
      "Un-vaccinated (declined)": """ covid_vax_declined AND NOT (covid_vax_1 OR covid_vax_2 OR covid_vax_3)""",
      "One vaccination": """ covid_vax_1 AND NOT covid_vax_2 """,
      "Two vaccinations": """ covid_vax_2 AND NOT covid_vax_3 """,
      "Three or more vaccinations": """ covid_vax_3 """
    },
    
    # first vaccine from during trials and up to covid infection date
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

  # latest vaccination date
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

  # CLINICAL CO-MORBIDITIES  ----
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
  diabetes=patients.with_these_clinical_events(
        diabetes_codes,
        on_or_before="covid_test_positive_date",
        returning="binary_flag",
        return_expectations={"incidence": 0.1, },
  ),
  chronic_cardiac_disease=patients.with_these_clinical_events(
        chronic_cardiac_dis_codes,
        on_or_before="covid_test_positive_date",
        returning="binary_flag",
        return_expectations={"incidence": 0.1, },
  ),
  hypertension=patients.with_these_clinical_events(
        hypertension_codes,
        on_or_before="covid_test_positive_date",
        returning="binary_flag",
        return_expectations={"incidence": 0.1, },
  ),
  chronic_respiratory_disease=patients.with_these_clinical_events(
        chronic_respiratory_dis_codes,
        on_or_before="covid_test_positive_date",
        returning="binary_flag",
        return_expectations={"incidence": 0.1, },
  ),
)
