
# --- IMPORT STATEMENTS ---

## Import code building blocks from cohort extractor package
from cohortextractor import (codelist, codelist_from_csv, combine_codelists)


# --- CODELISTS ---
# SUS-HES mabs
mabs_procedure_codes = codelist(
  ["X891", "X892"], system="opcs4"
)
# advanced decompensated liver cirrhosis
advanced_decompensated_cirrhosis_snomed_codes = codelist_from_csv(
    "codelists/opensafely-condition-advanced-decompensated-cirrhosis-of-the-liver.csv",
    system="snomed",
    column="code"
)
advanced_decompensated_cirrhosis_icd10_codes = codelist_from_csv(
    "codelists/opensafely-condition-advanced-decompensated-cirrhosis-of-the-liver-and-associated-conditions-icd-10.csv",
    system="icd10",
    column="code"
)
# ascitic drainage
ascitic_drainage_snomed_codes = codelist_from_csv(
    "codelists/opensafely-procedure-ascitic-drainage.csv",
    system="snomed",
    column="code"
)
# CKD 3-5
chronic_kidney_disease_stages_3_5_codes = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-ckd35.csv",
    system="snomed",
    column="code",
)
primis_ckd_stage = codelist_from_csv(
    "codelists/user-Louis-ckd-stage.csv",
    system="snomed",
    column="code",
    category_column="stage"
)
# dialysis
dialysis_codes = codelist_from_csv(
  "codelists/opensafely-dialysis.csv",
  system = "ctv3",
  column = "CTV3ID"
)
dialysis_icd10_codelist = codelist_from_csv(
    "codelists/ukrr-dialysis.csv",
    system="icd10",
    column="code"
)
dialysis_opcs4_codelist = codelist_from_csv(
    "codelists/ukrr-dialysis-opcs-4.csv",
    system="opcs4",
    column="code"
)
# kidney transplant
kidney_transplant_codes = codelist_from_csv(
  "codelists/opensafely-kidney-transplant.csv",
  system = "ctv3",
  column = "CTV3ID"
)
kidney_tx_icd10_codelist=codelist(["Z940"], system="icd10")
kidney_tx_opcs4_codelist = codelist_from_csv(
    "codelists/user-viyaasan-kidney-transplant-opcs-4.csv",
    system="opcs4",
    column="code"
)
# RRT
RRT_codelist = codelist_from_csv(
    "codelists/opensafely-renal-replacement-therapy.csv",
    system="ctv3",
    column="CTV3ID"
)
RRT_icd10_codelist = combine_codelists(
    kidney_tx_icd10_codelist,
    dialysis_icd10_codelist,
    codelist(["T861"], system="icd10")
)
RRT_opcs4_codelist = combine_codelists(
    kidney_tx_opcs4_codelist,
    dialysis_opcs4_codelist,
    codelist(["M023", "M026", "M027", "X412"], system="opcs4")
)
# blood creatinine
creatinine_codes_ctv3 = codelist(["XE2q5"], system="ctv3")
creatinine_codes_snomed = codelist_from_csv(
    "codelists/user-bangzheng-creatinine-value.csv", system="snomed", column="code"
)
creatinine_codes_short_snomed = codelist_from_csv(
    "codelists/user-bangzheng-creatinine-value-shortlist.csv", system="snomed", column="code"
)
# eGFR
eGFR_level_codelist = codelist_from_csv(
    "codelists/user-ss808-estimated-glomerular-filtration-rate-egfr-values.csv",
    system="snomed",
    column="code",
)
eGFR_short_level_codelist = codelist_from_csv(
    "codelists/user-bangzheng-egfr-value-shortlist.csv",
    system="snomed",
    column="code",
)
# Paxlovid interactions
drugs_do_not_use_codes = codelist_from_csv(
  "codelists/opensafely-sps-paxlovid-interactions-do-not-use-58f87823-dmd.csv", 
  system = "snomed", 
  column = "dmd_id"
)
drugs_consider_risk_codes = codelist_from_csv(
  "codelists/opensafely-nirmatrelvir-drug-interactions-3d3644f8-dmd.csv", 
  system = "snomed", 
  column = "dmd_id"
)

# Chronic cardiac disease
chronic_cardiac_dis_codes = codelist_from_csv(
    "codelists/opensafely-chronic-cardiac-disease-snomed.csv",
    system="snomed",
    column="id"
)
# Chronic respiratory disease
chronic_respiratory_dis_codes = codelist_from_csv(
    "codelists/opensafely-chronic-respiratory-disease-snomed.csv",
    system="snomed",
    column="id"
)
# Diabetes
diabetes_codes = codelist_from_csv(
    "codelists/opensafely-diabetes-snomed.csv",
    system="snomed",
    column="id"
)
# Hypertension
hypertension_codes = codelist_from_csv(
    "codelists/opensafely-hypertension-snomed.csv",
    system="snomed",
    column="id"
)

## ELIGIBILITY CRITERIA VARIABLES ----

### Onset of symptoms of COVID-19
covid_symptoms_snomed_codes = codelist_from_csv(
  "codelists/user-MillieGreen-covid-19-symptoms.csv",
  system = "snomed",
  column = "code",
)

### Require hospitalisation for COVID-19
covid_icd10_codes = codelist_from_csv(
  "codelists/opensafely-covid-identification.csv",
  system = "icd10",
  column = "icd10_code",
)

### Pregnancy
pregnancy_primis_codes = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-preg.csv",
  system = "snomed",
  column = "code",
)

### Pregnancy or delivery
pregdel_primis_codes = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-pregdel.csv",
  system = "snomed",
  column = "code",
)

### Weight
weight_opensafely_snomed_codes  = codelist_from_csv(
  "codelists/opensafely-weight-snomed.csv",
  system = "snomed",
  column = "code",
) 


## HIGH RISK GROUPS ----

### Down's syndrome
downs_syndrome_nhsd_snomed_codes = codelist_from_csv(
  "codelists/nhsd-downs-syndrome-snomed-ct.csv",
  system = "snomed",
  column = "code",
)

downs_syndrome_nhsd_icd10_codes = codelist_from_csv(
  "codelists/nhsd-downs-syndrome-icd-10.csv",
  system = "icd10",
  column = "code",
)

### Sickle cell disease
sickle_cell_disease_nhsd_snomed_codes = codelist_from_csv(
  "codelists/nhsd-sickle-spl-atriskv4-snomed-ct.csv",
  system = "snomed",
  column = "code",
)

sickle_cell_disease_nhsd_icd10_codes = codelist_from_csv(
  "codelists/nhsd-sickle-spl-hes-icd-10.csv",
  system = "icd10",
  column = "code",
)

### Solid cancer
non_haematological_cancer_opensafely_snomed_codes = codelist_from_csv(
  "codelists/opensafely-cancer-excluding-lung-and-haematological-snomed.csv",
  system = "snomed",
  column = "id",
)
non_haematological_cancer_opensafely_snomed_codes_new = codelist_from_csv(
  "codelists/user-bangzheng-cancer-excluding-lung-and-haematological-snomed-new.csv",
  system = "snomed",
  column = "code",
)
lung_cancer_opensafely_snomed_codes = codelist_from_csv(
  "codelists/opensafely-lung-cancer-snomed.csv", 
  system = "snomed", 
  column = "id"
)

chemotherapy_radiotherapy_opensafely_snomed_codes = codelist_from_csv(
  "codelists/opensafely-chemotherapy-or-radiotherapy-snomed.csv", 
  system = "snomed", 
  column = "id"
)

### Patients with a haematological diseases
haematopoietic_stem_cell_transplant_nhsd_snomed_codes = codelist_from_csv(
  "codelists/nhsd-haematopoietic-stem-cell-transplant-snomed.csv", 
  system = "snomed", 
  column = "code"
)

haematopoietic_stem_cell_transplant_nhsd_icd10_codes = codelist_from_csv(
  "codelists/nhsd-haematopoietic-stem-cell-transplant-icd-10.csv", 
  system = "icd10", 
  column = "code"
)

haematopoietic_stem_cell_transplant_nhsd_opcs4_codes = codelist_from_csv(
  "codelists/nhsd-haematopoietic-stem-cell-transplant-opcs4.csv", 
  system = "opcs4", 
  column = "code"
)

haematological_malignancies_nhsd_snomed_codes = codelist_from_csv(
  "codelists/nhsd-haematological-malignancies-snomed.csv",
  system = "snomed",
  column = "code"
)

haematological_malignancies_nhsd_icd10_codes = codelist_from_csv(
  "codelists/nhsd-haematological-malignancies-icd-10.csv", 
  system = "icd10", 
  column = "code"
)

### Patients with renal disease

#### CKD stage 5
ckd_stage_5_nhsd_snomed_codes = codelist_from_csv(
  "codelists/nhsd-ckd-stage-5-snomed-ct.csv", 
  system = "snomed", 
  column = "code"
)

ckd_stage_5_nhsd_icd10_codes = codelist_from_csv(
  "codelists/nhsd-ckd-stage-5-icd-10.csv", 
  system = "icd10", 
  column = "code"
)

### Patients with liver disease
liver_disease_nhsd_snomed_codes = codelist_from_csv(
  "codelists/nhsd-liver-cirrhosis.csv", 
  system = "snomed", 
  column = "code"
)

liver_disease_nhsd_icd10_codes = codelist_from_csv(
  "codelists/nhsd-liver-cirrhosis-icd-10.csv", 
  system = "icd10", 
  column = "code"
)

### Immune-mediated inflammatory disorders (IMID)
immunosuppresant_drugs_dmd_codes = codelist_from_csv(
  "codelists/nhsd-immunosuppresant-drugs-pra-dmd.csv", 
  system = "snomed", 
  column = "code"
)

immunosuppresant_drugs_snomed_codes = codelist_from_csv(
  "codelists/nhsd-immunosuppresant-drugs-pra-snomed.csv", 
  system = "snomed", 
  column = "code"
)

oral_steroid_drugs_dmd_codes = codelist_from_csv(
  "codelists/nhsd-oral-steroid-drugs-pra-dmd.csv",
  system = "snomed",
  column = "dmd_id",
)

oral_steroid_drugs_snomed_codes = codelist_from_csv(
  "codelists/nhsd-oral-steroid-drugs-snomed.csv", 
  system = "snomed", 
  column = "code"
)

### Primary immune deficiencies
immunosupression_nhsd_codes = codelist_from_csv(
  "codelists/nhsd-immunosupression-pcdcluster-snomed-ct.csv",
  system = "snomed",
  column = "code",
)
immunosupression_nhsd_codes_new = codelist_from_csv(
  "codelists/user-bangzheng-nhsd-immunosupression-pcdcluster-snomed-ct-new.csv",
  system = "snomed",
  column = "code",
)
## HIV/AIDs
hiv_aids_nhsd_snomed_codes = codelist_from_csv(
  "codelists/nhsd-hiv-aids-snomed.csv", 
  system = "snomed", 
  column = "code"
)

hiv_aids_nhsd_icd10_codes = codelist_from_csv(
  "codelists/nhsd-hiv-aids-icd10.csv", 
  system = "icd10", 
  column = "code"
)

## Solid organ transplant
solid_organ_transplant_codes = codelist_from_csv(
    "codelists/opensafely-solid-organ-transplantation-snomed.csv",
    system = "snomed",
    column = "id",
)

solid_organ_transplant_nhsd_snomed_codes = codelist_from_csv(
  "codelists/nhsd-transplant-spl-atriskv4-snomed-ct.csv",
  system = "snomed",
  column = "code",
)
solid_organ_transplant_nhsd_snomed_codes_new = codelist_from_csv(
  "codelists/user-bangzheng-nhsd-transplant-spl-atriskv4-snomed-ct-new.csv",
  system = "snomed",
  column = "code",
)
solid_organ_transplant_nhsd_opcs4_codes = codelist_from_csv(
  "codelists/nhsd-transplant-spl-hes-opcs4.csv", 
  system = "opcs4", 
  column = "code"
)

thymus_gland_transplant_nhsd_opcs4_codes = codelist_from_csv(
  "codelists/nhsd-transplant-thymus-gland-spl-hes-opcs4.csv", 
  system = "opcs4", 
  column = "code"
)

replacement_of_organ_transplant_nhsd_opcs4_codes = codelist_from_csv(
  "codelists/nhsd-transplant-replacement-of-organ-spl-hes-opcs4.csv", 
  system = "opcs4", 
  column = "code"
)

conjunctiva_transplant_nhsd_opcs4_codes = codelist_from_csv(
  "codelists/nhsd-transplant-conjunctiva-spl-hes-opcs4.csv", 
  system = "opcs4", 
  column = "code"
)

conjunctiva_y_codes_transplant_nhsd_opcs4_codes = codelist_from_csv(
  "codelists/nhsd-transplant-conjunctiva-y-codes-spl-hes-opcs4.csv", 
  system = "opcs4", 
  column = "code"
)

stomach_transplant_nhsd_opcs4_codes = codelist_from_csv(
  "codelists/nhsd-transplant-stomach-spl-hes-opcs4.csv", 
  system = "opcs4", 
  column = "code"
)

ileum_1_transplant_nhsd_opcs4_codes = codelist_from_csv(
  "codelists/nhsd-transplant-ileum_1-spl-hes-opcs4.csv", 
  system = "opcs4", 
  column = "code"
)

ileum_2_transplant_nhsd_opcs4_codes = codelist_from_csv(
  "codelists/nhsd-transplant-ileum_2-spl-hes-opcs4.csv", 
  system = "opcs4", 
  column = "code"
)

ileum_1_y_codes_transplant_nhsd_opcs4_codes = codelist_from_csv(
  "codelists/nhsd-transplant-ileum_1-y-codes-spl-hes-opcs4.csv", 
  system = "opcs4", 
  column = "code"
)

ileum_2_y_codes_transplant_nhsd_opcs4_codes = codelist_from_csv(
  "codelists/nhsd-transplant-ileum_2-y-codes-spl-hes-opcs4.csv", 
  system = "opcs4", 
  column = "code"
)

### Rare neurological conditions

#### Multiple sclerosis
multiple_sclerosis_nhsd_snomed_codes = codelist_from_csv(
  "codelists/nhsd-multiple-sclerosis-snomed-ct.csv",
  system = "snomed",
  column = "code",
)

multiple_sclerosis_nhsd_icd10_codes = codelist_from_csv(
  "codelists/nhsd-multiple-sclerosis.csv",
  system = "icd10",
  column = "code",
)

#### Motor neurone disease
motor_neurone_disease_nhsd_snomed_codes = codelist_from_csv(
  "codelists/nhsd-motor-neurone-disease-snomed-ct.csv",
  system = "snomed",
  column = "code",
)

motor_neurone_disease_nhsd_icd10_codes = codelist_from_csv(
  "codelists/nhsd-motor-neurone-disease-icd-10.csv",
  system = "icd10",
  column = "code",
)

#### Myasthenia gravis
myasthenia_gravis_nhsd_snomed_codes = codelist_from_csv(
  "codelists/nhsd-myasthenia-gravis-snomed-ct.csv",
  system = "snomed",
  column = "code",
)

myasthenia_gravis_nhsd_icd10_codes = codelist_from_csv(
  "codelists/nhsd-myasthenia-gravis.csv",
  system = "icd10",
  column = "code",
)

#### Huntington’s disease
huntingtons_disease_nhsd_snomed_codes = codelist_from_csv(
  "codelists/nhsd-huntingtons-snomed-ct.csv",
  system = "snomed",
  column = "code",
)

huntingtons_disease_nhsd_icd10_codes = codelist_from_csv(
  "codelists/nhsd-huntingtons.csv",
  system = "icd10",
  column = "code",
)  

## CLINICAL/DEMOGRAPHIC COVARIATES ----

### Ethnicity
ethnicity_primis_snomed_codes = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-eth2001.csv",
  system = "snomed",
  column = "code",
  category_column="grouping_6_id",
)


# OTHER COVARIATES ----
 
## Autism
autism_nhsd_snomed_codes = codelist_from_csv(
  "codelists/nhsd-primary-care-domain-refsets-autism_cod.csv",
  system = "snomed",
  column = "code",
)

## Care home 
care_home_primis_snomed_codes = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-longres.csv", 
    system = "snomed", 
    column = "code")

## Dementia
dementia_nhsd_snomed_codes = codelist_from_csv(
  "codelists/nhsd-primary-care-domain-refsets-dem_cod.csv", 
  system = "snomed", 
  column = "code",
)

## Housebound
housebound_opensafely_snomed_codes = codelist_from_csv(
    "codelists/opensafely-housebound.csv", 
    system = "snomed", 
    column = "code"
)

no_longer_housebound_opensafely_snomed_codes = codelist_from_csv(
    "codelists/opensafely-no-longer-housebound.csv", 
    system = "snomed", 
    column = "code"
)

## Learning disabilities
wider_ld_primis_snomed_codes = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-learndis.csv", 
    system = "snomed", 
    column = "code"
)

## Shielded
high_risk_primis_snomed_codes = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-shield.csv", 
    system = "snomed", 
    column = "code")

not_high_risk_primis_snomed_codes = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-nonshield.csv", 
    system = "snomed", 
    column = "code")
    
## Serious mental illness
serious_mental_illness_nhsd_snomed_codes = codelist_from_csv(
  "codelists/nhsd-primary-care-domain-refsets-mh_cod.csv",
  system = "snomed",
  column = "code",
)
    

## Vaccination declined
first_dose_declined = codelist_from_csv(
  "codelists/opensafely-covid-19-vaccination-first-dose-declined.csv",
  system = "snomed",
  column = "code",
)

second_dose_declined = codelist_from_csv(
  "codelists/opensafely-covid-19-vaccination-second-dose-declined.csv",
  system = "snomed",
  column = "code",
)

covid_vaccine_declined_codes = combine_codelists(
  first_dose_declined, second_dose_declined
)