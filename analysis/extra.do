EXTRAS

//Corrections 
//1. gen imid_on_drug_clean = 1 if imid_on_drug==1 & (imid_drug==1 | imid_drug_hcd==1)
//2. gen eligible_clean=1 if downs_syndrome_nhsd==1 |cancer_opensafely_snomed==1 |haematological_disease_nhsd==1 |ckd_stage_5_nhsd==1 |liver_disease_nhsd==1 | immunosupression_nhsd==1 |imid_on_drug_clean==1| hiv_aids_nhsd==1 |solid_organ_transplant_nhsd==1|neurological_disease_nhsd==1 