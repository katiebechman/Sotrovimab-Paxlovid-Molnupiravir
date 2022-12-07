from cohortextractor import StudyDefinition, patients, codelist, codelist_from_csv, combine_codelists, filter_codes_by_category

## Define study population and variables
from codelists import *
start_date = "2021-12-16"
end_date = date.today().isoformat()

study = StudyDefinition(
    ## Configure the expectations framework
    default_expectations={
        "date": {"earliest": "2021-11-01", "latest": "today"},
        "rate": "uniform",
        "incidence": 0.05,
    },

    ## Define index date
    index_date = start_date,
    
    ## Population
    population=patients.satisfying(
        """
        age >=18 AND age <= 110
        AND NOT has_died
        AND registered_treated 
        AND (sotrovimab_covid_therapeutics OR paxlovid_covid_therapeutics OR molnupiravir_covid_therapeutics)
        """,
  ),

  registered_treated = patients.registered_as_of("date_treated"), 
