### create sqlite

library(obolite)
date <- '20240628'
name <- "Disease Ontology"
url <- "https://github.com/DiseaseOntology/HumanDiseaseOntology/blob/main/src/ontology/HumanDO.obo"

create_sqlite("Downloads/HumanDO.obo", "HDO.sqlite", name, date, url)

