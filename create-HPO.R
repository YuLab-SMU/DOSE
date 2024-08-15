## Human Phenotype Ontology
## repo: https://github.com/obophenotype/human-phenotype-ontology
## release: https://github.com/obophenotype/human-phenotype-ontology/releases


pg <- read.delim("HPO/phenotype_to_genes.txt")
hpo2gene <- pg[, c("hpo_id", "ncbi_gene_id")]
hpo2gene <- na.omit(unique(hpo2gene)) |> setNames(c("id", "gene"))

# updated date
date <- 20240813


library(obolite)
create_sqlite("HPO/hp.obo", "HPO.sqlite", 
    name = "Human Phenotype Ontology",
    date = date,
    url = 'https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2024-08-13/hp.obo', 
    ont2gene = hpo2gene
    )


