## Human Phenotype Ontology
## repo: https://github.com/obophenotype/human-phenotype-ontology
## release: https://github.com/obophenotype/human-phenotype-ontology/releases


source("utils.R")

res <- get_release(repo = "obophenotype/human-phenotype-ontology",
            files = c("phenotype_to_genes.txt", "hp.obo"),
            dir = "HPO"
        )

pg <- read.delim("HPO/phenotype_to_genes.txt")
hpo2gene <- pg[, c("hpo_id", "ncbi_gene_id")]
hpo2gene <- na.omit(unique(hpo2gene)) |> setNames(c("id", "gene"))

# updated date
date <- res$date


source("loadobolite.r")

create_sqlite("HPO/hp.obo", "HPO.sqlite", 
    name = "Human Phenotype Ontology",
    date = date,
    url = res$files[["hp.obo"]], 
    ont2gene = hpo2gene
    )


