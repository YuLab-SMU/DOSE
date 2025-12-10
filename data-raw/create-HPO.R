## Human Phenotype Ontology
## repo: https://github.com/obophenotype/human-phenotype-ontology
## release: https://github.com/obophenotype/human-phenotype-ontology/releases

source("utils.R")

print("Processing HPO data...")

print("1. Downloading HPO files...")
res <- get_release(repo = "obophenotype/human-phenotype-ontology",
            files = c("phenotype_to_genes.txt", "hp.obo"),
            dir = "HPO"
        )

print("2. Parsing Phenotype-to-Gene file...")
pg <- read.delim("HPO/phenotype_to_genes.txt")
hpo2gene <- pg[, c("hpo_id", "ncbi_gene_id")]
hpo2gene <- na.omit(unique(hpo2gene)) |> setNames(c("id", "gene"))

# updated date
date <- res$date


print("3. Parsing OBO file...")
source("loadobolite.r")

create_sqlite("HPO/hp.obo", "HPO.sqlite", 
    name = "Human Phenotype Ontology",
    date = date,
    url = res$files[["hp.obo"]], 
    ont2gene = hpo2gene
    )



print("Finish HPO data creation.")
cat("\n\n")
