## Mouse Phenotype Ontology
## repo: https://github.com/mgijax/mammalian-phenotype-ontology
## release: https://github.com/mgijax/mammalian-phenotype-ontology/releases

# Mammalian Phenotype (MP)-Mouse Developmental Anatomy (EMAPA) Mappings
# file: MP_EMAPA.rpt
# downloaded at: https://www.informatics.jax.org/downloads/reports/index.html

print("Processing MPO data...")

source("utils.R")

print("1. Downloading MPO files...")
download_file("https://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt", "MPO/HMD_HumanPhenotype.rpt")
download_file("https://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt", "MPO/MGI_Gene_Model_Coord.rpt")

res <- get_release(repo = "mgijax/mammalian-phenotype-ontology",
            files = c("mp.obo"),
            dir = "MPO"
        )

read.rpt <- function(file, header=FALSE, ...) {
    read.delim(file, header=header, ...)
}

print("2. Parsing HMD_HumanPhenotype.rpt file...")
pg <- read.rpt("MPO/HMD_HumanPhenotype.rpt")
gene2mpo <- pg[, 4:5]
gene2mpo <- gene2mpo[gene2mpo[,2] != "",]
mpo2mgi <- strsplit(gene2mpo[,2], split=", ") |>
    setNames(gene2mpo[,1]) |>
    stack() |>
    setNames(c("id", "mgi")) |>
    unique()

print("3. Parsing MGI_Gene_Model_Coord.rpt file...")
x <- readr::read_tsv("MPO/MGI_Gene_Model_Coord.rpt")
mgi2eg <- x[,c("1. MGI accession id", "6. Entrez gene id")]
names(mgi2eg) <- c("mgi", "gene")
mpo2gene <- merge(mpo2mgi, mgi2eg, by='mgi')
mpo2gene <- unique(mpo2gene[, -1])

# updated date
date <- res$date

print("4. Parsing OBO file...")
source("loadobolite.r")

create_sqlite("MPO/mp.obo", "MPO.sqlite", 
    name = "Mouse Phenotype Ontology",
    date = date,
    url = res$files[["mp.obo"]], 
    ont2gene = mpo2gene
    )


print("Finish MPO data creation.")
