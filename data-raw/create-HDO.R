
# downloaded from https://www.alliancegenome.org/downloads
# access date: 2024-08-13

print("Processing HDO data...")

if (!dir.exists("HDO")) dir.create("HDO")

library(wget)
wget_set()

print("1. Downloading HDO files...")
download.file("https://fms.alliancegenome.org/download/DISEASE-ALLIANCE_HUMAN.tsv.gz", 
    "HDO/DISEASE-ALLIANCE_HUMAN.tsv.gz") # , method = "curl", extra = "-L")

print("2. Parsing DISEASE-ALLIANCE_HUMAN.tsv.gz file...")
x <- read.delim(gzfile("HDO/DISEASE-ALLIANCE_HUMAN.tsv.gz"), comment.char="#")
x <- x[, c("DBObjectSymbol", "DOID")] |> setNames(c("SYMBOL", "id"))

library(clusterProfiler)
library(org.Hs.eg.db)

name2eg <- bitr(x[, 1], "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)
ont2gene <- merge(x, name2eg, by = 'SYMBOL')
ont2gene <- ont2gene[, -1] |> setNames(c("id", "gene"))

print("3. Parsing OBO file...")
source("loadobolite.r")

url <- "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/HumanDO.obo"
download.file(url, "HDO/HumanDO.obo")

if (file.info("HDO/HumanDO.obo")$size < 100 * 1024) {
    stop("HDO/HumanDO.obo file too small. Download failed.")
}

source("get-remote-file-date.r")
date <- get_remote_file_date("DiseaseOntology", "HumanDiseaseOntology", "src/ontology/HumanDO.obo")

create_sqlite("HDO/HumanDO.obo",
              "HDO.sqlite",
              name = "Human Disease Ontology",
              date,
              url,
              ont2gene = ont2gene)

print("Finish HDO data creation.")
cat("\n\n")
