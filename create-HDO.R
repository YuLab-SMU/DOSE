
# downloaded from https://www.alliancegenome.org/downloads
# access date: 2024-08-13

x <- read.delim(gzfile("HDO/DISEASE-ALLIANCE_HUMAN.tsv.gz"), comment.char="#")
x <- x[, c("DBObjectSymbol", "DOID")] |> setNames(c("SYMBOL", "id"))

library(clusterProfiler)
library(org.Hs.eg.db)

name2eg <- bitr(x[, 1], "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)
ont2gene <- merge(x, name2eg, by = 'SYMBOL')
ont2gene <- ont2gene[, -1] |> setNames(c("id", "gene"))

library(obolite)
date <- '20240731'

url <- "https://github.com/DiseaseOntology/HumanDiseaseOntology/blob/main/src/ontology/HumanDO.obo"

create_sqlite("HDO/HumanDO.obo", 
    "HDO.sqlite", 
    name = "Human Disease Ontology", 
    date, 
    url, 
    ont2gene = ont2gene)


