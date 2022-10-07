## do_rif.human.txt was downloaded from
## http://projects.bioinformatics.northwestern.edu/do_rif/
##
do.rif <- read.delim2(gzfile("do_rif.human.txt.gz"), sep="\t",
                      stringsAsFactors=F, header=F)
eg.do1 <- do.rif[,c(1,5)]
colnames(eg.do1) <- c("eg", "doid")

## IDMappings.txt from
## http://doa.nubic.northwestern.edu/pages/download.php
domapping <- read.delim(gzfile("IDMappings.txt.gz"), stringsAsFactors=F)
eg.do2 <- domapping[,c(2,1)]
colnames(eg.do2) <- c("eg", "doid")
eg.do2$doid <- paste("DOID:", eg.do2$doid, sep="")

eg.do <- rbind(eg.do2, eg.do1)
eg.do <- unique(eg.do)



## update in 2022_7_21
library(data.table)
anno <- fread("DISEASE-ALLIANCE_HUMAN.tsv.gz")[, c("DBObjectSymbol", "DOID")]
class(anno) <- "data.frame"
library(clusterProfiler)
library(org.Hs.eg.db)
anno_bitr <- bitr(anno[, 1], "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)
anno_bitr <- anno_bitr[!duplicated(anno_bitr[, 1]), ]
rownames(anno_bitr) <- anno_bitr[, 1]
anno[, 1] <- anno_bitr[anno[, 1], 2]
anno <- anno[!is.na(anno[, 1]), ]
colnames(anno) <- c("eg", "doid")
## add old data
load("olddata\\DO2EG.rda")
anno_old <- stack(DO2EG)
colnames(anno_old) <- c("eg", "doid")
anno <- rbind(anno, anno_old)
anno <- unique(anno)
eg.do <- anno
DOSE:::rebuildAnnoData.internal(eg.do)

# DOIC.rda 
# for new slot `geneAnno` of object `GOSemSimDATA`
DOSE:::build_dodata()
