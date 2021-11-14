# download from http://ncg.kcl.ac.uk/download.php
# NCG 6.0: All cancer genes -> •List of 2372 cancer genes and supporting literature
# NCG 7.0: List of all 3347 cancer drivers and their annotation and supporting evidence
#x=read.delim("NCG6_cancergenes.tsv", stringsAsFactor=F)
x=read.delim("NCG6_cancergenes.tsv", stringsAsFactor=F, encoding = "latin1")
path2gene <- x[, c("cancer_type", "entrez")]
path2gene <- path2gene[path2gene[,1] != '',]

## gene2name <- x[, c("entrez", "symbol")]

path2name=NULL

.NCG_DOSE_Env <- DOSE:::build_Anno(path2gene, path2name)

NCG_EXTID2PATHID = get("EXTID2PATHID", envir=.NCG_DOSE_Env)
NCG_PATHID2EXTID = get("PATHID2EXTID", envir = .NCG_DOSE_Env)
NCG_PATHID2NAME = get("PATHID2NAME", envir = .NCG_DOSE_Env)

save(NCG_EXTID2PATHID, file = "NCG_EXTID2PATHID.rda", compress='xz')
save(NCG_PATHID2EXTID, file="NCG_PATHID2EXTID.rda", compress='xz')
save(NCG_PATHID2NAME, file="NCG_PATHID2NAME.rda", compress='xz')





