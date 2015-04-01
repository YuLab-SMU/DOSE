x=read.delim("annotationData/DOSE/cancergenes_4.0_20150331.txt", stringsAsFactor=F)
path2gene <- x[, c("cancer_type", "entrez")]
path2gene <- path2gene[path2gene[,1] != '',]

## gene2name <- x[, c("entrez", "symbol")]


NCG_DOSE_Env <- clusterProfiler:::buildAnno(path2gene, path2name)


save(NCG_DOSE_Env, file="NCG_DOSE_ENV.rda")



