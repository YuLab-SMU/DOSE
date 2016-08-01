x <- read.delim("all_gene_disease_associations.tsv", comment.char="#", stringsAsFactor=F)
d2n <- unique(x[, c(4, 5)])
d2g <- unique(x[, c(4, 1)])

.DGN_DOSE_Env <- clusterProfiler:::build_Anno(d2g, d2n)


save(.DGN_DOSE_Env, file="DGN_DOSE_Env.rda", compress='xz')



y <- read.delim("all_snps_sentences_pubmeds.tsv", comment.char="#", stringsAsFactor=F)
d2n <- unique(y[, c(5, 6)])
d2s <- unique(y[, c(5, 1)])


.VDGN_DOSE_Env <- clusterProfiler:::build_Anno(d2s, d2n)

save(.VDGN_DOSE_Env, file="VDGN_DOSE_Env.rda", compress='xz')
