do.rif <- read.delim2(gzfile("do_rif.human.txt.gz"), sep="\t",
                      stringsAsFactors=F, header=F)
eg.do1 <- do.rif[,c(1,5)]
colnames(eg.do1) <- c("eg", "doid")

domapping <- read.delim(gzfile("IDMappings.txt.gz"), stringsAsFactors=F)
eg.do2 <- domapping[,c(2,1)]
colnames(eg.do2) <- c("eg", "doid")
eg.do2$doid <- paste("DOID:", eg.do2$doid, sep="")

eg.do <- rbind(eg.do2, eg.do1)
eg.do <- unique(eg.do)


DOSE:::rebuildAnnoData.internal(eg.do)
