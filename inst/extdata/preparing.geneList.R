require(breastCancerMAINZ)
data(mainz)

require("hgu133a.db")

require(siggenes)

clmainz=pData(mainz)$grade

dd <- exprs(mainz)
g1 <- dd[,clmainz == 1]
g3 <- dd[,clmainz == 3]
geneList <- exp(rowMeans(g3))/exp(rowMeans(g1))
geneList <- sort(geneList, decreasing=TRUE)
geneList <- log(geneList, base=2)

eg <- mget(names(geneList), hgu133aENTREZID, ifnotfound=NA)

gg <- data.frame(probe=names(geneList), val = geneList)
eg.df <- data.frame(probe=names(eg), eg=unlist(eg))

xx <- merge(gg, eg.df, by.x="probe", by.y="probe")
xx <- xx[,-1]
xx <- unique(xx)
xx <- xx[!is.na(xx[,2]),]

require(plyr)
yy <- ddply(xx, .(eg), function(x) data.frame(val=mean(x$val)))

geneList <- yy$val
names(geneList) <- yy$eg
geneList <- sort(geneList, decreasing=TRUE)

save(geneList, file="geneList.rda")
