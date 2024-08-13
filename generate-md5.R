ff <- list.files(pattern = ".sqlite$")
md5 <- vapply(ff, function(f) digest::digest(f, algo='md5', file=TRUE), character(1))
cat(sprintf("%s\t%s\n", ff, md5), file="md5.txt", sep="")

