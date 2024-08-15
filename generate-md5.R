ff <- list.files(pattern = ".sqlite$")
md5 <- vapply(ff, function(f) digest::digest(f, algo='md5', file=TRUE), character(1))

if (file.exists("md5.txt")) {
    x <- read.delim("md5.txt", header=F)
    oldmd5 <- setNames(x[,2], x[,1])
    md5 <- c(md5, oldmd5[!names(oldmd5) %in% names(md5)])
}


cat(sprintf("%s\t%s\n", ff, md5), file="md5.txt", sep="")

sapply(ff, R.utils::gzip, overwrite = TRUE)

