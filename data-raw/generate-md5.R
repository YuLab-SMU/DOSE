ff <- list.files(pattern = ".sqlite$")
ff2 <- list.files(pattern = "\\.(tsv|gson)\\.gz$")
md5 <- vapply(ff, function(f) digest::digest(f, algo='md5', file=TRUE), character(1))
md5_tsv <- vapply(ff2, function(f) digest::digest(f, algo='md5', file=TRUE), character(1))
md5 <- c(md5, md5_tsv)

if (file.exists("md5.txt")) {
    x <- read.table("md5.txt")
    oldmd5 <- setNames(x[,2], x[,1])
    md5 <- c(md5, oldmd5[!names(oldmd5) %in% names(md5)])
}


cat(sprintf("%s\t%s\n", names(md5), md5), file="md5.txt", sep="")

sapply(ff, R.utils::gzip, overwrite = TRUE)

