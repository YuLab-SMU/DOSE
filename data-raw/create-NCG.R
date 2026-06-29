## Network of Cancer Genes (NCG)
## repo: https://www.network-cancer-genes.org/
## download: https://www.network-cancer-genes.org/download.php
##
## The download endpoint returns a TSV via POST with columns:
##   entrez, symbol, pubmed_id, type, organ_system, primary_site,
##   cancer_type, method, coding_status, ...
##
## We extract cancer_type + entrez (gene ID) pairs for enrichment analysis.

print("Processing NCG data...")

print("1. Downloading NCG data...")
url <- "http://www.network-cancer-genes.org/download.php"

response <- httr::POST(url, body = list(downloadcancergenes = "Download"))
resp_text <- httr::content(response, as = "text", encoding = "UTF-8")

if (nchar(resp_text) < 100) {
    stop("NCG download failed: response too short.")
}

print("2. Parsing TSV...")
x <- read.delim(text = resp_text, stringsAsFactors = FALSE, quote = "")
cat(sprintf("  Total rows: %d\n", nrow(x)))
cat(sprintf("  Columns: %s\n", paste(names(x), collapse = ", ")))

print("3. Extracting cancer_type + entrez...")
ncg <- x[, c("cancer_type", "entrez")]
ncg$cancer_type <- gsub('"', '', ncg$cancer_type, fixed = TRUE)
ncg <- ncg[ncg$cancer_type != "", ]
ncg <- ncg[!is.na(ncg$entrez), ]
ncg <- unique(ncg)

cat(sprintf("  After filtering: %d gene-disease pairs\n", nrow(ncg)))
cat(sprintf("  Unique cancer types: %d\n", length(unique(ncg$cancer_type))))
cat(sprintf("  Unique genes: %d\n", length(unique(ncg$entrez))))

print("4. Saving NCG.tsv...")
outfile <- "NCG.tsv"
write.table(ncg, outfile, sep = "\t", row.names = FALSE)

cat(sprintf("\nDone. Output: %s\n", outfile))
cat(sprintf("File size: %s\n", file.info(outfile)$size))
