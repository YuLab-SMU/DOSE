## Network of Cancer Genes (NCG)
## repo: https://www.network-cancer-genes.org/
## download: https://www.network-cancer-genes.org/download.php
##
## The download endpoint returns a TSV via POST with columns:
##   entrez, symbol, pubmed_id, type, organ_system, primary_site,
##   cancer_type, method, coding_status, ...
##
## We extract cancer_type + entrez (gene ID) pairs and save them as GSON
## for enrichment analysis.

print("Processing NCG data...")

print("1. Downloading NCG data...")
url <- "http://www.network-cancer-genes.org/download.php"

response <- httr::POST(url, body = list(downloadcancergenes = "Download"))
resp_text <- httr::content(response, as = "text", encoding = "UTF-8")

if (nchar(resp_text) < 100) {
    stop("NCG download failed: response too short.")
}

print("2. Fetching NCG metadata...")
home_url <- "http://www.network-cancer-genes.org/"
home_response <- httr::GET(home_url)
home_text <- httr::content(home_response, as = "text", encoding = "UTF-8")
ncg_version <- regmatches(home_text, gregexpr("NCG[0-9]+(\\.[0-9]+)+", home_text))[[1]][1]
if (length(ncg_version) == 0 || is.na(ncg_version) || ncg_version == "") {
    stop("NCG version not found on homepage.")
}
ncg_accessed_date <- as.character(Sys.Date())
cat(sprintf("  Version: %s\n", ncg_version))
cat(sprintf("  Accessed date: %s\n", ncg_accessed_date))

print("3. Parsing TSV...")
x <- read.delim(text = resp_text, stringsAsFactors = FALSE, quote = "")
cat(sprintf("  Total rows: %d\n", nrow(x)))
cat(sprintf("  Columns: %s\n", paste(names(x), collapse = ", ")))

print("4. Extracting cancer_type + entrez...")
ncg <- x[, c("cancer_type", "entrez")]
ncg$cancer_type <- gsub('"', '', ncg$cancer_type, fixed = TRUE)
ncg <- ncg[ncg$cancer_type != "", ]
ncg <- ncg[!is.na(ncg$entrez), ]
ncg <- unique(ncg)

cat(sprintf("  After filtering: %d gene-disease pairs\n", nrow(ncg)))
cat(sprintf("  Unique cancer types: %d\n", length(unique(ncg$cancer_type))))
cat(sprintf("  Unique genes: %d\n", length(unique(ncg$entrez))))

print("5. Building NCG GSON...")
PATHID2EXTID <- split(as.character(ncg$entrez), as.character(ncg$cancer_type))

gsid2gene <- stack(PATHID2EXTID)
colnames(gsid2gene) <- c("gene", "gsid")
gsid2gene <- gsid2gene[, c("gsid", "gene")]

gsid2name <- data.frame(gsid = unique(ncg$cancer_type),
                        name = unique(ncg$cancer_type),
                        stringsAsFactors = FALSE)

ncg_gson <- gson::gson(gsid2gene = gsid2gene,
                       gsid2name = gsid2name,
                       species = "Homo sapiens",
                       gsname = "NCG",
                       keytype = "ENTREZID",
                       version = ncg_version,
                       accessed_date = ncg_accessed_date)

print("6. Saving NCG.gson.gz...")
outfile <- "NCG.gson"
gson::write.gson(ncg_gson, outfile)
R.utils::gzip(outfile, overwrite = TRUE)
outfile_gz <- paste0(outfile, ".gz")

cat(sprintf("\nDone. Output: %s\n", outfile_gz))
cat(sprintf("File size: %s\n", file.info(outfile_gz)$size))
