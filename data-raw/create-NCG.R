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

## The NCG server (UCL) is frequently slow or briefly unreachable, which used
## to fail the whole weekly data refresh on a single connect timeout.
## Retry with a raised timeout and linear backoff: 30s, 60s, 90s.
fetch_ncg <- function(url, body = NULL, attempts = 4) {
    last_error <- NULL
    for (attempt in seq_len(attempts)) {
        res <- tryCatch({
            if (is.null(body)) {
                httr::GET(url, httr::timeout(120), httr::config(connecttimeout = 60))
            } else {
                httr::POST(url, body = body,
                           httr::timeout(120), httr::config(connecttimeout = 60))
            }
        }, error = function(e) {
            last_error <<- conditionMessage(e)
            NULL
        })

        txt <- NULL
        if (!is.null(res) && httr::status_code(res) < 400) {
            txt <- tryCatch(httr::content(res, as = "text", encoding = "UTF-8"),
                            error = function(e) NULL)
        }
        if (!is.null(txt) && !is.na(txt) && nchar(txt) >= 100) {
            return(txt)
        }

        if (attempt < attempts) {
            wait <- 30 * attempt
            message(sprintf("  attempt %d/%d failed (%s); retrying in %ds ...",
                            attempt, attempts,
                            if (is.null(last_error)) "unusable response" else last_error,
                            wait))
            Sys.sleep(wait)
        }
    }
    NULL
}

print("Processing NCG data...")

print("1. Downloading NCG data...")
url <- "http://www.network-cancer-genes.org/download.php"

resp_text <- fetch_ncg(url, body = list(downloadcancergenes = "Download"))

if (is.null(resp_text)) {
    stop("NCG download failed: ", url, " unreachable after 4 attempts. ",
         "The previously published NCG.gson.gz is left untouched.")
}

print("2. Fetching NCG metadata...")
home_url <- "http://www.network-cancer-genes.org/"
home_text <- fetch_ncg(home_url, attempts = 2)
if (is.null(home_text)) {
    stop("NCG metadata fetch failed: homepage unreachable.")
}
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
