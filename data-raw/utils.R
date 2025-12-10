
library(jsonlite)

get_release <- function(repo, files, dir = ".") {
    # 1. 获取最新 release 信息
    api_url <- paste0("https://api.github.com/repos/", repo, "/releases/latest")
    release <- jsonlite::fromJSON(api_url)

    # 2. 提取发布日期 (格式 YYYYMMDD)
    date <- format(as.Date(release$published_at), "%Y%m%d")

    # 3. 找到目标文件的下载链接
    assets <- release$assets
    
    # 4. 下载文件
    if (!dir.exists(dir)) dir.create(dir)

    dl_files <- list()
    for (f in files) {
        url <- assets$browser_download_url[assets$name == f]
        if (length(url) > 0) {
            message("Downloading ", f, "...")
            download.file(url, file.path(dir, f), mode = "wb")
            dl_files[[f]] <- url
        } else {
            warning(f, " not found in release assets.")
        }
    }

    # 5. 返回结果
    return(list(
        files = dl_files,
        date = date
    ))
}


download_file <- function(url, outfile) {
    if (!file.exists(dirname(outfile))) {
        dir.create(dirname(outfile))
    }
    
    message("Downloading ", outfile, "...")
    download.file(url, outfile, mode = "wb")
}
