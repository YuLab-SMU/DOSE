get_remote_file_date <- function(owner, repo, path) {
  url <- paste0("https://api.github.com/repos/", owner, "/", repo, 
                "/commits?path=", path, "&per_page=1")
  d <- jsonlite::fromJSON(url)
  # 转换为 YYYYMMDD 格式
  format(as.Date(d$commit$committer$date), "%Y%m%d") 
}
