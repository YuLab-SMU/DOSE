
parse_do <- function(obofile) {
  x = readLines(obofile)
  
  start = grep("^\\[Term\\]", x)
  end <- c(start[-1] -1, length(x))

    res <- lapply(seq_along(start), function(i) {
    extract_do_item(x[start[i]:end[i]])
  })
  
  j <- vapply(res, is.null, logical(1))
  message(paste0(sum(j), '/', length(j)), " obsolete terms found.") 
  res <- res[!j]

  doinfo <- lapply(res, function(x) x$do) %>% do.call('rbind', .) %>% as.data.frame()

  rel <- lapply(res, function(x) x$relation) %>% do.call('rbind', .) 

  return(list(doinfo = doinfo, rel = rel))  
}


extract_do_item <- function(item) {
  i <- grep('^\\[Typedef\\]', item)
  if (length(i) > 0) {
    item <- item[-(i[1]:length(item))]
  }
  ## is_obsolete: true
  useless <- get_do_info(item, '^is_obsolete:')
  if (!is.na(useless)) return(NULL)
  
  id <- get_do_info(item, "^id:")
  name <- get_do_info(item, "^name:")
  def <- get_do_info(item, "^def:") 
  def <- sub('\\"', "", def)
  def <- sub('\\".*', "", def)
  
  isa <- get_do_info(item, '^is_a:')
  isa <- sub("\\s*!.*", "", isa)
  res <- list(do=c(id=id, name=name, def=def),
              relationship = data.frame(id=id, parent=isa))
}

get_do_info <- function(item, pattern) {
  i <- grep(pattern, item)
  if (length(i) == 0) return(NA)
  
  sub("\\s*", "", 
      sub(pattern, "", item[i])
  )
}


#x <- parse_do('HumanDO.obo')

