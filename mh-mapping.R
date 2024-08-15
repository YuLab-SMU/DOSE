

# Mouse-Human Ontology Mapping Initiative (MHMI)
#
# https://github.com/mapping-commons/mh_mapping_initiative
#



dir <- 'mh_mapping_initiative/mappings'


read_mpo2hpo <- function(dir) {
    ff <- list.files(path=dir, pattern="^mp_hp", full.names=TRUE)

    res <- lapply(ff, function(f) {
            x <- read.delim(f, comment.char="#")
            x <- x[, c("subject_id", "object_id")]
            x <- x[x[,1] != "sssom:NoTermFound", ]
            x <- x[x[,2] != "sssom:NoTermFound", ]
            return(x)
        }) |> 
        rbindlist() |>
        unique() |>
        setNames(c("mpo_id", "hpo_id"))

    if (!all(grepl("^MP:", res[,1]))) {
        message("Not all ID in column 1 is started with MP, please double check...")
    }
    
    return(res)
}


mpo2hpo <- read_mpo2hpo(dir)

head(mpo2hpo)


## 这些映射的，后面再处理，然后比如说，小鼠的富集结果，就可以把这个映射到人的信息，用桑基图来呈现。

## 映射关系还有很多，后面再处理。如mpo2do, mpo2omim


