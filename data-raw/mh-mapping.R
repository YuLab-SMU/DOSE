

# Mouse-Human Ontology Mapping Initiative (MHMI)
#
# https://github.com/mapping-commons/mh_mapping_initiative
#


read.sssom <- function(file) {
    read.delim(file, comment.char = "#")
}



dir <- 'mh_mapping_initiative/mappings'


hpo2do <- read.sssom(file.path(dir, 'hp_doid_pistoia.sssom.tsv'))
    
hpo2do <- hpo2do[, c(1,3)] |>
    setNames(c("HPO", "DO")) |>
    unique()

dim(hpo2do)
head(hpo2do)    


hpo2omim <- read.sssom('phenotype.hpoa')[, c("hpo_id", "database_id")] |>
    setNames(c("hpo_id", "omim_id"))
do2omim <- read.sssom("OMIMinDO.tsv")[, c("id", "xrefs")] |>
    setNames(c("do_id", "omim_id"))
do2omim[,2] <- sprintf("O%s", do2omim[,2])
head(do2omim)

# too large, the relationship may not true
hpo2do2 <- merge(hpo2omim, do2omim, by='omim_id')[,c("hpo_id", "do_id")] |> unique()





read_mpo2hpo <- function(dir) {
    ff <- list.files(path=dir, pattern="^mp_hp", full.names=TRUE)

    res <- lapply(ff, function(f) {
            x <- read.sssom(f)
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


