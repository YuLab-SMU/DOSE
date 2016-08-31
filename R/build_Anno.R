build_Anno <- function(path2gene, path2name) {
    if (!exists(".Anno_clusterProfiler_Env", envir = .GlobalEnv)) {
        pos <- 1
        envir <- as.environment(pos) 
        assign(".Anno_clusterProfiler_Env", new.env(), envir = envir)
    }
    Anno_clusterProfiler_Env <- get(".Anno_clusterProfiler_Env", envir= .GlobalEnv)

    path2gene <- as.data.frame(path2gene) 
    path2gene <- path2gene[!is.na(path2gene[,1]), ]
    path2gene <- path2gene[!is.na(path2gene[,2]), ]
    path2gene <- unique(path2gene)
    
    PATHID2EXTID <- split(as.character(path2gene[,2]), as.character(path2gene[,1]))
    EXTID2PATHID <- split(as.character(path2gene[,1]), as.character(path2gene[,2]))
    
    assign("PATHID2EXTID", PATHID2EXTID, envir = Anno_clusterProfiler_Env)
    assign("EXTID2PATHID", EXTID2PATHID, envir = Anno_clusterProfiler_Env)

    if ( missing(path2name) || is.null(path2name) || is.na(path2name)) {
        assign("PATHID2NAME", NULL, envir = Anno_clusterProfiler_Env)
    } else {
        path2name <- as.data.frame(path2name)
        path2name <- path2name[!is.na(path2name[,1]), ]
        path2name <- path2name[!is.na(path2name[,2]), ]
	path2name <- unique(path2name)
	PATH2NAME <- as.character(path2name[,2])
	names(PATH2NAME) <- as.character(path2name[,1]) 
        assign("PATHID2NAME", PATH2NAME, envir = Anno_clusterProfiler_Env)
    }
    return(Anno_clusterProfiler_Env)
}
