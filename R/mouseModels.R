#' Rank mouse disease models from human query genes
#'
#' @param gene A character vector of human Entrez gene IDs.
#' @param disease_model_file Path to an `MGI_DiseaseMouseModel.rpt`-format file.
#' @param homology_file Path to a `HOM_MouseHumanSequence.rpt`-format file.
#' @param gene_pheno_file Path to an `MGI_GenePheno.rpt`-format file.
#' @param top Maximum number of mouse-model targets to return.
#'
#' @return A `doseInterpretResult` with `target_type = "mouse_model"`.
#' @export
rankMouseModels <- function(gene,
                            disease_model_file,
                            homology_file,
                            gene_pheno_file,
                            top = 50) {
    .validate_entrez_ids(gene)

    if (length(top) != 1L || is.na(top) || !is.numeric(top) || top < 1) {
        stop("`top` must be one positive number.", call. = FALSE)
    }

    disease_models <- .read_mgi_disease_mouse_model(disease_model_file)
    homology <- .build_mouse_human_ortholog_map(.read_hom_mouse_human_sequence(homology_file))
    gene_pheno <- .read_mgi_gene_pheno(gene_pheno_file)

    disease_models <- disease_models[is.na(disease_models$not_model) | !nzchar(disease_models$not_model), , drop = FALSE]
    query_map <- homology[homology$human_entrez_gene_id %in% as.character(gene), , drop = FALSE]

    if (!nrow(query_map)) {
        return(.empty_mouse_model_result(
            gene = gene,
            top = top
        ))
    }

    candidate_models <- disease_models[
        disease_models$marker_mgi_id %in% query_map$mouse_mgi_id,
        ,
        drop = FALSE
    ]

    if (!nrow(candidate_models)) {
        return(.empty_mouse_model_result(
            gene = gene,
            top = top
        ))
    }

    candidate_models$.contract_target_key <- .mouse_model_target_key(candidate_models)
    split_candidates <- split(candidate_models, candidate_models$.contract_target_key)
    model_parts <- lapply(split_candidates, function(rows) {
        support_map <- query_map[query_map$mouse_mgi_id %in% rows$marker_mgi_id, , drop = FALSE]
        .mouse_model_group_to_contract(rows, support_map, gene_pheno)
    })

    model_parts <- Filter(function(x) !is.null(x$result), model_parts)
    if (!length(model_parts)) {
        return(.empty_mouse_model_result(
            gene = gene,
            top = top
        ))
    }

    result_df <- do.call(rbind, lapply(model_parts, function(part) part$result))
    evidence_df <- do.call(rbind, lapply(model_parts, function(part) part$evidence))
    rownames(result_df) <- NULL
    rownames(evidence_df) <- NULL

    ord <- order(-result_df$score, result_df$target_name, result_df$target_id)
    result_df <- result_df[ord, , drop = FALSE]
    result_df <- result_df[seq_len(min(nrow(result_df), top)), , drop = FALSE]
    result_df$rank <- seq_len(nrow(result_df))

    keep_ids <- result_df$target_id
    evidence_df <- evidence_df[evidence_df$target_id %in% keep_ids, , drop = FALSE]
    rownames(evidence_df) <- NULL

    doseInterpretResult(
        result = result_df,
        evidence = evidence_df,
        query = .build_query_df(
            input = "gene",
            organism = "human",
            query_id = "q1"
        ),
        sources = data.frame(
            source = c("MGI_DiseaseMouseModel", "HOM_MouseHumanSequence", "MGI_GenePheno"),
            version = rep("supplied_file", 3),
            checksum = rep(NA_character_, 3),
            stringsAsFactors = FALSE
        ),
        parameters = list(
            input = "gene",
            organism = "human",
            target = "mouse_model",
            top = top,
            disease_model_file = normalizePath(disease_model_file, mustWork = TRUE),
            homology_file = normalizePath(homology_file, mustWork = TRUE),
            gene_pheno_file = normalizePath(gene_pheno_file, mustWork = TRUE)
        ),
        explanation = .empty_explanation()
    )
}

.normalize_contract_names <- function(x) {
    colnames(x) <- gsub("^_+|_+$", "", gsub("[^a-z0-9]+", "_", tolower(colnames(x))))
    x
}

.read_mgi_disease_mouse_model <- function(file) {
    if (!file.exists(file)) {
        stop("`disease_model_file` does not exist.", call. = FALSE)
    }

    x <- utils::read.delim(
        file,
        comment.char = "#",
        header = FALSE,
        fill = TRUE,
        stringsAsFactors = FALSE
    )
    colnames(x) <- c(
        "disease_name",
        "disease_id",
        "not_model",
        "allele_pairs",
        "strain_background",
        "allele_symbol",
        "allele_mgi_id",
        "allele_reference_count",
        "allele_repository_id",
        "allele_rrid",
        "marker_symbol",
        "marker_mgi_id",
        "gene_repository_id"
    )
    x
}

.read_hom_mouse_human_sequence <- function(file) {
    if (!file.exists(file)) {
        stop("`homology_file` does not exist.", call. = FALSE)
    }

    x <- utils::read.delim(
        file,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    .normalize_contract_names(x)
}

.read_mgi_gene_pheno <- function(file) {
    if (!file.exists(file)) {
        stop("`gene_pheno_file` does not exist.", call. = FALSE)
    }

    x <- utils::read.delim(
        file,
        header = FALSE,
        fill = TRUE,
        stringsAsFactors = FALSE
    )
    colnames(x) <- c(
        "genotype_symbol",
        "allele_symbol",
        "allele_mgi_id",
        "strain_background",
        "mp_id",
        "reference_id",
        "marker_mgi_id",
        "genotype_mgi_id"
    )
    x
}

.build_mouse_human_ortholog_map <- function(homology) {
    mouse <- homology[
        homology$common_organism_name == "mouse, laboratory",
        c("db_class_key", "symbol", "entrezgene_id", "mouse_mgi_id"),
        drop = FALSE
    ]
    human <- homology[
        homology$common_organism_name == "human",
        c("db_class_key", "symbol", "entrezgene_id"),
        drop = FALSE
    ]

    out <- merge(mouse, human, by = "db_class_key", suffixes = c("_mouse", "_human"))
    colnames(out) <- c(
        "db_class_key",
        "mouse_symbol",
        "mouse_entrez_gene_id",
        "mouse_mgi_id",
        "human_symbol",
        "human_entrez_gene_id"
    )
    out
}

.mouse_model_target_key <- function(x) {
    genotype_key <- trimws(as.character(x$allele_pairs))
    genotype_key[is.na(genotype_key)] <- ""

    allele_id_key <- trimws(as.character(x$allele_mgi_id))
    allele_id_key[is.na(allele_id_key)] <- ""

    allele_symbol_key <- trimws(as.character(x$allele_symbol))
    allele_symbol_key[is.na(allele_symbol_key)] <- ""

    marker_key <- trimws(as.character(x$marker_mgi_id))
    marker_key[is.na(marker_key)] <- ""

    fallback <- ifelse(
        nzchar(allele_id_key),
        allele_id_key,
        ifelse(nzchar(allele_symbol_key), allele_symbol_key, marker_key)
    )

    ifelse(
        nzchar(genotype_key),
        paste0("genotype::", genotype_key),
        paste0("allele::", fallback)
    )
}

.slug_mouse_model <- function(x) {
    gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", x))
}

.empty_mouse_model_result <- function(gene, top) {
    doseInterpretResult(
        result = .empty_result_df(),
        evidence = .empty_evidence_df(),
        query = .build_query_df(
            input = "gene",
            organism = "human",
            query_id = "q1"
        ),
        sources = data.frame(
            source = c("MGI_DiseaseMouseModel", "HOM_MouseHumanSequence", "MGI_GenePheno"),
            version = rep("supplied_file", 3),
            checksum = rep(NA_character_, 3),
            stringsAsFactors = FALSE
        ),
        parameters = list(
            input = "gene",
            organism = "human",
            target = "mouse_model",
            top = top,
            query_gene = as.character(gene)
        ),
        explanation = .empty_explanation()
    )
}

.mouse_model_group_to_contract <- function(rows, support_map, gene_pheno) {
    genotype_hits <- gene_pheno[gene_pheno$genotype_symbol %in% rows$allele_pairs, , drop = FALSE]
    allele_hits <- gene_pheno[gene_pheno$allele_mgi_id %in% rows$allele_mgi_id, , drop = FALSE]

    if (nrow(genotype_hits)) {
        phenotype_rows <- genotype_hits
        target_id <- unique(genotype_hits$genotype_mgi_id)[1]
        target_name <- unique(as.character(genotype_hits$genotype_symbol))
        target_name <- target_name[!is.na(target_name) & nzchar(target_name)]
        target_name <- if (length(target_name)) target_name[[1]] else rows$allele_pairs[[1]]
    } else {
        phenotype_rows <- allele_hits
        target_id <- sprintf(
            "ALLELE_MODEL:%s",
            .slug_mouse_model(rows$allele_mgi_id[[1]])
        )
        target_name <- rows$allele_symbol[[1]]
    }

    if (nrow(phenotype_rows)) {
        phenotype_rows <- unique(phenotype_rows[, c(
            "genotype_symbol",
            "allele_symbol",
            "allele_mgi_id",
            "strain_background",
            "mp_id",
            "reference_id",
            "marker_mgi_id",
            "genotype_mgi_id"
        ), drop = FALSE])
        rownames(phenotype_rows) <- NULL
    }

    ortholog_rows <- unique(support_map[, c(
        "db_class_key",
        "human_entrez_gene_id",
        "human_symbol",
        "mouse_symbol",
        "mouse_mgi_id"
    )])

    disease_labels <- unique(sprintf("%s (%s)", rows$disease_name, rows$disease_id))
    ortholog_count <- nrow(ortholog_rows)

    result_df <- data.frame(
        query_id = "q1",
        target_id = target_id,
        target_name = target_name,
        target_type = "mouse_model",
        rank = 1L,
        score = ortholog_count,
        score_type = "ortholog_support_count",
        score_direction = "higher",
        pvalue = NA_real_,
        p.adjust = NA_real_,
        evidence_count = ortholog_count + max(1L, nrow(phenotype_rows)),
        top_evidence_type = if (ortholog_count) "ortholog" else "phenotype_profile",
        source_count = 3L,
        stringsAsFactors = FALSE
    )

    ortholog_evidence <- data.frame(
        evidence_id = sprintf(
            "%s::ortholog::%s",
            target_id,
            ortholog_rows$human_entrez_gene_id
        ),
        query_id = "q1",
        target_id = target_id,
        target_type = "mouse_model",
        evidence_type = "ortholog",
        direction = "support",
        feature_id = ortholog_rows$human_entrez_gene_id,
        feature_name = ortholog_rows$human_symbol,
        component_score = rep(1, nrow(ortholog_rows)),
        source = "HOM_MouseHumanSequence",
        source_record_id = as.character(ortholog_rows$db_class_key),
        evidence_path_id = sprintf("%s::ortholog", target_id),
        derived_from_evidence_id = rep(NA_character_, nrow(ortholog_rows)),
        reference_id = rep(NA_character_, nrow(ortholog_rows)),
        note = rep(
            sprintf(
                "Supports mouse marker %s through the supplied mouse-human ortholog class.",
                ortholog_rows$mouse_symbol[[1]]
            ),
            nrow(ortholog_rows)
        ),
        stringsAsFactors = FALSE
    )

    phenotype_evidence <- if (nrow(phenotype_rows)) {
        data.frame(
            evidence_id = sprintf(
                "%s::phenotype_profile::%s::%s",
                target_id,
                phenotype_rows$mp_id,
                ifelse(
                    is.na(phenotype_rows$reference_id) | !nzchar(phenotype_rows$reference_id),
                    "na",
                    phenotype_rows$reference_id
                )
            ),
            query_id = "q1",
            target_id = target_id,
            target_type = "mouse_model",
            evidence_type = "phenotype_profile",
            direction = "support",
            feature_id = phenotype_rows$mp_id,
            feature_name = phenotype_rows$mp_id,
            component_score = rep(NA_real_, nrow(phenotype_rows)),
            source = "MGI_GenePheno",
            source_record_id = phenotype_rows$genotype_mgi_id,
            evidence_path_id = sprintf("%s::phenotype_profile", target_id),
            derived_from_evidence_id = rep(NA_character_, nrow(phenotype_rows)),
            reference_id = phenotype_rows$reference_id,
            note = rep(
                sprintf(
                    "Model phenotype support for %s.",
                    paste(disease_labels, collapse = "; ")
                ),
                nrow(phenotype_rows)
            ),
            stringsAsFactors = FALSE
        )
    } else {
        data.frame(
            evidence_id = sprintf("%s::phenotype_profile::missing", target_id),
            query_id = "q1",
            target_id = target_id,
            target_type = "mouse_model",
            evidence_type = "phenotype_profile",
            direction = "missing",
            feature_id = NA_character_,
            feature_name = NA_character_,
            component_score = NA_real_,
            source = "MGI_GenePheno",
            source_record_id = target_id,
            evidence_path_id = sprintf("%s::phenotype_profile", target_id),
            derived_from_evidence_id = NA_character_,
            reference_id = NA_character_,
            note = sprintf(
                "No matching model-side MP phenotype profile was found for %s.",
                target_name
            ),
            stringsAsFactors = FALSE
        )
    }

    list(
        result = result_df,
        evidence = rbind(ortholog_evidence, phenotype_evidence)
    )
}
