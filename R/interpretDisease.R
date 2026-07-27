#' Canonical interpretation result class
#'
#' `doseInterpretResult` is the single public result contract for disease
#' interpretation workflows in DOSE.
#'
#' @slot result A data.frame of ranked interpretation targets.
#' @slot evidence A data.frame of long-form evidence rows linked to targets.
#' @slot query A data.frame describing the input query or queries.
#' @slot sources A data.frame describing source versions and checksums.
#' @slot parameters A named list of analysis parameters.
#' @slot explanation A named list storing optional explanation output.
#'
#' @name doseInterpretResult-class
#' @exportClass doseInterpretResult
NULL

.dose_result_fields <- c(
    "query_id", "target_id", "target_name", "target_type", "rank",
    "score", "score_type", "score_direction", "pvalue", "p.adjust",
    "evidence_count", "top_evidence_type", "source_count"
)

.dose_evidence_fields <- c(
    "evidence_id", "query_id", "target_id", "target_type", "evidence_type",
    "direction", "feature_id", "feature_name", "component_score",
    "source", "source_record_id", "evidence_path_id",
    "derived_from_evidence_id", "reference_id", "note"
)

.dose_query_fields <- c("query_id", "input_type", "organism")
.dose_source_fields <- c("source", "version", "checksum")
.dose_evidence_direction <- c("support", "missing", "conflict", "ambiguous")

.empty_result_df <- function() {
    data.frame(
        query_id = character(),
        target_id = character(),
        target_name = character(),
        target_type = character(),
        rank = integer(),
        score = numeric(),
        score_type = character(),
        score_direction = character(),
        pvalue = numeric(),
        p.adjust = numeric(),
        evidence_count = integer(),
        top_evidence_type = character(),
        source_count = integer(),
        stringsAsFactors = FALSE
    )
}

.empty_evidence_df <- function() {
    data.frame(
        evidence_id = character(),
        query_id = character(),
        target_id = character(),
        target_type = character(),
        evidence_type = character(),
        direction = character(),
        feature_id = character(),
        feature_name = character(),
        component_score = numeric(),
        source = character(),
        source_record_id = character(),
        evidence_path_id = character(),
        derived_from_evidence_id = character(),
        reference_id = character(),
        note = character(),
        stringsAsFactors = FALSE
    )
}

.empty_query_df <- function() {
    data.frame(
        query_id = character(),
        input_type = character(),
        organism = character(),
        stringsAsFactors = FALSE
    )
}

.empty_sources_df <- function() {
    data.frame(
        source = character(),
        version = character(),
        checksum = character(),
        stringsAsFactors = FALSE
    )
}

.normalize_contract_df <- function(x, empty) {
    if (is.null(x)) return(empty())

    x <- as.data.frame(x, stringsAsFactors = FALSE)
    if (nrow(x) == 0L) return(empty())
    x
}

.missing_contract_fields <- function(x, required) {
    setdiff(required, colnames(x))
}

.result_key <- function(x) {
    paste(x$query_id, x$target_type, x$target_id, sep = "\r")
}

methods::setClass(
    "doseInterpretResult",
    slots = c(
        result = "data.frame",
        evidence = "data.frame",
        query = "data.frame",
        sources = "data.frame",
        parameters = "list",
        explanation = "list"
    ),
    validity = function(object) {
        missing_result <- .missing_contract_fields(object@result, .dose_result_fields)
        if (length(missing_result)) {
            return(
                sprintf(
                    "`result` is missing required field(s): %s.",
                    paste(missing_result, collapse = ", ")
                )
            )
        }

        missing_evidence <- .missing_contract_fields(object@evidence, .dose_evidence_fields)
        if (length(missing_evidence)) {
            return(
                sprintf(
                    "`evidence` is missing required field(s): %s.",
                    paste(missing_evidence, collapse = ", ")
                )
            )
        }

        missing_query <- .missing_contract_fields(object@query, .dose_query_fields)
        if (length(missing_query)) {
            return(
                sprintf(
                    "`query` is missing required field(s): %s.",
                    paste(missing_query, collapse = ", ")
                )
            )
        }

        missing_sources <- .missing_contract_fields(object@sources, .dose_source_fields)
        if (length(missing_sources)) {
            return(
                sprintf(
                    "`sources` is missing required field(s): %s.",
                    paste(missing_sources, collapse = ", ")
                )
            )
        }

        if (anyDuplicated(.result_key(object@result))) {
            return(
                paste0(
                    "Each result row must be uniquely identified by query_id, ",
                    "target_type, and target_id."
                )
            )
        }

        evidence_id <- object@evidence$evidence_id
        if (anyDuplicated(evidence_id[nzchar(evidence_id)])) {
            return("`evidence_id` values must be unique.")
        }

        direction <- stats::na.omit(unique(object@evidence$direction))
        invalid_direction <- setdiff(direction, .dose_evidence_direction)
        if (length(invalid_direction)) {
            return(
                sprintf(
                    "Unsupported evidence `direction` value(s): %s.",
                    paste(invalid_direction, collapse = ", ")
                )
            )
        }

        if (nrow(object@evidence) > 0L) {
            if (nrow(object@result) == 0L) {
                return("Evidence rows require at least one corresponding result row.")
            }

            result_keys <- .result_key(object@result)
            evidence_keys <- .result_key(object@evidence)
            if (!all(evidence_keys %in% result_keys)) {
                return(
                    paste0(
                        "Every evidence row must link to a valid result row via ",
                        "query_id, target_type, and target_id."
                    )
                )
            }
        }

        TRUE
    }
)

#' Construct a canonical interpretation result
#'
#' @param result Ranked target rows.
#' @param evidence Long-form evidence rows.
#' @param query Query metadata rows.
#' @param sources Source metadata rows.
#' @param parameters Named list of analysis parameters.
#' @param explanation Named list of explanation output.
#'
#' @return A validated `doseInterpretResult` object.
#' @export
doseInterpretResult <- function(result = NULL,
                                evidence = NULL,
                                query = NULL,
                                sources = NULL,
                                parameters = list(),
                                explanation = list()) {
    object <- methods::new(
        "doseInterpretResult",
        result = .normalize_contract_df(result, .empty_result_df),
        evidence = .normalize_contract_df(evidence, .empty_evidence_df),
        query = .normalize_contract_df(query, .empty_query_df),
        sources = .normalize_contract_df(sources, .empty_sources_df),
        parameters = parameters,
        explanation = explanation
    )
    methods::validObject(object)
    object
}

#' Filter evidence rows from a canonical interpretation result
#'
#' @param x A `doseInterpretResult` object.
#' @param target_id Optional target identifier to keep.
#' @param type Optional evidence type to keep.
#' @param source Optional evidence source to keep.
#' @param direction Optional evidence direction to keep.
#'
#' @return A filtered evidence data.frame.
#' @export
methods::setGeneric(
    "evidence",
    function(x, target_id = NULL, type = NULL, source = NULL, direction = NULL) {
        methods::standardGeneric("evidence")
    }
)

#' @rdname evidence
#' @export
methods::setMethod(
    "evidence",
    "doseInterpretResult",
    function(x, target_id = NULL, type = NULL, source = NULL, direction = NULL) {
        out <- x@evidence
        if (!is.null(target_id)) out <- out[out$target_id %in% target_id, , drop = FALSE]
        if (!is.null(type)) out <- out[out$evidence_type %in% type, , drop = FALSE]
        if (!is.null(source)) out <- out[out$source %in% source, , drop = FALSE]
        if (!is.null(direction)) out <- out[out$direction %in% direction, , drop = FALSE]
        rownames(out) <- NULL
        out
    }
)

#' Coerce a canonical interpretation result to a data.frame
#'
#' @param x A `doseInterpretResult` object.
#' @param row.names,optional,... Passed to [base::as.data.frame()].
#'
#' @return A data.frame of ranked result rows.
#' @rdname doseInterpretResult-dataframe
#' @export
methods::setMethod(
    "as.data.frame",
    "doseInterpretResult",
    function(x, row.names = NULL, optional = FALSE, ...) {
        x@result
    }
)

#' Summarize a canonical interpretation result
#'
#' @param object A `doseInterpretResult` object.
#' @param ... Reserved for generic compatibility.
#'
#' @return A named list summarizing queries, targets, evidence, and sources.
#' @rdname doseInterpretResult-summary
#' @export
methods::setMethod(
    "summary",
    "doseInterpretResult",
    function(object, ...) {
        result <- object@result
        list(
            query_count = length(unique(object@query$query_id)),
            target_count = nrow(result),
            evidence_count = nrow(object@evidence),
            top_targets = utils::head(result[, c("target_id", "target_name", "score"), drop = FALSE], 5L),
            evidence_types = sort(unique(object@evidence$evidence_type)),
            sources = sort(unique(object@sources$source))
        )
    }
)

#' Show a canonical interpretation result
#'
#' @param object A `doseInterpretResult` object.
#'
#' @return The object is shown to the console and returned invisibly.
#' @rdname doseInterpretResult-show
#' @export
methods::setMethod(
    "show",
    "doseInterpretResult",
    function(object) {
        info <- summary(object)
        cat("doseInterpretResult\n")
        cat("  queries :", info$query_count, "\n")
        cat("  targets :", info$target_count, "\n")
        cat("  evidence:", info$evidence_count, "\n")
        if (nrow(info$top_targets) > 0L) {
            cat("  top target:", info$top_targets$target_name[[1]], "(", info$top_targets$target_id[[1]], ")\n")
        }
        if (length(info$sources) > 0L) {
            cat("  sources :", paste(info$sources, collapse = ", "), "\n")
        }
        invisible(object)
    }
)

.infer_interpret_input <- function(x) {
    if (is.list(x) && !is.data.frame(x) && !is.atomic(x)) {
        return("gene_set_list")
    }

    if (is.numeric(x) && !is.null(names(x))) {
        return("ranked_gene")
    }

    if (is.character(x)) {
        phenotype_idx <- grepl("^(HP:|MP:)", x)
        if (all(phenotype_idx)) return("phenotype")
        if (any(phenotype_idx)) {
            stop(
                paste0(
                    "Mixed identifier systems are not supported in `interpretDisease()`. ",
                    "Use only gene IDs or only phenotype IDs."
                ),
                call. = FALSE
            )
        }
        return("gene")
    }

    stop(
        paste0(
            "`interpretDisease()` expects a character gene vector, a named numeric ",
            "ranked gene vector, or a named list of gene vectors."
        ),
        call. = FALSE
    )
}

.validate_gene_set_list <- function(x) {
    if (!length(x)) {
        stop("`x` must contain at least one named gene set.", call. = FALSE)
    }
    if (is.null(names(x)) || any(!nzchar(names(x)))) {
        stop("`x` must be a named list of gene sets.", call. = FALSE)
    }
    if (anyDuplicated(names(x))) {
        stop("`x` must have unique gene-set names.", call. = FALSE)
    }

    invisible(lapply(x, function(ids) .validate_entrez_ids(ids, "gene set")))
}

.resolve_interpret_ontology <- function(ontology, organism) {
    if (identical(ontology, "auto")) {
        if (organism != "human") {
            stop(
                paste0(
                    "Cross-species disease interpretation from mouse input is not yet enabled. ",
                    "The first public path supports human gene input against HDO only."
                ),
                call. = FALSE
            )
        }
        return("HDO")
    }

    info <- .resolve_ontology_organism(ontology, organism)
    if (info$ontology != "HDO") {
        stop(
            paste0(
                "Only ontology = 'HDO' is currently enabled for `interpretDisease()`. ",
                "Phenotype-driven and mouse-model targets remain gated."
            ),
            call. = FALSE
        )
    }

    info$ontology
}

.validate_interpret_request <- function(x,
                                        input,
                                        organism,
                                        target,
                                        ontology,
                                        id_type) {
    inferred_input <- .infer_interpret_input(x)
    if (input == "auto") {
        input <- inferred_input
    } else if (!identical(input, inferred_input)) {
        stop(
            sprintf(
                "Declared `input = '%s'` is incompatible with the supplied object; inferred input is '%s'.",
                input,
                inferred_input
            ),
            call. = FALSE
        )
    }

    if (input == "phenotype") {
        stop(
            paste0(
                "`interpretDisease()` reserves phenotype input, but this path is not yet enabled ",
                "until the phenotype mapping, coverage, and provenance gates pass."
            ),
            call. = FALSE
        )
    }

    ontology <- .resolve_interpret_ontology(ontology, organism)

    if (!id_type %in% c("auto", "ENTREZID")) {
        stop(
            "`id_type` must be 'auto' or 'ENTREZID' for the current interpretation contract.",
            call. = FALSE
        )
    }

    if (input == "gene") {
        .validate_entrez_ids(x)
    } else if (input == "ranked_gene") {
        .validate_ranked_gene_list(x)
    } else if (input == "gene_set_list") {
        .validate_gene_set_list(x)
    }

    list(input = input, ontology = ontology)
}

.split_gene_ids <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(character())
    unique(unlist(strsplit(x, "/", fixed = TRUE), use.names = FALSE))
}

.build_sources_df <- function(ontology) {
    data.frame(
        source = ontology,
        version = "unknown",
        checksum = NA_character_,
        stringsAsFactors = FALSE
    )
}

.build_query_df <- function(input, organism, query_id = "q1", query_name = NULL) {
    out <- data.frame(
        query_id = query_id,
        input_type = input,
        organism = organism,
        stringsAsFactors = FALSE
    )

    if (!is.null(query_name)) {
        out$query_name <- query_name
    }

    out
}

.extract_cross_species_args <- function(extra_args, target) {
    required <- c("disease_model_file", "homology_file", "gene_pheno_file")
    arg_names <- names(extra_args)
    if (is.null(arg_names)) {
        arg_names <- rep("", length(extra_args))
        names(extra_args) <- arg_names
    }

    keep_idx <- !nzchar(arg_names) | !(arg_names %in% required)
    analysis_args <- extra_args[keep_idx]

    if (!(target %in% c("mouse_model", "both"))) {
        return(list(files = NULL, analysis_args = analysis_args))
    }

    missing <- required[vapply(
        required,
        function(name) is.null(extra_args[[name]]),
        logical(1)
    )]
    if (length(missing)) {
        stop(
            paste0(
                "`interpretDisease(..., target = \"", target, "\")` requires the supplied cross-species contract files: ",
                paste(required, collapse = ", "),
                "."
            ),
            call. = FALSE
        )
    }

    files <- extra_args[required]
    invalid <- required[vapply(
        files,
        function(path) !(is.character(path) && length(path) == 1L && !is.na(path) && nzchar(path)),
        logical(1)
    )]
    if (length(invalid)) {
        stop(
            sprintf(
                "`%s` must be a length-1 character path.",
                invalid[[1]]
            ),
            call. = FALSE
        )
    }

    list(files = files, analysis_args = analysis_args)
}

.extract_llm_args <- function(extra_args, explain) {
    llm_fields <- c("provider", "model", "language", "style", "generation_params", "llm_adapter")
    arg_names <- names(extra_args)
    if (is.null(arg_names)) {
        arg_names <- rep("", length(extra_args))
        names(extra_args) <- arg_names
    }

    llm_idx <- nzchar(arg_names) & arg_names %in% llm_fields
    llm_args <- extra_args[llm_idx]
    analysis_args <- extra_args[!llm_idx]

    if (!identical(explain, "llm")) {
        return(list(llm_args = list(), extra_args = analysis_args))
    }

    list(llm_args = llm_args, extra_args = analysis_args)
}

.empty_explanation <- function() {
    list(method = "none", text = character(), metadata = list())
}

.retag_interpret_object <- function(x,
                                    input,
                                    organism,
                                    query_id = "q1",
                                    query_name = NULL,
                                    target = NULL,
                                    query_gene = NULL) {
    original_query_id <- x@query$query_id[[1]]

    x@query <- .build_query_df(
        input = input,
        organism = organism,
        query_id = query_id,
        query_name = query_name
    )

    if (nrow(x@result)) {
        x@result$query_id <- query_id
    }
    if (nrow(x@evidence)) {
        x@evidence$query_id <- query_id
        if (!identical(query_id, original_query_id)) {
            x@evidence$evidence_id <- sprintf("%s::%s", query_id, x@evidence$evidence_id)
            x@evidence$evidence_path_id <- sprintf("%s::%s", query_id, x@evidence$evidence_path_id)

            keep_idx <- !is.na(x@evidence$derived_from_evidence_id) & nzchar(x@evidence$derived_from_evidence_id)
            x@evidence$derived_from_evidence_id[keep_idx] <- sprintf(
                "%s::%s",
                query_id,
                x@evidence$derived_from_evidence_id[keep_idx]
            )
        }
    }

    x@parameters$input <- input
    x@parameters$organism <- organism
    if (!is.null(target)) {
        x@parameters$target <- target
    }
    if (!is.null(query_gene)) {
        x@parameters$query_gene <- as.character(query_gene)
    }

    methods::validObject(x)
    x
}

.genes_for_mouse_model_input <- function(x, input, query_id = NULL) {
    if (identical(input, "gene")) {
        return(as.character(x))
    }

    if (identical(input, "ranked_gene")) {
        return(as.character(names(x)))
    }

    if (identical(input, "gene_set_list")) {
        return(as.character(x[[query_id]]))
    }

    stop("Mouse-model interpretation expects gene-derived input.", call. = FALSE)
}

.combine_interpret_objects <- function(parts, parameters) {
    parts <- Filter(Negate(is.null), parts)
    if (!length(parts)) {
        return(doseInterpretResult(parameters = parameters, explanation = .empty_explanation()))
    }

    result_df <- do.call(rbind, lapply(parts, function(part) part@result))
    evidence_df <- do.call(rbind, lapply(parts, function(part) part@evidence))
    query_df <- unique(do.call(rbind, lapply(parts, function(part) part@query)))
    sources_df <- unique(do.call(rbind, lapply(parts, function(part) part@sources)))

    result_df <- .normalize_contract_df(result_df, .empty_result_df)
    evidence_df <- .normalize_contract_df(evidence_df, .empty_evidence_df)
    query_df <- .normalize_contract_df(query_df, .empty_query_df)
    sources_df <- .normalize_contract_df(sources_df, .empty_sources_df)

    if (nrow(result_df)) {
        query_levels <- unique(query_df$query_id)
        target_levels <- c("disease", "mouse_model")
        ord <- order(
            match(result_df$query_id, query_levels),
            match(result_df$target_type, target_levels),
            result_df$rank,
            result_df$target_name,
            result_df$target_id
        )
        result_df <- result_df[ord, , drop = FALSE]
        rownames(result_df) <- NULL
    }

    if (nrow(evidence_df)) {
        query_levels <- unique(query_df$query_id)
        target_levels <- c("disease", "mouse_model")
        ord <- order(
            match(evidence_df$query_id, query_levels),
            match(evidence_df$target_type, target_levels),
            evidence_df$target_id,
            evidence_df$evidence_type,
            evidence_df$evidence_id
        )
        evidence_df <- evidence_df[ord, , drop = FALSE]
        rownames(evidence_df) <- NULL
    }

    rownames(query_df) <- NULL
    rownames(sources_df) <- NULL

    doseInterpretResult(
        result = result_df,
        evidence = evidence_df,
        query = query_df,
        sources = sources_df,
        parameters = parameters,
        explanation = .empty_explanation()
    )
}

.slice_interpret_object <- function(x, query_id) {
    doseInterpretResult(
        result = x@result[x@result$query_id == query_id, , drop = FALSE],
        evidence = x@evidence[x@evidence$query_id == query_id, , drop = FALSE],
        query = x@query[x@query$query_id == query_id, , drop = FALSE],
        sources = x@sources,
        parameters = x@parameters,
        explanation = .empty_explanation()
    )
}

.serialize_llm_payload <- function(payload) {
    paste(utils::capture.output(dput(payload)), collapse = "\n")
}

.md5_text <- function(text) {
    path <- tempfile("dose-llm-", fileext = ".txt")
    on.exit(unlink(path), add = TRUE)
    writeLines(enc2utf8(text), path, useBytes = TRUE)
    unname(tools::md5sum(path)[[1]])
}

.build_llm_payload <- function(x, language, style) {
    list(
        query = x@query,
        result = x@result,
        evidence = x@evidence,
        sources = x@sources,
        parameters = x@parameters,
        explanation_request = list(
            language = language,
            style = style
        )
    )
}

.build_llm_prompt <- function(payload_text, language, style) {
    paste(
        "You may only use the evidence supplied below.",
        "Do not add external facts.",
        "Do not upgrade weak evidence to strong evidence.",
        "State uncertainty when evidence is incomplete or ambiguous.",
        sprintf("Write the explanation in language: %s.", language),
        sprintf("Use this style: %s.", style),
        "Return a structured object with `text`, `cited_evidence_ids`, `unsupported_claims`, and `metadata`.",
        "",
        payload_text,
        sep = "\n"
    )
}

.default_llm_adapter <- function(request) {
    if (!requireNamespace("aisdk", quietly = TRUE)) {
        stop(
            paste0(
                "`explainDisease(..., method = \"llm\")` requires the optional `aisdk` package ",
                "or an explicit `llm_adapter`."
            ),
            call. = FALSE
        )
    }

    ns <- asNamespace("aisdk")
    fun_name <- c("generateText", "generate_text")
    fun_name <- fun_name[vapply(fun_name, exists, logical(1), envir = ns, inherits = FALSE)]
    if (!length(fun_name)) {
        stop(
            "The installed `aisdk` package does not expose a supported text-generation entry point.",
            call. = FALSE
        )
    }

    fun <- get(fun_name[[1]], envir = ns, inherits = FALSE)
    args <- c(
        list(
            prompt = request$prompt,
            provider = request$provider,
            model = request$model
        ),
        request$generation_params
    )
    do.call(fun, args)
}

.normalize_llm_response <- function(response, valid_evidence_ids) {
    if (!is.list(response)) {
        stop("The LLM adapter must return a list response.", call. = FALSE)
    }

    if (is.null(response$text) || !is.character(response$text) || length(response$text) != 1L || !nzchar(response$text)) {
        stop("The LLM adapter response must include one non-empty `text` field.", call. = FALSE)
    }

    cited_evidence_ids <- response$cited_evidence_ids
    if (is.null(cited_evidence_ids)) {
        cited_evidence_ids <- character()
    }
    cited_evidence_ids <- unique(stats::na.omit(as.character(cited_evidence_ids)))

    invalid_ids <- setdiff(cited_evidence_ids, valid_evidence_ids)
    if (length(invalid_ids)) {
        stop(
            sprintf(
                "The LLM adapter cited unsupported evidence_id value(s): %s.",
                paste(invalid_ids, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    unsupported_claims <- response$unsupported_claims
    if (is.null(unsupported_claims)) {
        unsupported_claims <- character()
    }
    unsupported_claims <- unique(stats::na.omit(as.character(unsupported_claims)))
    unsupported_claims <- unsupported_claims[nzchar(unsupported_claims)]
    if (length(unsupported_claims)) {
        stop(
            sprintf(
                "The LLM adapter returned unsupported claims: %s.",
                paste(unsupported_claims, collapse = " | ")
            ),
            call. = FALSE
        )
    }

    if (length(valid_evidence_ids) && !length(cited_evidence_ids)) {
        stop(
            "The LLM adapter response must cite at least one `evidence_id` when evidence rows are available.",
            call. = FALSE
        )
    }

    list(
        text = response$text,
        cited_evidence_ids = cited_evidence_ids,
        unsupported_claims = unsupported_claims,
        metadata = if (is.null(response$metadata)) list() else response$metadata
    )
}

.run_llm_explanation <- function(x,
                                 provider = NULL,
                                 model = NULL,
                                 language = "en",
                                 style = "brief",
                                 generation_params = list(),
                                 llm_adapter = NULL) {
    if (!is.list(generation_params)) {
        stop("`generation_params` must be a list.", call. = FALSE)
    }

    if (is.null(llm_adapter) &&
        is.null(getOption("DOSE.llm_adapter")) &&
        (is.null(provider) || is.null(model))) {
        stop(
            paste0(
                "`explainDisease(..., method = \"llm\")` requires the optional `aisdk` package ",
                "with configured `provider` and `model`, or an explicit `llm_adapter`."
            ),
            call. = FALSE
        )
    }

    adapter <- llm_adapter
    if (is.null(adapter)) {
        adapter <- getOption("DOSE.llm_adapter")
    }
    if (is.null(adapter)) {
        adapter <- .default_llm_adapter
    }
    if (!is.function(adapter)) {
        stop("`llm_adapter` must be a function.", call. = FALSE)
    }

    query_ids <- unique(x@query$query_id)
    text <- stats::setNames(character(length(query_ids)), query_ids)
    query_meta <- vector("list", length(query_ids))
    names(query_meta) <- query_ids

    for (query_id in query_ids) {
        slice <- .slice_interpret_object(x, query_id)
        payload <- .build_llm_payload(slice, language = language, style = style)
        payload_text <- .serialize_llm_payload(payload)
        prompt <- .build_llm_prompt(payload_text, language = language, style = style)

        request <- list(
            prompt = prompt,
            payload = payload,
            provider = provider,
            model = model,
            generation_params = generation_params
        )

        response <- tryCatch(
            adapter(request),
            error = function(e) {
                stop(
                    sprintf("The LLM adapter failed for query '%s': %s", query_id, conditionMessage(e)),
                    call. = FALSE
                )
            }
        )

        normalized <- .normalize_llm_response(
            response = response,
            valid_evidence_ids = slice@evidence$evidence_id
        )
        response_text <- .serialize_llm_payload(response)

        text[[query_id]] <- normalized$text
        query_meta[[query_id]] <- list(
            provider = provider,
            model = model,
            evidence_hash = .md5_text(paste(slice@evidence$evidence_id, collapse = "\n")),
            payload_hash = .md5_text(payload_text),
            prompt_hash = .md5_text(prompt),
            response_hash = .md5_text(response_text),
            generation_params = generation_params,
            response_metadata = normalized$metadata,
            cited_evidence_ids = normalized$cited_evidence_ids
        )
    }

    x@explanation <- list(
        method = "llm",
        text = text,
        metadata = list(
            provider = provider,
            model = model,
            language = language,
            style = style,
            queries = query_meta
        )
    )
    methods::validObject(x)
    x
}

.collapse_template_features <- function(feature_id, n = 3L) {
    feature_id <- unique(stats::na.omit(as.character(feature_id)))
    feature_id <- feature_id[nzchar(feature_id)]
    if (!length(feature_id)) {
        return("no named features")
    }

    shown <- utils::head(feature_id, n)
    out <- paste(shown, collapse = ", ")
    if (length(feature_id) > n) {
        out <- paste0(out, ", ...")
    }
    out
}

.template_explanation_text <- function(x, query_id, max_features = 3L) {
    query_tbl <- x@query[x@query$query_id == query_id, , drop = FALSE]
    result_tbl <- x@result[x@result$query_id == query_id, , drop = FALSE]
    source_label <- paste(unique(x@sources$source), collapse = ", ")

    if (!nrow(result_tbl)) {
        return(
            sprintf(
                "Query '%s' returns no ranked target after the current thresholds, so there is no evidence-grounded summary to report from %s.",
                query_id,
                source_label
            )
        )
    }

    top <- result_tbl[1, , drop = FALSE]
    top_evidence <- x@evidence[
        x@evidence$query_id == query_id &
            x@evidence$target_id == top$target_id,
        ,
        drop = FALSE
    ]

    intro <- sprintf(
        "Query '%s' prioritizes %s (%s) from %s at rank %s (score %.3f; adjusted p = %s).",
        query_id,
        top$target_name,
        top$target_id,
        source_label,
        top$rank,
        top$score,
        format(top$p.adjust, digits = 3, scientific = TRUE)
    )

    evidence_line <- if (!nrow(top_evidence)) {
        "No traceable evidence rows are attached to the top target yet."
    } else {
        evidence_types <- paste(unique(top_evidence$evidence_type), collapse = ", ")
        feature_summary <- .collapse_template_features(top_evidence$feature_id, max_features)
        sprintf(
            "It is currently backed by %s %s evidence row(s), led by %s.",
            nrow(top_evidence),
            evidence_types,
            feature_summary
        )
    }

    other_targets <- result_tbl[-1, , drop = FALSE]
    tail_line <- if (!nrow(other_targets)) {
        "No additional ranked targets are attached to this query."
    } else {
        other_labels <- utils::head(
            sprintf("%s (%s)", other_targets$target_name, other_targets$target_id),
            2L
        )
        sprintf(
            "Other ranked targets include %s.",
            paste(other_labels, collapse = "; ")
        )
    }

    if ("query_name" %in% colnames(query_tbl) &&
        nzchar(query_tbl$query_name[[1]]) &&
        !identical(query_tbl$query_name[[1]], query_id)) {
        intro <- sub(
            sprintf("Query '%s'", query_id),
            sprintf("Query '%s' (%s)", query_id, query_tbl$query_name[[1]]),
            intro,
            fixed = TRUE
        )
    }

    paste(intro, evidence_line, tail_line)
}

.enrich_result_to_interpret <- function(res,
                                        input,
                                        organism,
                                        ontology,
                                        target,
                                        id_type,
                                        orthology,
                                        explain,
                                        top,
                                        query_id = "q1",
                                        query_name = NULL) {
    result_df <- .empty_result_df()
    evidence_df <- .empty_evidence_df()

    if (!is.null(res)) {
        direct <- as.data.frame(res)
        if (nrow(direct) > 0L) {
            direct <- direct[seq_len(min(nrow(direct), top)), , drop = FALSE]

            result_df <- data.frame(
                query_id = rep(query_id, nrow(direct)),
                target_id = as.character(direct$ID),
                target_name = as.character(direct$Description),
                target_type = rep("disease", nrow(direct)),
                rank = seq_len(nrow(direct)),
                score = as.numeric(direct$FoldEnrichment),
                score_type = rep("fold_enrichment", nrow(direct)),
                score_direction = rep("higher", nrow(direct)),
                pvalue = as.numeric(direct$pvalue),
                p.adjust = as.numeric(direct$p.adjust),
                evidence_count = as.integer(direct$Count),
                top_evidence_type = rep("gene", nrow(direct)),
                source_count = rep(1L, nrow(direct)),
                stringsAsFactors = FALSE
            )

            evidence_rows <- lapply(seq_len(nrow(direct)), function(i) {
                genes <- .split_gene_ids(direct$geneID[[i]])
                if (!length(genes)) return(NULL)

                data.frame(
                    evidence_id = sprintf(
                        "%s::%s::gene::%s",
                        query_id,
                        as.character(direct$ID[[i]]),
                        genes
                    ),
                    query_id = rep(query_id, length(genes)),
                    target_id = rep(as.character(direct$ID[[i]]), length(genes)),
                    target_type = rep("disease", length(genes)),
                    evidence_type = rep("gene", length(genes)),
                    direction = rep("support", length(genes)),
                    feature_id = genes,
                    feature_name = genes,
                    component_score = rep(NA_real_, length(genes)),
                    source = rep(ontology, length(genes)),
                    source_record_id = rep(as.character(direct$ID[[i]]), length(genes)),
                    evidence_path_id = rep(
                        sprintf("%s::gene_overlap", as.character(direct$ID[[i]])),
                        length(genes)
                    ),
                    derived_from_evidence_id = rep(NA_character_, length(genes)),
                    reference_id = rep(NA_character_, length(genes)),
                    note = rep("Input gene overlaps ontology annotation.", length(genes)),
                    stringsAsFactors = FALSE
                )
            })

            evidence_rows <- Filter(Negate(is.null), evidence_rows)
            if (length(evidence_rows)) {
                evidence_df <- do.call(rbind, evidence_rows)
                rownames(evidence_df) <- NULL
            }
        }
    }

    doseInterpretResult(
        result = result_df,
        evidence = evidence_df,
        query = .build_query_df(
            input = input,
            organism = organism,
            query_id = query_id,
            query_name = query_name
        ),
        sources = .build_sources_df(ontology),
        parameters = list(
            input = input,
            organism = organism,
            target = target,
            ontology = ontology,
            id_type = id_type,
            orthology = orthology,
            explain = explain,
            top = top
        ),
        explanation = .empty_explanation()
    )
}

.gsea_result_to_interpret <- function(res,
                                      input,
                                      organism,
                                      ontology,
                                      target,
                                      id_type,
                                      orthology,
                                      explain,
                                      top,
                                      geneList,
                                      query_id = "q1",
                                      query_name = NULL) {
    result_df <- .empty_result_df()
    evidence_df <- .empty_evidence_df()

    if (!is.null(res)) {
        direct <- as.data.frame(res)
        if (nrow(direct) > 0L) {
            direct <- direct[seq_len(min(nrow(direct), top)), , drop = FALSE]

            result_df <- data.frame(
                query_id = rep(query_id, nrow(direct)),
                target_id = as.character(direct$ID),
                target_name = as.character(direct$Description),
                target_type = rep("disease", nrow(direct)),
                rank = seq_len(nrow(direct)),
                score = as.numeric(direct$NES),
                score_type = rep("normalized_enrichment_score", nrow(direct)),
                score_direction = rep("higher", nrow(direct)),
                pvalue = as.numeric(direct$pvalue),
                p.adjust = as.numeric(direct$p.adjust),
                evidence_count = vapply(
                    direct$core_enrichment,
                    function(x) length(.split_gene_ids(x)),
                    integer(1)
                ),
                top_evidence_type = rep("ranked_gene", nrow(direct)),
                source_count = rep(1L, nrow(direct)),
                stringsAsFactors = FALSE
            )

            evidence_rows <- lapply(seq_len(nrow(direct)), function(i) {
                genes <- .split_gene_ids(direct$core_enrichment[[i]])
                if (!length(genes)) return(NULL)

                gene_scores <- unname(geneList[genes])

                data.frame(
                    evidence_id = sprintf(
                        "%s::%s::ranked_gene::%s",
                        query_id,
                        as.character(direct$ID[[i]]),
                        genes
                    ),
                    query_id = rep(query_id, length(genes)),
                    target_id = rep(as.character(direct$ID[[i]]), length(genes)),
                    target_type = rep("disease", length(genes)),
                    evidence_type = rep("ranked_gene", length(genes)),
                    direction = rep("support", length(genes)),
                    feature_id = genes,
                    feature_name = genes,
                    component_score = as.numeric(gene_scores),
                    source = rep(ontology, length(genes)),
                    source_record_id = rep(as.character(direct$ID[[i]]), length(genes)),
                    evidence_path_id = rep(
                        sprintf("%s::leading_edge", as.character(direct$ID[[i]])),
                        length(genes)
                    ),
                    derived_from_evidence_id = rep(NA_character_, length(genes)),
                    reference_id = rep(NA_character_, length(genes)),
                    note = rep("Leading-edge ranked gene support.", length(genes)),
                    stringsAsFactors = FALSE
                )
            })

            evidence_rows <- Filter(Negate(is.null), evidence_rows)
            if (length(evidence_rows)) {
                evidence_df <- do.call(rbind, evidence_rows)
                rownames(evidence_df) <- NULL
            }
        }
    }

    doseInterpretResult(
        result = result_df,
        evidence = evidence_df,
        query = .build_query_df(
            input = input,
            organism = organism,
            query_id = query_id,
            query_name = query_name
        ),
        sources = .build_sources_df(ontology),
        parameters = list(
            input = input,
            organism = organism,
            target = target,
            ontology = ontology,
            id_type = id_type,
            orthology = orthology,
            explain = explain,
            top = top
        ),
        explanation = .empty_explanation()
    )
}

#' High-level disease interpretation entry point
#'
#' This function defines the public interpretation contract and performs
#' validation before the first disease-ranking implementation is enabled.
#'
#' @param x Query input.
#' @param input One of `"auto"`, `"gene"`, `"ranked_gene"`,
#'   `"gene_set_list"`, or `"phenotype"`.
#' @param organism One of `"human"` or `"mouse"`.
#' @param target One of `"disease"`, `"mouse_model"`, or `"both"`.
#' @param ontology Ontology to use. `"auto"` currently resolves to `"HDO"` for
#'   the supported human disease path.
#' @param id_type Identifier type. The current implementation accepts `"auto"`
#'   and `"ENTREZID"`.
#' @param orthology One of `"one_to_one"` or `"all"`.
#' @param explain One of `"none"`, `"template"`, or `"llm"`.
#' @param top Maximum number of targets to return.
#' @param ... Additional arguments passed to the underlying enrichment
#'   functions. Cross-species `target = "mouse_model"` and `target = "both"`
#'   also require `disease_model_file`, `homology_file`, and `gene_pheno_file`.
#'
#' @return A `doseInterpretResult` once the ranking path is implemented.
#' @export
interpretDisease <- function(
    x,
    input = c("auto", "gene", "ranked_gene", "gene_set_list", "phenotype"),
    organism = c("human", "mouse"),
    target = c("disease", "mouse_model", "both"),
    ontology = "auto",
    id_type = "auto",
    orthology = c("one_to_one", "all"),
    explain = c("none", "template", "llm"),
    top = 50,
    ...
) {
    extra_args <- list(...)
    input <- match.arg(input)
    organism <- match.arg(organism)
    target <- match.arg(target)
    orthology <- match.arg(orthology)
    explain <- match.arg(explain)

    if (length(top) != 1L || is.na(top) || !is.numeric(top) || top < 1) {
        stop("`top` must be one positive number.", call. = FALSE)
    }

    request <- .validate_interpret_request(
        x = x,
        input = input,
        organism = organism,
        target = target,
        ontology = ontology,
        id_type = id_type
    )

    input <- request$input
    ontology <- request$ontology
    llm_args <- .extract_llm_args(extra_args, explain)
    extra_args <- llm_args$extra_args
    resolved_args <- .extract_cross_species_args(extra_args, target)
    cross_species_files <- resolved_args$files
    analysis_args <- resolved_args$analysis_args

    if (!(input %in% c("gene", "ranked_gene", "gene_set_list") &&
          identical(organism, "human"))) {
        stop(
            paste0(
                "The current interpretation implementation supports human gene vectors, ranked gene vectors, ",
                "and named gene-set lists from human input only. ",
                "Phenotype input and mouse-source queries remain gated."
            ),
            call. = FALSE
        )
    }

    if (identical(input, "gene")) {
        disease_out <- NULL
        model_out <- NULL

        if (target %in% c("disease", "both")) {
            res <- do.call(
                enrichDisease,
                c(
                    list(
                        gene = x,
                        organism = organism,
                        ontology = ontology
                    ),
                    analysis_args
                )
            )

            disease_out <- .enrich_result_to_interpret(
                res = res,
                input = input,
                organism = organism,
                ontology = ontology,
                target = "disease",
                id_type = id_type,
                orthology = orthology,
                explain = explain,
                top = top
            )
        }

        if (target %in% c("mouse_model", "both")) {
            model_out <- rankMouseModels(
                gene = .genes_for_mouse_model_input(x, input),
                disease_model_file = cross_species_files$disease_model_file,
                homology_file = cross_species_files$homology_file,
                gene_pheno_file = cross_species_files$gene_pheno_file,
                top = top
            )
            model_out <- .retag_interpret_object(
                model_out,
                input = input,
                organism = organism,
                target = "mouse_model",
                query_gene = .genes_for_mouse_model_input(x, input)
            )
        }

        out <- if (identical(target, "disease")) {
            disease_out
        } else if (identical(target, "mouse_model")) {
            model_out
        } else {
            .combine_interpret_objects(
                list(disease_out, model_out),
                parameters = c(
                    list(
                        input = input,
                        organism = organism,
                        target = target,
                        ontology = ontology,
                        id_type = id_type,
                        orthology = orthology,
                        explain = explain,
                        top = top
                    ),
                    lapply(cross_species_files, normalizePath, mustWork = TRUE)
                )
            )
        }

        if (identical(explain, "template")) {
            return(explainDisease(out, method = "template"))
        }
        if (identical(explain, "llm")) {
            return(do.call(
                explainDisease,
                c(
                    list(x = out, method = "llm"),
                    llm_args$llm_args
                )
            ))
        }

        return(out)
    }

    if (identical(input, "gene_set_list")) {
        parts <- lapply(names(x), function(query_id) {
            disease_part <- NULL
            model_part <- NULL

            if (target %in% c("disease", "both")) {
                res <- do.call(
                    enrichDisease,
                    c(
                        list(
                            gene = x[[query_id]],
                            organism = organism,
                            ontology = ontology
                        ),
                        analysis_args
                    )
                )

                disease_part <- .enrich_result_to_interpret(
                    res = res,
                    input = input,
                    organism = organism,
                    ontology = ontology,
                    target = "disease",
                    id_type = id_type,
                    orthology = orthology,
                    explain = explain,
                    top = top,
                    query_id = query_id,
                    query_name = query_id
                )
            }

            if (target %in% c("mouse_model", "both")) {
                model_part <- rankMouseModels(
                    gene = .genes_for_mouse_model_input(x, input, query_id = query_id),
                    disease_model_file = cross_species_files$disease_model_file,
                    homology_file = cross_species_files$homology_file,
                    gene_pheno_file = cross_species_files$gene_pheno_file,
                    top = top
                )
                model_part <- .retag_interpret_object(
                    model_part,
                    input = input,
                    organism = organism,
                    query_id = query_id,
                    query_name = query_id,
                    target = "mouse_model",
                    query_gene = .genes_for_mouse_model_input(x, input, query_id = query_id)
                )
            }

            if (identical(target, "disease")) {
                return(disease_part)
            }
            if (identical(target, "mouse_model")) {
                return(model_part)
            }

            .combine_interpret_objects(
                list(disease_part, model_part),
                parameters = list()
            )
        })

        out <- .combine_interpret_objects(
            parts,
            parameters = c(
                list(
                    input = input,
                    organism = organism,
                    target = target,
                    ontology = ontology,
                    id_type = id_type,
                    orthology = orthology,
                    explain = explain,
                    top = top
                ),
                if (is.null(cross_species_files)) list() else lapply(cross_species_files, normalizePath, mustWork = TRUE)
            )
        )

        if (identical(explain, "template")) {
            return(explainDisease(out, method = "template"))
        }
        if (identical(explain, "llm")) {
            return(do.call(
                explainDisease,
                c(
                    list(x = out, method = "llm"),
                    llm_args$llm_args
                )
            ))
        }

        return(out)
    }

    disease_out <- NULL
    model_out <- NULL

    if (target %in% c("disease", "both")) {
        res <- do.call(
            gseDisease,
            c(
                list(
                    geneList = x,
                    organism = organism,
                    ontology = ontology
                ),
                analysis_args
            )
        )

        disease_out <- .gsea_result_to_interpret(
            res = res,
            input = input,
            organism = organism,
            ontology = ontology,
            target = "disease",
            id_type = id_type,
            orthology = orthology,
            explain = explain,
            top = top,
            geneList = x
        )
    }

    if (target %in% c("mouse_model", "both")) {
        model_out <- rankMouseModels(
            gene = .genes_for_mouse_model_input(x, input),
            disease_model_file = cross_species_files$disease_model_file,
            homology_file = cross_species_files$homology_file,
            gene_pheno_file = cross_species_files$gene_pheno_file,
            top = top
        )
        model_out <- .retag_interpret_object(
            model_out,
            input = input,
            organism = organism,
            target = "mouse_model",
            query_gene = .genes_for_mouse_model_input(x, input)
        )
    }

    out <- if (identical(target, "disease")) {
        disease_out
    } else if (identical(target, "mouse_model")) {
        model_out
    } else {
        .combine_interpret_objects(
            list(disease_out, model_out),
            parameters = c(
                list(
                    input = input,
                    organism = organism,
                    target = target,
                    ontology = ontology,
                    id_type = id_type,
                    orthology = orthology,
                    explain = explain,
                    top = top
                ),
                lapply(cross_species_files, normalizePath, mustWork = TRUE)
            )
        )
    }

    if (identical(explain, "template")) {
        return(explainDisease(out, method = "template"))
    }
    if (identical(explain, "llm")) {
        return(do.call(
            explainDisease,
            c(
                list(x = out, method = "llm"),
                llm_args$llm_args
            )
        ))
    }

    out
}

#' Add evidence-grounded explanation text to a canonical interpretation result
#'
#' @param x A `doseInterpretResult` object.
#' @param method Explanation method. `"template"` is available offline.
#' @param max_features Maximum number of evidence features to mention per query.
#' @param provider Optional provider name passed to the LLM adapter.
#' @param model Optional model name passed to the LLM adapter.
#' @param language Explanation language for the optional LLM adapter.
#' @param style Explanation style for the optional LLM adapter.
#' @param generation_params Optional generation parameters passed to the LLM
#'   adapter.
#' @param llm_adapter Optional adapter function for deterministic testing or
#'   custom provider integration.
#'
#' @return A `doseInterpretResult` with the `explanation` slot populated.
#' @export
explainDisease <- function(x,
                           method = c("template", "llm"),
                           max_features = 3L,
                           provider = NULL,
                           model = NULL,
                           language = c("en", "zh"),
                           style = c("brief", "research", "cautious"),
                           generation_params = list(),
                           llm_adapter = NULL) {
    method <- match.arg(method)
    language <- match.arg(language)
    style <- match.arg(style)

    if (!methods::is(x, "doseInterpretResult")) {
        stop("`x` must be a doseInterpretResult object.", call. = FALSE)
    }

    if (length(max_features) != 1L || is.na(max_features) || max_features < 1L) {
        stop("`max_features` must be one positive number.", call. = FALSE)
    }

    if (identical(method, "llm")) {
        return(.run_llm_explanation(
            x = x,
            provider = provider,
            model = model,
            language = language,
            style = style,
            generation_params = generation_params,
            llm_adapter = llm_adapter
        ))
    }

    query_ids <- unique(x@query$query_id)
    text <- stats::setNames(
        vapply(
            query_ids,
            function(query_id) .template_explanation_text(x, query_id, max_features),
            character(1)
        ),
        query_ids
    )

    x@explanation <- list(
        method = method,
        text = text,
        metadata = list(max_features = as.integer(max_features))
    )
    methods::validObject(x)
    x
}
