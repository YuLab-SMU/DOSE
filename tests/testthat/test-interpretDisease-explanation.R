test_that("template explanation works for ranked gene interpretation", {
    data(geneList, package = "DOSE")

    set.seed(1)
    x <- interpretDisease(
        geneList,
        organism = "human",
        ontology = "HDO",
        explain = "template",
        verbose = FALSE
    )

    top <- as.data.frame(x)[1, , drop = FALSE]
    top_evidence <- evidence(x, target_id = top$target_id)[1, , drop = FALSE]

    expect_equal(x@explanation$method, "template")
    expect_true(grepl(top$target_name[[1]], x@explanation$text[[1]], fixed = TRUE))
    expect_true(grepl(top$target_id[[1]], x@explanation$text[[1]], fixed = TRUE))
    expect_true(grepl(top_evidence$feature_id[[1]], x@explanation$text[[1]], fixed = TRUE))
})

test_that("template explanation returns one summary per query for gene-set lists", {
    gene_sets <- list(
        cluster_a = c("1", "2", "9", "10"),
        cluster_b = c("1", "2", "9", "10", "11")
    )

    x <- interpretDisease(
        gene_sets,
        organism = "human",
        ontology = "HDO",
        explain = "template"
    )

    expect_equal(sort(names(x@explanation$text)), sort(names(gene_sets)))
    expect_true(all(nzchar(unname(x@explanation$text))))
})

test_that("template explanation surfaces missing evidence explicitly", {
    x <- doseInterpretResult(
        result = data.frame(
            query_id = "q_missing",
            target_id = "DOID:399",
            target_name = "tuberculosis",
            target_type = "disease",
            rank = 1L,
            score = 2,
            score_type = "fold_enrichment",
            score_direction = "higher",
            pvalue = 0.001,
            p.adjust = 0.01,
            evidence_count = 0L,
            top_evidence_type = "gene",
            source_count = 1L,
            stringsAsFactors = FALSE
        ),
        evidence = data.frame(
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
        ),
        query = data.frame(
            query_id = "q_missing",
            input_type = "gene",
            organism = "human",
            stringsAsFactors = FALSE
        ),
        sources = data.frame(
            source = "HDO",
            version = "unknown",
            checksum = NA_character_,
            stringsAsFactors = FALSE
        )
    )

    x <- explainDisease(x, method = "template")

    expect_match(
        x@explanation$text[["q_missing"]],
        "No traceable evidence rows are attached to the top target yet.",
        fixed = TRUE
    )
})

test_that("llm explanation adapter requires an optional provider bridge", {
    x <- interpretDisease(c("1", "2", "9", "10"))

    expect_error(
        explainDisease(x, method = "llm"),
        "requires the optional `aisdk` package"
    )
})
