test_that("interpretDisease maps ranked gene GSEA into the canonical contract", {
    data(geneList, package = "DOSE")

    set.seed(1)
    direct <- as.data.frame(gseDisease(
        geneList,
        organism = "human",
        ontology = "HDO",
        verbose = FALSE
    ))

    set.seed(1)
    x <- interpretDisease(
        geneList,
        organism = "human",
        ontology = "HDO",
        verbose = FALSE
    )
    result_tbl <- as.data.frame(x)

    expect_s4_class(x, "doseInterpretResult")
    expect_true(nrow(result_tbl) > 0)
    expect_equal(result_tbl$target_id, direct$ID[seq_len(nrow(result_tbl))])
    expect_equal(result_tbl$target_name, direct$Description[seq_len(nrow(result_tbl))])
    expect_equal(result_tbl$score, direct$NES[seq_len(nrow(result_tbl))])
    expect_equal(result_tbl$p.adjust, direct$p.adjust[seq_len(nrow(result_tbl))])
    expect_true(all(result_tbl$target_type == "disease"))
    expect_true(all(result_tbl$top_evidence_type == "ranked_gene"))

    evidence_tbl <- evidence(x)
    expect_true(nrow(evidence_tbl) > 0)
    expect_true(all(evidence_tbl$evidence_type == "ranked_gene"))
    expect_true(all(evidence_tbl$target_id %in% result_tbl$target_id))
    expect_true(all(evidence_tbl$feature_id %in% names(geneList)))
    expect_true(all(is.finite(evidence_tbl$component_score)))
})

test_that("ranked gene interpretation can return an empty canonical object", {
    data(geneList, package = "DOSE")

    set.seed(1)
    x <- interpretDisease(
        geneList,
        organism = "human",
        ontology = "HDO",
        pvalueCutoff = 0,
        verbose = FALSE
    )

    expect_s4_class(x, "doseInterpretResult")
    expect_equal(nrow(as.data.frame(x)), 0)
    expect_equal(nrow(evidence(x)), 0)
    expect_equal(x@query$input_type[[1]], "ranked_gene")
})
