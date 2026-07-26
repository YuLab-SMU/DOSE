test_that("interpretDisease maps human gene enrichment into the canonical contract", {
    genes <- c("1", "2", "9", "10")

    direct <- as.data.frame(enrichDisease(genes, organism = "human", ontology = "HDO"))
    x <- interpretDisease(genes)
    result_tbl <- as.data.frame(x)

    expect_s4_class(x, "doseInterpretResult")
    expect_true(nrow(result_tbl) > 0)
    expect_equal(result_tbl$target_id, direct$ID[seq_len(nrow(result_tbl))])
    expect_equal(result_tbl$target_name, direct$Description[seq_len(nrow(result_tbl))])
    expect_equal(result_tbl$p.adjust, direct$p.adjust[seq_len(nrow(result_tbl))])
    expect_true(all(result_tbl$target_type == "disease"))
    expect_true(all(result_tbl$top_evidence_type == "gene"))

    evidence_tbl <- evidence(x)
    expect_true(nrow(evidence_tbl) > 0)
    expect_true(all(evidence_tbl$evidence_type == "gene"))
    expect_true(all(evidence_tbl$target_id %in% result_tbl$target_id))
    expect_true(all(evidence_tbl$feature_id %in% genes))
    expect_true(all(evidence_tbl$source == "HDO"))
})

test_that("interpretDisease returns an empty canonical object when enrichment has no hits", {
    x <- interpretDisease(
        c("1", "2", "9", "10"),
        pvalueCutoff = 0,
        qvalueCutoff = 0
    )

    expect_s4_class(x, "doseInterpretResult")
    expect_equal(nrow(as.data.frame(x)), 0)
    expect_equal(nrow(evidence(x)), 0)
    expect_equal(x@query$input_type[[1]], "gene")
    expect_equal(x@sources$source[[1]], "HDO")
})

test_that("explanation remains gated in the first human gene interpretation slice", {
    expect_error(
        interpretDisease(c("1", "2", "9", "10"), explain = "template"),
        "explain != \"none\".*not yet enabled"
    )
})
