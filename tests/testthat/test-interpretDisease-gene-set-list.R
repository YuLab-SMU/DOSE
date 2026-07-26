test_that("interpretDisease maps a named gene-set list into one canonical object", {
    gene_sets <- list(
        cluster_a = c("1", "2", "9", "10"),
        cluster_b = c("1", "2", "9", "10", "11")
    )

    x <- interpretDisease(
        gene_sets,
        organism = "human",
        ontology = "HDO"
    )
    result_tbl <- as.data.frame(x)
    evidence_tbl <- evidence(x)

    expect_s4_class(x, "doseInterpretResult")
    expect_equal(sort(x@query$query_id), sort(names(gene_sets)))
    expect_equal(sort(x@query$query_name), sort(names(gene_sets)))
    expect_true(all(x@query$input_type == "gene_set_list"))
    expect_equal(nrow(x@sources), 1)
    expect_true(all(result_tbl$query_id %in% names(gene_sets)))
    expect_true(all(evidence_tbl$query_id %in% names(gene_sets)))

    for (query_id in names(gene_sets)) {
        direct <- interpretDisease(
            gene_sets[[query_id]],
            organism = "human",
            ontology = "HDO"
        )
        direct_result <- as.data.frame(direct)
        direct_evidence <- evidence(direct)

        query_result <- result_tbl[result_tbl$query_id == query_id, , drop = FALSE]
        query_evidence <- evidence_tbl[evidence_tbl$query_id == query_id, , drop = FALSE]

        expect_equal(query_result$target_id, direct_result$target_id)
        expect_equal(query_result$target_name, direct_result$target_name)
        expect_equal(query_result$p.adjust, direct_result$p.adjust)
        expect_equal(nrow(query_evidence), nrow(direct_evidence))
        expect_true(all(query_evidence$feature_id %in% gene_sets[[query_id]]))
        expect_true(all(query_evidence$target_id %in% query_result$target_id))
    }
})

test_that("gene-set-list interpretation preserves query metadata when all sets miss", {
    gene_sets <- list(
        cluster_a = c("1", "2", "9", "10"),
        cluster_b = c("1", "2", "9", "10", "11")
    )

    x <- interpretDisease(
        gene_sets,
        organism = "human",
        ontology = "HDO",
        pvalueCutoff = 0,
        qvalueCutoff = 0
    )

    expect_s4_class(x, "doseInterpretResult")
    expect_equal(nrow(as.data.frame(x)), 0)
    expect_equal(nrow(evidence(x)), 0)
    expect_equal(sort(x@query$query_id), sort(names(gene_sets)))
    expect_equal(sort(x@query$query_name), sort(names(gene_sets)))
    expect_true(all(x@query$input_type == "gene_set_list"))
})
