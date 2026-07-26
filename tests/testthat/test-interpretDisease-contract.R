test_that("doseInterpretResult can be constructed and inspected offline", {
    result <- data.frame(
        query_id = "q1",
        target_id = "DOID:1",
        target_name = "Example disease",
        target_type = "disease",
        rank = 1L,
        score = 3.5,
        score_type = "ora_score",
        score_direction = "higher",
        pvalue = 0.01,
        p.adjust = 0.02,
        evidence_count = 1L,
        top_evidence_type = "gene",
        source_count = 1L,
        stringsAsFactors = FALSE
    )

    evidence_rows <- data.frame(
        evidence_id = "ev1",
        query_id = "q1",
        target_id = "DOID:1",
        target_type = "disease",
        evidence_type = "gene",
        direction = "support",
        feature_id = "1",
        feature_name = "A1BG",
        component_score = 1,
        source = "HDO",
        source_record_id = "rec1",
        evidence_path_id = "path1",
        derived_from_evidence_id = NA_character_,
        reference_id = NA_character_,
        note = NA_character_,
        stringsAsFactors = FALSE
    )

    query <- data.frame(
        query_id = "q1",
        input_type = "gene",
        organism = "human",
        stringsAsFactors = FALSE
    )

    sources <- data.frame(
        source = "HDO",
        version = "unknown",
        checksum = "na",
        stringsAsFactors = FALSE
    )

    x <- doseInterpretResult(
        result = result,
        evidence = evidence_rows,
        query = query,
        sources = sources,
        parameters = list(ontology = "HDO"),
        explanation = list()
    )

    expect_s4_class(x, "doseInterpretResult")
    expect_equal(nrow(as.data.frame(x)), 1)
    expect_equal(nrow(evidence(x)), 1)
    expect_equal(nrow(evidence(x, target_id = "DOID:1", type = "gene")), 1)

    summary_out <- summary(x)
    expect_equal(summary_out$query_count, 1)
    expect_equal(summary_out$target_count, 1)
    expect_equal(summary_out$evidence_count, 1)

    expect_output(show(x), "doseInterpretResult")
})

test_that("doseInterpretResult rejects broken result or evidence links", {
    result <- data.frame(
        query_id = c("q1", "q1"),
        target_id = c("DOID:1", "DOID:1"),
        target_name = c("Disease", "Disease"),
        target_type = c("disease", "disease"),
        rank = c(1L, 1L),
        score = c(1, 1),
        score_type = c("ora_score", "ora_score"),
        score_direction = c("higher", "higher"),
        pvalue = c(0.1, 0.1),
        p.adjust = c(0.2, 0.2),
        evidence_count = c(1L, 1L),
        top_evidence_type = c("gene", "gene"),
        source_count = c(1L, 1L),
        stringsAsFactors = FALSE
    )

    expect_error(
        doseInterpretResult(result = result),
        "uniquely identified"
    )

    result <- result[1, , drop = FALSE]
    evidence_rows <- data.frame(
        evidence_id = "ev1",
        query_id = "q1",
        target_id = "DOID:2",
        target_type = "disease",
        evidence_type = "gene",
        direction = "support",
        feature_id = "1",
        feature_name = "A1BG",
        component_score = 1,
        source = "HDO",
        source_record_id = "rec1",
        evidence_path_id = "path1",
        derived_from_evidence_id = NA_character_,
        reference_id = NA_character_,
        note = NA_character_,
        stringsAsFactors = FALSE
    )

    expect_error(
        doseInterpretResult(result = result, evidence = evidence_rows),
        "link to a valid result row"
    )
})

test_that("interpretDisease reserves gated paths and validates inputs early", {
    expect_error(
        interpretDisease(c("HP:0000118", "HP:0001250")),
        "phenotype input.*not yet enabled"
    )

    expect_error(
        interpretDisease(c("1", "2"), target = "both"),
        "target = \"both\".*not yet enabled"
    )

    expect_error(
        interpretDisease(c("1", "2"), organism = "mouse"),
        "Cross-species disease interpretation from mouse input is not yet enabled"
    )

    expect_error(
        interpretDisease(c("1", "2"), organism = "mouse", ontology = "HDO"),
        "incompatible.*explicit ortholog conversion"
    )

    expect_error(
        interpretDisease(c("1", "TP53")),
        "keytype = 'ENTREZID'"
    )
})
