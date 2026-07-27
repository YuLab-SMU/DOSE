fixture_path <- function(name) {
    testthat::test_path("fixtures", "mouse-models", name)
}

test_that("rankMouseModels returns canonical mouse-model targets from a real report slice", {
    x <- rankMouseModels(
        gene = "2737",
        disease_model_file = fixture_path("MGI_DiseaseMouseModel.gli3.rpt"),
        homology_file = fixture_path("HOM_MouseHumanSequence.gli3.rpt"),
        gene_pheno_file = fixture_path("MGI_GenePheno.gli3.rpt"),
        top = 10
    )

    result_tbl <- as.data.frame(x)
    evidence_tbl <- evidence(x)

    expect_s4_class(x, "doseInterpretResult")
    expect_equal(nrow(result_tbl), 2)
    expect_true(all(result_tbl$target_type == "mouse_model"))
    expect_equal(
        sort(result_tbl$target_id),
        sort(c("MGI:2166944", "MGI:3700824"))
    )
    expect_equal(
        sort(result_tbl$target_name),
        sort(c("Gli3<Xt-J>/Gli3<Xt-J>", "Gli3<tm1Urt>/Gli3<tm1Urt>"))
    )
    expect_true(all(result_tbl$score == 1))
    expect_true(all(result_tbl$score_type == "ortholog_support_count"))

    expect_true(all(evidence_tbl$target_id %in% result_tbl$target_id))
    expect_true(all(evidence_tbl$target_type == "mouse_model"))
    expect_equal(sort(unique(evidence_tbl$evidence_type)), c("ortholog", "phenotype_profile"))
    expect_equal(nrow(evidence(x, type = "ortholog")), 2)
    expect_equal(nrow(evidence(x, type = "phenotype_profile")), 4)
    expect_true(all(evidence(x, type = "ortholog")$feature_id == "2737"))
    expect_true(all(evidence(x, type = "ortholog")$feature_name == "GLI3"))
    expect_identical(anyDuplicated(evidence_tbl$evidence_id), 0L)
})

test_that("rankMouseModels keeps phenotype profiles as evidence rather than targets", {
    x <- rankMouseModels(
        gene = "2737",
        disease_model_file = fixture_path("MGI_DiseaseMouseModel.gli3.rpt"),
        homology_file = fixture_path("HOM_MouseHumanSequence.gli3.rpt"),
        gene_pheno_file = fixture_path("MGI_GenePheno.gli3.rpt")
    )

    result_tbl <- as.data.frame(x)
    phenotype_tbl <- evidence(x, type = "phenotype_profile")

    expect_false(any(grepl("^MP:", result_tbl$target_id)))
    expect_true(all(grepl("^MP:", phenotype_tbl$feature_id)))
    expect_true(all(phenotype_tbl$direction == "support"))
})

test_that("rankMouseModels returns an empty canonical object when no ortholog support is present", {
    x <- rankMouseModels(
        gene = "999999",
        disease_model_file = fixture_path("MGI_DiseaseMouseModel.gli3.rpt"),
        homology_file = fixture_path("HOM_MouseHumanSequence.gli3.rpt"),
        gene_pheno_file = fixture_path("MGI_GenePheno.gli3.rpt")
    )

    expect_s4_class(x, "doseInterpretResult")
    expect_equal(nrow(as.data.frame(x)), 0)
    expect_equal(nrow(evidence(x)), 0)
    expect_equal(x@parameters$target, "mouse_model")
})
