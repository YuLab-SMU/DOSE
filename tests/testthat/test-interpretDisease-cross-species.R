cross_species_fixture_path <- function(name) {
    testthat::test_path("fixtures", "mouse-models", name)
}

mock_mouse_model_result <- function() {
    doseInterpretResult(
        result = data.frame(
            query_id = "q1",
            target_id = "MGI:3700824",
            target_name = "Gli3<tm1Urt>/Gli3<tm1Urt>",
            target_type = "mouse_model",
            rank = 1L,
            score = 1,
            score_type = "ortholog_support_count",
            score_direction = "higher",
            pvalue = NA_real_,
            p.adjust = NA_real_,
            evidence_count = 2L,
            top_evidence_type = "ortholog",
            source_count = 3L,
            stringsAsFactors = FALSE
        ),
        evidence = data.frame(
            evidence_id = c(
                "MGI:3700824::ortholog::2737",
                "MGI:3700824::phenotype_profile::MP:0002135::11978771"
            ),
            query_id = c("q1", "q1"),
            target_id = c("MGI:3700824", "MGI:3700824"),
            target_type = c("mouse_model", "mouse_model"),
            evidence_type = c("ortholog", "phenotype_profile"),
            direction = c("support", "support"),
            feature_id = c("2737", "MP:0002135"),
            feature_name = c("GLI3", "MP:0002135"),
            component_score = c(1, NA_real_),
            source = c("HOM_MouseHumanSequence", "MGI_GenePheno"),
            source_record_id = c("51812157", "MGI:3700824"),
            evidence_path_id = c(
                "MGI:3700824::ortholog",
                "MGI:3700824::phenotype_profile"
            ),
            derived_from_evidence_id = c(NA_character_, NA_character_),
            reference_id = c(NA_character_, "11978771"),
            note = c(
                "Supports mouse marker Gli3 through the supplied mouse-human ortholog class.",
                "Model phenotype support for Pallister-Hall syndrome (DOID:9248)."
            ),
            stringsAsFactors = FALSE
        ),
        query = data.frame(
            query_id = "q1",
            input_type = "gene",
            organism = "human",
            stringsAsFactors = FALSE
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
            top = 50
        ),
        explanation = list(method = "none", text = character(), metadata = list())
    )
}

test_that("cross-species targets require explicit contract files", {
    expect_error(
        interpretDisease(c("1", "2"), target = "mouse_model"),
        "requires the supplied cross-species contract files"
    )

    expect_error(
        interpretDisease(c("1", "2"), target = "both"),
        "requires the supplied cross-species contract files"
    )
})

test_that("interpretDisease combines disease and mouse-model targets without collapsing target type", {
    observed_gene <- NULL

    local_mocked_bindings(
        enrichDisease = function(gene, organism, ontology, ...) {
            data.frame(
                ID = "DOID:9256",
                Description = "Pallister-Hall syndrome",
                FoldEnrichment = 4.2,
                pvalue = 0.001,
                p.adjust = 0.01,
                Count = 1L,
                geneID = "2737",
                stringsAsFactors = FALSE
            )
        },
        rankMouseModels = function(gene,
                                   disease_model_file,
                                   homology_file,
                                   gene_pheno_file,
                                   top = 50) {
            observed_gene <<- gene
            mock_mouse_model_result()
        },
        .package = "DOSE"
    )

    x <- interpretDisease(
        c("2737", "9999"),
        target = "both",
        disease_model_file = cross_species_fixture_path("MGI_DiseaseMouseModel.gli3.rpt"),
        homology_file = cross_species_fixture_path("HOM_MouseHumanSequence.gli3.rpt"),
        gene_pheno_file = cross_species_fixture_path("MGI_GenePheno.gli3.rpt")
    )

    result_tbl <- as.data.frame(x)
    evidence_tbl <- evidence(x)

    expect_equal(observed_gene, c("2737", "9999"))
    expect_equal(result_tbl$target_type, c("disease", "mouse_model"))
    expect_equal(result_tbl$rank, c(1L, 1L))
    expect_true(all(c("DOID:9256", "MGI:3700824") %in% result_tbl$target_id))
    expect_true(all(c("gene", "ortholog", "phenotype_profile") %in% evidence_tbl$evidence_type))
    expect_true(all(c("HDO", "MGI_DiseaseMouseModel", "HOM_MouseHumanSequence", "MGI_GenePheno") %in% x@sources$source))
})

test_that("ranked gene mouse-model interpretation uses gene names as the cross-species query", {
    observed_gene <- NULL

    local_mocked_bindings(
        rankMouseModels = function(gene,
                                   disease_model_file,
                                   homology_file,
                                   gene_pheno_file,
                                   top = 50) {
            observed_gene <<- gene
            mock_mouse_model_result()
        },
        .package = "DOSE"
    )

    x <- interpretDisease(
        c(`2737` = 2, `9999` = 1),
        target = "mouse_model",
        disease_model_file = cross_species_fixture_path("MGI_DiseaseMouseModel.gli3.rpt"),
        homology_file = cross_species_fixture_path("HOM_MouseHumanSequence.gli3.rpt"),
        gene_pheno_file = cross_species_fixture_path("MGI_GenePheno.gli3.rpt")
    )

    expect_equal(observed_gene, c("2737", "9999"))
    expect_true(all(as.data.frame(x)$target_type == "mouse_model"))
    expect_equal(x@query$input_type[[1]], "ranked_gene")
})

test_that("gene-set-list both preserves query metadata across disease and model slices", {
    local_mocked_bindings(
        enrichDisease = function(gene, organism, ontology, ...) {
            data.frame(
                ID = "DOID:9256",
                Description = "Pallister-Hall syndrome",
                FoldEnrichment = 4.2,
                pvalue = 0.001,
                p.adjust = 0.01,
                Count = 1L,
                geneID = gene[[1]],
                stringsAsFactors = FALSE
            )
        },
        rankMouseModels = function(gene,
                                   disease_model_file,
                                   homology_file,
                                   gene_pheno_file,
                                   top = 50) {
            mock_mouse_model_result()
        },
        .package = "DOSE"
    )

    x <- interpretDisease(
        list(cluster_a = c("2737", "1"), cluster_b = c("2737", "2")),
        target = "both",
        disease_model_file = cross_species_fixture_path("MGI_DiseaseMouseModel.gli3.rpt"),
        homology_file = cross_species_fixture_path("HOM_MouseHumanSequence.gli3.rpt"),
        gene_pheno_file = cross_species_fixture_path("MGI_GenePheno.gli3.rpt")
    )

    result_tbl <- as.data.frame(x)

    expect_equal(sort(unique(result_tbl$query_id)), c("cluster_a", "cluster_b"))
    expect_true(all(c("disease", "mouse_model") %in% result_tbl$target_type))
    expect_equal(sort(x@query$query_name), c("cluster_a", "cluster_b"))
    expect_true(all(x@query$input_type == "gene_set_list"))
})
