test_that("ontology determines the annotation species", {
    expect_equal(
        DOSE:::.resolve_ontology_organism("HDO"),
        list(
            ontology = "HDO",
            organism = "hsa",
            species = "Homo sapiens",
            keytype = "ENTREZID"
        )
    )
    expect_equal(DOSE:::.resolve_ontology_organism("MPO")$organism, "mm")
    expect_equal(DOSE:::.resolve_ontology_organism("MPO", "mouse")$organism, "mm")
    expect_equal(DOSE:::.resolve_ontology_organism("MPO", "mmu")$organism, "mm")
})

test_that("incompatible ontology and organism fail before analysis", {
    expect_error(
        enrichDO("1", ont = "MPO", organism = "hsa"),
        "incompatible.*explicit ortholog conversion"
    )
    expect_error(
        gseDO(c(`1` = 2, `2` = 1), ont = "MPO", organism = "human"),
        "incompatible.*explicit ortholog conversion"
    )
    expect_error(
        gene2DO("1", ont = "HPO", organism = "mouse"),
        "incompatible.*explicit ortholog conversion"
    )
})

test_that("gene inputs enforce the documented Entrez key type", {
    expect_error(enrichDO("TP53"), "keytype = 'ENTREZID'")
    expect_error(geneSim(c("1", "TP53")), "keytype = 'ENTREZID'")
    expect_error(
        gseDO(c(TP53 = 2, BRCA1 = 1)),
        "keytype = 'ENTREZID'"
    )
    expect_error(
        gseDO(c(`1` = 2, `1` = 1)),
        "unique Entrez gene ID names"
    )
})

test_that("ontology gson metadata records the annotation species", {
    cache <- new.env(parent = emptyenv())
    local_mocked_bindings(
        get_dose_env = function() cache,
        get_ont2allgene = function(...) list(`MP:0000001` = "1"),
        .package = "DOSE"
    )
    local_mocked_bindings(
        get_onto_data = function(...) {
            data.frame(id = "MP:0000001", name = "phenotype")
        },
        .package = "GOSemSim"
    )
    local_mocked_bindings(
        gson = function(...) list(...),
        .package = "gson"
    )

    result <- DOSE:::get_dose_data("MPO")
    expect_equal(result$species, "Mus musculus")
    expect_equal(result$keytype, "ENTREZID")
})

test_that("cluster similarity uses the requested ontology for gene lookup", {
    observed <- character()
    local_mocked_bindings(
        gene2DO = function(gene, organism, ont) {
            observed <<- c(observed, ont)
            "MP:0000001"
        },
        doseSim = function(...) matrix(1, nrow = 1, ncol = 1),
        combineScores = function(...) 1,
        .package = "DOSE"
    )

    expect_equal(clusterSim("1", "2", ont = "MPO"), 1)
    expect_equal(observed, c("MPO", "MPO"))
})

test_that("multi-cluster similarity uses the requested ontology for gene lookup", {
    observed <- character()
    local_mocked_bindings(
        gene2DO = function(gene, organism, ont) {
            observed <<- c(observed, ont)
            "MP:0000001"
        },
        doseSim = function(...) matrix(1, nrow = 1, ncol = 1),
        combineScores = function(...) 1,
        .package = "DOSE"
    )

    result <- mclusterSim(list(a = "1", b = "2"), ont = "MPO")
    expect_equal(dim(result), c(2, 2))
    expect_equal(observed, c("MPO", "MPO"))
})

test_that("NCG data uses the shared package cache", {
    cache <- new.env(parent = emptyenv())
    sentinel <- structure(list(), class = "cached_ncg")
    assign(".NCG_DOSE_GSON", sentinel, envir = cache)
    local_mocked_bindings(
        get_dose_env = function() cache,
        .package = "DOSE"
    )

    expect_identical(DOSE:::get_NCG_data(), sentinel)
})
