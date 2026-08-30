test_that("GSEA entry points expose a seed argument defaulting to FALSE", {
    expect_equal(formals(DOSE::gseDO)$seed, FALSE)
    expect_equal(formals(DOSE::gseNCG)$seed, FALSE)
    expect_equal(formals(DOSE:::gseDisease)$seed, FALSE)
})

test_that("gseDisease forwards seed to enrichit::gsea_gson", {
    captured <- NULL
    local_mocked_bindings(
        gsea_gson = function(...) {
            captured <<- list(...)
            NULL
        },
        get_anno_data = function(ontology) structure(list(), class = "GSON"),
        .package = "DOSE"
    )

    geneList <- c("7157" = 2, "7272" = 1, "1956" = -1)

    DOSE:::gseDisease(geneList = geneList, ontology = "HDO", seed = 42, verbose = FALSE)
    expect_equal(captured$seed, 42)

    DOSE::gseDO(geneList = geneList, ont = "HDO", seed = 99, verbose = FALSE)
    expect_equal(captured$seed, 99)

    DOSE::gseNCG(geneList = geneList, seed = 7, verbose = FALSE)
    expect_equal(captured$seed, 7)
})
