library(DOSE)

context("doSim")

test_that("doSim", {
    res <- sapply(c("Wang", "Lin", "Jiang", "Resnik", "Rel"), function(mm) doSim("DOID:1002", "DOID:10003", measure=mm))
    expect_true(res >= 0 && res <= 1)
})
