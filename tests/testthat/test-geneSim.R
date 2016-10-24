library(DOSE)

context("geneSim")

test_that("geneSim", {
    res <-  sapply(c("Wang", "Lin", "Jiang", "Resnik", "Rel"), function(mm) geneSim('2524', '3070', measure=mm))
    expect_true(res >=0 && res <=1)
})

