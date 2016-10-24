library(DOSE)

context("gene2DO")

test_that("gene2DO", {
    
    expect_equal(DOSE:::gene2DO('3'),
                 DOSE:::gene2DO(3))

})


