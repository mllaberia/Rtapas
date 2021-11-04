test_that("Association matrix works", {
  expect_equal(ass_m(cbind(c('h1','h1','h2'), c('s1','s2','s2'))),
               matrix(c(1,0,1,1), ncol = 2, nrow = 2,
                      dimnames = list(c('h1','h2'), c('s1','s2'))))
})
