test_that("Specifying baseline values works", {
  est1 <- multe(score ~ treatment | school, data = star)$est
  est2 <- multe(score ~ treatment | school, data = star, base_val = "aide")$est
  expect_equal(est1[1,],est2[1,]*(-1))
})

test_that("Both numerical and string treatment variable works", {
  est1 <- multe(score ~ with(
    star, as.integer(factor(treatment, levels = unique(treatment)))
    ) | school, data = star)
  est2 <- multe(score ~ treatment | school, data = star)
  expect_equal(est1$est,est2$est)
  expect_equal(est1$se_po,est2$se_po)
  expect_equal(est1$se_or,est2$se_or)
})

test_that("Both numerical and string control variable works", {
  est1 <- multe(score ~ treatment | as.character(school), data = star)
  est2 <- multe(score ~ treatment | school, data = star)
  expect_equal(est1$est,est2$est)
  expect_equal(est1$se_po,est2$se_po)
  expect_equal(est1$se_or,est2$se_or)
})
