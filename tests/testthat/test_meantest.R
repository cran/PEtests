test_that(
  "Check whether the mean tests work properly",
  {
    n1 = 100; n2 = 200; pp = 500
    X = matrix(rnorm(n1*pp), n1, pp)
    Y = matrix(rnorm(n2*pp), n2, pp)
    output <- meantest(X,Y)
    expect_equal(names(output), c("method", "stat", "pval"))
    expect_no_error(meantest(X,Y,method='clx'))
    expect_no_error(meantest(X,Y,method='cq'))
    expect_no_error(meantest(X,Y,method='pe.cauchy'))
    expect_no_error(meantest(X,Y,method='pe.comp'))
    expect_no_error(meantest(X,Y,method='pe.fisher'))
  }
)
