test_that(
  "Check whether the covariance tests work properly",
  {
    n1 = 100; n2 = 200; pp = 500
    X = matrix(rnorm(n1*pp), n1, pp)
    Y = matrix(rnorm(n2*pp), n2, pp)
    output <- covtest(X,Y)
    expect_equal(names(output), c("method", "stat", "pval"))
    expect_no_error(covtest(X,Y,method='clx'))
    expect_no_error(covtest(X,Y,method='lc'))
    expect_no_error(covtest(X,Y,method='pe.cauchy'))
    expect_no_error(covtest(X,Y,method='pe.comp'))
    expect_no_error(covtest(X,Y,method='pe.fisher'))
  }
)
