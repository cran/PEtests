test_that(
  "Check whether the simultaneous tests work properly",
  {
    n1 = 100; n2 = 200; pp = 500
    X = matrix(rnorm(n1*pp), n1, pp)
    Y = matrix(rnorm(n2*pp), n2, pp)
    output <- simultest(X,Y)
    expect_equal(names(output), c("method", "stat", "pval"))
    expect_no_error(simultest(X,Y,method='cauchy'))
    expect_no_error(simultest(X,Y,method='chisq'))
    expect_no_error(simultest(X,Y,method='fisher'))
    expect_no_error(simultest(X,Y,method='pe.cauchy'))
    expect_no_error(simultest(X,Y,method='pe.chisq'))
    expect_no_error(simultest(X,Y,method='pe.fisher'))
  }
)
