#' Power-Enhanced (PE) Tests for High-Dimensional Data
#'
#' The package implements several two-sample power-enhanced mean tests,
#' covariance tests, and simultaneous tests on mean vectors and covariance matrices
#' for high-dimensional data.
#'
#' There are three main functions:\cr
#' \code{\link{covtest}}\cr
#' \code{\link{meantest}}\cr
#' \code{\link{simultest}}
#'
#' @name PEtests-package
#' @docType package
#' @references
#' Chen, S. X. and Qin, Y. L. (2010). A two-sample test for high-dimensional data
#' with applications to gene-set testing.
#' \emph{Annals of Statistics}, 38(2):808–835.
#' \doi{10.1214/09-AOS716}
#'
#' Cai, T. T., Liu, W., and Xia, Y. (2013). Two-sample covariance matrix testing
#' and support recovery in high-dimensional and sparse settings.
#' \emph{Journal of the American Statistical Association}, 108(501):265–277.
#' \doi{10.1080/01621459.2012.758041}
#'
#' Cai, T. T., Liu, W., and Xia, Y. (2014). Two-sample test of high dimensional
#' means under dependence. \emph{Journal of the Royal Statistical Society:
#' Series B: Statistical Methodology}, 76(2):349–372.
#' \doi{10.1111/rssb.12034}
#'
#' Li, J. and Chen, S. X. (2012). Two sample tests for high-dimensional
#' covariance matrices. \emph{The Annals of Statistics}, 40(2):908–940.
#' \doi{10.1214/12-AOS993}
#'
#' Yu, X., Li, D., and Xue, L. (2022). Fisher’s combined probability test
#' for high-dimensional covariance matrices. \emph{Journal of the American
#' Statistical Association}, (in press):1–14.
#' \doi{10.1080/01621459.2022.2126781}
#'
#'
#' Yu, X., Li, D., Xue, L., and Li, R. (2022). Power-enhanced simultaneous test
#' of high-dimensional mean vectors and covariance matrices with application
#' to gene-set testing. \emph{Journal of the American Statistical Association},
#' (in press):1–14.
#' \doi{10.1080/01621459.2022.2061354}
#'
#' @examples
#'
#' n1 = 100; n2 = 100; pp = 500
#' set.seed(1)
#' X = matrix(rnorm(n1*pp), nrow=n1, ncol=pp)
#' Y = matrix(rnorm(n2*pp), nrow=n2, ncol=pp)
#' covtest(X, Y)
#' meantest(X, Y)
#' simultest(X, Y)
#'
NULL



