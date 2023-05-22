#' Two-sample simultaneous test using Fisher's combination
#' @description
#' This function implements the two-sample simultaneous test on high-dimensional
#' mean vectors and covariance matrices using Fisher's combination.
#' Suppose \eqn{\{\mathbf{X}_1, \ldots, \mathbf{X}_{n_1}\}} are i.i.d.
#' copies of \eqn{\mathbf{X}}, and \eqn{\{\mathbf{Y}_1, \ldots, \mathbf{Y}_{n_2}\}}
#' are i.i.d. copies of \eqn{\mathbf{Y}}.
#' Let \eqn{p_{CQ}} and \eqn{p_{LC}} denote the \eqn{p}-values associated with
#' the \eqn{l_2}-norm-based mean test proposed in Chen and Qin (2010)
#' (see \code{\link{meantest.cq}} for details)
#' and the \eqn{l_2}-norm-based covariance test proposed in Li and Chen (2012)
#' (see \code{\link{covtest.lc}} for details),
#' respectively.
#' The simultaneous test statistic via Fisher's combination is defined as
#' \deqn{J_{n_1, n_2} = -2\log(p_{CQ}) -2\log(p_{LC}).}
#' It has been proved that with some regularity conditions, under the null hypothesis
#' \eqn{H_0: \boldsymbol{\mu}_1 = \boldsymbol{\mu}_2 \ \text{ and }
#' \ \mathbf{\Sigma}_1 = \mathbf{\Sigma}_2},
#' the two tests are asymptotically independent as \eqn{n_1, n_2, p\rightarrow \infty},
#' and therefore \eqn{J_{n_1,n_2}} asymptotically converges in distribution to
#' a \eqn{\chi_4^2} distribution.
#' The asymptotic \eqn{p}-value is obtained by
#' \deqn{p\text{-value} = 1-F_{\chi_4^2}(J_{n_1,n_2}),}
#' where \eqn{F_{\chi_4^2}(\cdot)} is the cdf of the \eqn{\chi_4^2} distribution.
#' @import stats
#' @usage
#' simultest.fisher(dataX,dataY)
#' @param dataX an \eqn{n_1} by \eqn{p} data matrix
#' @param dataY an \eqn{n_2} by \eqn{p} data matrix
#' @return
#' `stat` the value of test statistic
#'
#' `pval` the p-value for the test.
#' @export
#' @references
#' Chen, S. X. and Qin, Y. L. (2010). A two-sample test for high-dimensional data
#' with applications to gene-set testing.
#' \emph{Annals of Statistics}, 38(2):808–835.
#'
#' Li, J. and Chen, S. X. (2012). Two sample tests for high-dimensional
#' covariance matrices. \emph{The Annals of Statistics}, 40(2):908–940.
#'
#' Yu, X., Li, D., Xue, L., and Li, R. (2022). Power-enhanced simultaneous test
#' of high-dimensional mean vectors and covariance matrices with application
#' to gene-set testing. \emph{Journal of the American Statistical Association},
#' (in press):1–14.
#' @examples
#' n1 = 100; n2 = 100; pp = 500
#' set.seed(1)
#' X = matrix(rnorm(n1*pp), nrow=n1, ncol=pp)
#' Y = matrix(rnorm(n2*pp), nrow=n2, ncol=pp)
#' simultest.fisher(X,Y)

simultest.fisher <- function(dataX,dataY)
{
  X=dataX;Y=dataY
  n1=nrow(X);n2=nrow(Y);p=ncol(X);p2=ncol(Y)
  if(p!= p2)
  {
    stop(" The data dimensions ncol(dataX) and ncol(dataY) do not match!")
  }

  if(p <= 30)
  {
    warning(paste0("These methods are designed for high-dimensional data!
                   The data dimension (p=", p, ") is small in the input data,
                   which may results in an inflated Type-I error rate."))
  }
  output_mean_CQ = meantest.cq(X,Y)
  output_cov_LC  = covtest.lc(X,Y)
  pval_mean_CQ   = output_mean_CQ$pval
  pval_cov_LC    = output_cov_LC$pval
  stat_simu_fisher = -2*log(pval_mean_CQ)-2*log(pval_cov_LC)
  pval_simu_fisher = 1-pchisq(stat_simu_fisher, df=4)
  return(list(stat = stat_simu_fisher,
              pval = pval_simu_fisher))
}



