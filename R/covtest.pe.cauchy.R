#' Two-sample PE covariance test for high-dimensional data via Cauchy combination
#' @description
#' This function implements the two-sample PE covariance test via
#' Cauchy combination.
#' Suppose \eqn{\{\mathbf{X}_1, \ldots, \mathbf{X}_{n_1}\}} are i.i.d.
#' copies of \eqn{\mathbf{X}}, and \eqn{\{\mathbf{Y}_1, \ldots, \mathbf{Y}_{n_2}\}}
#' are i.i.d. copies of \eqn{\mathbf{Y}}.
#' Let \eqn{p_{LC}} and \eqn{p_{CLX}} denote the \eqn{p}-values associated with
#' the \eqn{l_2}-norm-based covariance test (see \code{\link{covtest.lc}} for details)
#' and the \eqn{l_\infty}-norm-based covariance test
#' (see \code{\link{covtest.clx}} for details), respectively.
#' The PE covariance test via Cauchy combination is defined as
#' \deqn{T_{Cauchy} = \frac{1}{2}\tan((0.5-p_{LC})\pi) + \frac{1}{2}\tan((0.5-p_{CLX})\pi).}
#' It has been proved that with some regularity conditions, under the null hypothesis
#' \eqn{H_{0c}: \mathbf{\Sigma}_1 = \mathbf{\Sigma}_2,}
#' the two tests are asymptotically independent as \eqn{n_1, n_2, p\rightarrow \infty},
#' and therefore \eqn{T_{Cauchy}} asymptotically converges in distribution to a standard Cauchy distribution.
#' The asymptotic \eqn{p}-value is obtained by
#' \deqn{p\text{-value} = 1-F_{Cauchy}(T_{Cauchy}),}
#' where \eqn{F_{Cauchy}(\cdot)} is the cdf of the standard Cauchy distribution.
#' @import stats
#' @usage
#' covtest.pe.cauchy(dataX,dataY)
#' @param dataX an \eqn{n_1} by \eqn{p} data matrix
#' @param dataY an \eqn{n_2} by \eqn{p} data matrix
#' @return
#' `stat` the value of test statistic
#'
#' `pval` the p-value for the test.
#' @export
#' @references
#' Yu, X., Li, D., and Xue, L. (2022). Fisher’s combined probability test
#' for high-dimensional covariance matrices. \emph{Journal of the American
#' Statistical Association}, (in press):1–14.
#' @examples
#' n1 = 100; n2 = 100; pp = 500
#' set.seed(1)
#' X = matrix(rnorm(n1*pp), nrow=n1, ncol=pp)
#' Y = matrix(rnorm(n2*pp), nrow=n2, ncol=pp)
#' covtest.pe.cauchy(X,Y)

covtest.pe.cauchy <- function(dataX,dataY)
{
  X=dataX; Y=dataY
  n1=nrow(X);n2=nrow(Y);p=ncol(X); p2=ncol(Y)
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
  output_cov_LC  = covtest.lc(X,Y)
  output_cov_CLX = covtest.clx(X,Y)
  pval_cov_LC  = output_cov_LC$pval
  pval_cov_CLX = output_cov_CLX$pval
  stat_cov_Cauchy =0.5*tan((0.5-pval_cov_LC)*pi) + 0.5*tan((0.5-pval_cov_CLX)*pi)
  pval_cov_Cauchy = 1-pcauchy(stat_cov_Cauchy)
  return(list(stat = stat_cov_Cauchy,
              pval = pval_cov_Cauchy))
}
