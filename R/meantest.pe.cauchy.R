#' Two-sample PE mean test for high-dimensional data via Cauchy combination
#' @description
#' This function implements the two-sample PE covariance test via
#' Cauchy combination.
#' Suppose \eqn{\{\mathbf{X}_1, \ldots, \mathbf{X}_{n_1}\}} are i.i.d.
#' copies of \eqn{\mathbf{X}}, and \eqn{\{\mathbf{Y}_1, \ldots, \mathbf{Y}_{n_2}\}}
#' are i.i.d. copies of \eqn{\mathbf{Y}}.
#' Let \eqn{p_{CQ}} and \eqn{p_{CLX}} denote the \eqn{p}-values associated with
#' the \eqn{l_2}-norm-based covariance test (see \code{\link{meantest.cq}} for details)
#' and the \eqn{l_\infty}-norm-based covariance test
#' (see \code{\link{meantest.clx}} for details), respectively.
#' The PE covariance test via Cauchy combination is defined as
#' \deqn{M_{Cauchy} = \frac{1}{2}\tan((0.5-p_{CQ})\pi) + \frac{1}{2}\tan((0.5-p_{CLX})\pi).}
#' It has been proved that with some regularity conditions, under the null hypothesis
#' \eqn{H_{0m}: \boldsymbol{\mu}_1 = \boldsymbol{\mu}_2,}
#' the two tests are asymptotically independent as \eqn{n_1, n_2, p\rightarrow \infty},
#' and therefore \eqn{M_{Cauchy}} asymptotically converges in distribution to a standard Cauchy distribution.
#' The asymptotic \eqn{p}-value is obtained by
#' \deqn{p\text{-value} = 1-F_{Cauchy}(M_{Cauchy}),}
#' where \eqn{F_{Cauchy}(\cdot)} is the cdf of the standard Cauchy distribution.
#' @import stats
#' @usage
#' meantest.pe.cauchy(dataX,dataY)
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
#' Cai, T. T., Liu, W., and Xia, Y. (2014). Two-sample test of high dimensional
#'  means under dependence. \emph{Journal of the Royal Statistical Society:
#'  Series B: Statistical Methodology}, 76(2):349–372.
#' @examples
#' n1 = 100; n2 = 100; pp = 500
#' set.seed(1)
#' X = matrix(rnorm(n1*pp), nrow=n1, ncol=pp)
#' Y = matrix(rnorm(n2*pp), nrow=n2, ncol=pp)
#' meantest.pe.cauchy(X,Y)

meantest.pe.cauchy <- function(dataX,dataY)
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
  n1=nrow(X);n2=nrow(Y);p=ncol(X)
  output_mean_CQ  = meantest.cq(X,Y)
  output_mean_CLX = meantest.clx(X,Y)
  pval_mean_CQ  = output_mean_CQ$pval
  pval_mean_CLX = output_mean_CLX$pval
  stat_mean_Cauchy =0.5*tan((0.5-pval_mean_CQ)*pi) + 0.5*tan((0.5-pval_mean_CLX)*pi)
  pval_mean_Cauchy = 1-pcauchy(stat_mean_Cauchy)
  return(list(stat = stat_mean_Cauchy,
              pval = pval_mean_Cauchy))
}
