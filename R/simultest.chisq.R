#' Two-sample simultaneous test using chi-squared approximation
#' @description
#' This function implements the two-sample simultaneous test on high-dimensional
#' mean vectors and covariance matrices using chi-squared approximation.
#' Suppose \eqn{\{\mathbf{X}_1, \ldots, \mathbf{X}_{n_1}\}} are i.i.d.
#' copies of \eqn{\mathbf{X}}, and \eqn{\{\mathbf{Y}_1, \ldots, \mathbf{Y}_{n_2}\}}
#' are i.i.d. copies of \eqn{\mathbf{Y}}.
#' Let \eqn{M_{CQ}/\hat\sigma_{M_{CQ}}} denote
#' the \eqn{l_2}-norm-based mean test statistic proposed in Chen and Qin (2010)
#' (see \code{\link{meantest.cq}} for details),
#' and let \eqn{T_{LC}/\hat\sigma_{T_{LC}}}
#' denote the \eqn{l_2}-norm-based covariance test statistic
#' proposed in Li and Chen (2012) (see \code{\link{covtest.lc}} for details).
#' The simultaneous test statistic via chi-squared approximation is defined as
#' \deqn{S_{n_1, n_2} = M_{CQ}^2/\hat\sigma^2_{M_{CQ}} + T_{LC}^2/\hat\sigma^2_{T_{LC}}.}
#' It has been proved that with some regularity conditions, under the null hypothesis
#' \eqn{H_0: \boldsymbol{\mu}_1 = \boldsymbol{\mu}_2 \ \text{ and }
#' \ \mathbf{\Sigma}_1 = \mathbf{\Sigma}_2},
#' the two tests are asymptotically independent as \eqn{n_1, n_2, p\rightarrow \infty},
#' and therefore \eqn{S_{n_1,n_2}} asymptotically converges in distribution to
#' a \eqn{\chi_2^2} distribution.
#' The asymptotic \eqn{p}-value is obtained by
#' \deqn{p\text{-value} = 1-F_{\chi_2^2}(S_{n_1,n_2}),}
#' where \eqn{F_{\chi_2^2}(\cdot)} is the cdf of the \eqn{\chi_2^2} distribution.
#' @import stats
#' @usage
#' simultest.chisq(dataX,dataY)
#' @param dataX n1 by p data matrix
#' @param dataY n2 by p data matrix
#' @return
#' `stat` the value of test statistic
#'
#' `pval` the p-value for the test.
#' @export
#' @references
#' Yu, X., Li, D., Xue, L., and Li, R. (2022). Power-enhanced simultaneous test
#' of high-dimensional mean vectors and covariance matrices with application
#' to gene-set testing. \emph{Journal of the American Statistical Association},
#' (in press):1–14.
#' @examples
#' n1 = 100; n2 = 100; pp = 500
#' set.seed(1)
#' X = matrix(rnorm(n1*pp), nrow=n1, ncol=pp)
#' Y = matrix(rnorm(n2*pp), nrow=n2, ncol=pp)
#' simultest.chisq(X,Y)

simultest.chisq <- function(dataX,dataY)
{
  X=dataX; Y=dataY
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
  stat_mean_CQ   = output_mean_CQ$stat
  stat_cov_LC    = output_cov_LC$stat
  stat_simu_chisq = stat_mean_CQ^2 + stat_cov_LC^2
  pval_simu_chisq = 1-pchisq(stat_simu_chisq, df=2)
  return(list(stat = stat_simu_chisq,
              pval = pval_simu_chisq))
}



