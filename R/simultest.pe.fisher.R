#' Two-sample PE simultaneous test using Fisher's combination
#' @description
#' This function implements the two-sample PE simultaneous test on
#' high-dimensional mean vectors and covariance matrices using Fisher's combination.
#' Suppose \eqn{\{\mathbf{X}_1, \ldots, \mathbf{X}_{n_1}\}} are i.i.d.
#' copies of \eqn{\mathbf{X}}, and \eqn{\{\mathbf{Y}_1, \ldots, \mathbf{Y}_{n_2}\}}
#' are i.i.d. copies of \eqn{\mathbf{Y}}.
#' Let \eqn{M_{PE}} and \eqn{T_{PE}} denote
#' the PE mean test statistic and PE covariance test statistic, respectively.
#' (see \code{\link{meantest.pe.comp}}
#' and \code{\link{covtest.pe.comp}} for details).
#' Let \eqn{p_{m}} and \eqn{p_{c}} denote their respective \eqn{p}-values.
#' The PE simultaneous test statistic via Fisher's combination is defined as
#' \deqn{J_{PE} = -2\log(p_{m})-2\log(p_{c}).}
#' It has been proved that with some regularity conditions, under the null hypothesis
#' \eqn{H_0: \boldsymbol{\mu}_1 = \boldsymbol{\mu}_2 \ \text{ and }
#' \ \mathbf{\Sigma}_1 = \mathbf{\Sigma}_2},
#' the two tests are asymptotically independent as \eqn{n_1, n_2, p\rightarrow \infty},
#' and therefore \eqn{J_{PE}} asymptotically converges in distribution to
#' a \eqn{\chi_4^2} distribution.
#' The asymptotic \eqn{p}-value is obtained by
#' \deqn{p\text{-value} = 1-F_{\chi_4^2}(J_{PE}),}
#' where \eqn{F_{\chi_4^2}(\cdot)} is the cdf of the \eqn{\chi_4^2} distribution.
#' @import stats
#' @usage
#' simultest.pe.fisher(dataX,dataY,delta_mean=NULL,delta_cov=NULL)
#' @param dataX an \eqn{n_1} by \eqn{p} data matrix
#' @param dataY an \eqn{n_2} by \eqn{p} data matrix
#' @param delta_mean a scalar; the thresholding value used in the construction of
#' the PE component for mean test; see \code{\link{meantest.pe.comp}} for details.
#' @param delta_cov a scalar; the thresholding value used in the construction of
#' the PE component for covariance test; see \code{\link{covtest.pe.comp}} for details.
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
#' simultest.pe.fisher(X,Y)

simultest.pe.fisher <- function(dataX,dataY,delta_mean=NULL,delta_cov=NULL)
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
  if(is.null(delta_mean))
  {
    delta_mean = 2*log(log(n1+n2))*log(p)
  }
  if(is.null(delta_cov))
  {
    delta_cov = 4*log(log(n1+n2))*log(p)
  }
  output_mean_pe = meantest.pe.comp(X,Y,delta_mean)
  output_cov_pe  = covtest.pe.comp(X,Y,delta_cov)
  pval_mean_pe   = output_mean_pe$pval
  pval_cov_pe    = output_cov_pe$pval
  stat_simu_pe_fisher = -2*log(pval_mean_pe)-2*log(pval_cov_pe)
  pval_simu_pe_fisher = 1-pchisq(stat_simu_pe_fisher, df=4)
  return(list(stat = stat_simu_pe_fisher,
              pval = pval_simu_pe_fisher))
}



