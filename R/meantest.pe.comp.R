#' Two-sample PE mean test for high-dimensional data via PE component
#' @description
#' This function implements the two-sample PE mean via the
#' construction of the PE component. Let \eqn{M_{CQ}/\hat\sigma_{M_{CQ}}}
#' denote the \eqn{l_2}-norm-based mean test statistic
#' (see \code{\link{meantest.cq}} for details).
#' The PE component is constructed by
#' \deqn{J_m = \sqrt{p}\sum_{i=1}^p M_i\widehat\nu^{-1/2}_i
#' \mathcal{I}\{ \sqrt{2}M_i\widehat\nu^{-1/2}_i + 1 >  \delta_{mean} \}, }
#' where \eqn{\delta_{mean}} is a threshold for the screening procedure,
#' recommended to take the value of \eqn{\delta_{mean}=2\log(\log (n_1+n_2))\log p}.
#' The explicit forms of \eqn{M_{i}} and \eqn{\widehat\nu_{j}}
#' can be found in Section 3.1 of Yu et al. (2022).
#' The PE covariance test statistic is defined as
#' \deqn{M_{PE}=M_{CQ}/\hat\sigma_{M_{CQ}}+J_m.}
#' With some regularity conditions, under the null hypothesis
#' \eqn{H_{0m}:  \boldsymbol{\mu}_1 = \boldsymbol{\mu}_2},
#' the test statistic \eqn{M_{PE}} converges in distribution to
#' a standard normal distribution as \eqn{n_1, n_2, p \rightarrow \infty}.
#' The asymptotic  \eqn{p}-value is obtained by
#' \deqn{p\text{-value}= 1-\Phi(M_{PE}),}
#' where \eqn{\Phi(\cdot)} is the cdf of the standard normal distribution.
#' @import stats
#' @usage
#' meantest.pe.comp(dataX,dataY,delta=NULL)
#' @param dataX an \eqn{n_1} by \eqn{p} data matrix
#' @param dataY an \eqn{n_2} by \eqn{p} data matrix
#' @param delta a scalar; the thresholding value used in the construction of
#' the PE component. If not specified, the function uses a default value
#' \eqn{\delta_{mean}=2\log(\log (n_1+n_2))\log p}.
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
#' meantest.pe.comp(X,Y)

meantest.pe.comp <- function(dataX,dataY,delta=NULL)
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

  if(is.null(delta))
  {
    delta = 2*log(log(n1+n2))*log(p)
  }
  xbar = apply(X, 2, mean)
  x2bar = apply(X**2, 2, mean)
  ybar = apply(Y, 2, mean)
  y2bar = apply(Y**2, 2, mean)
  Ak = n1*n1*xbar*xbar - n1*x2bar
  Bk = n2*n2*ybar*ybar - n2*y2bar
  Ck = n1*n2*xbar*ybar
  Tk = Ak/(n1*(n1-1)) + Bk/(n2*(n2-1)) - 2*Ck/(n1*n2)

  sig1 = diag(cov(X))
  sig2 = diag(cov(Y))
  sigk = 2*sig1*sig1/(n1*(n1-1)) + 2*sig2*sig2/(n2*(n2-1)) + 4*sig1*sig2/(n1*n2)
  sdk = sqrt(sigk)

  stat_normalized = Tk/sdk
  idx = sqrt(2)*stat_normalized + 1 > delta
  J0 = sqrt(p)*sum(stat_normalized[idx])

  stat_mean_CQ = meantest.cq(X,Y)$stat
  stat_mean_PE = stat_mean_CQ + J0
  pval_mean_PE = 1-pnorm(stat_mean_PE)

  return(list(stat = stat_mean_PE,
              pval = pval_mean_PE))
}


