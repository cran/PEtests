#' Two-sample PE covariance test for high-dimensional data via PE component
#' @description
#' This function implements the two-sample PE covariance test via the
#' construction of the PE component. Let \eqn{T_{LC}/\hat\sigma_{T_{LC}}}
#' denote the \eqn{l_2}-norm-based covariance test statistic
#' (see \code{\link{covtest.lc}} for details).
#' The PE component is constructed by
#' \deqn{J_c=\sqrt{p}\sum_{i=1}^p\sum_{j=1}^p T_{ij}\widehat\xi^{-1/2}_{ij}
#' \mathcal{I}\{ \sqrt{2}T_{ij}\widehat\xi^{-1/2}_{ij} +1 >  \delta_{cov} \}, }
#' where \eqn{\delta_{cov}} is a threshold for the screening procedure,
#' recommended to take the value of \eqn{\delta_{cov}=4\log(\log (n_1+n_2))\log p}.
#' The explicit forms of \eqn{T_{ij}} and \eqn{\widehat\xi_{ij}}
#' can be found in Section 3.2 of Yu et al. (2022).
#' The PE covariance test statistic is defined as
#' \deqn{T_{PE}=T_{LC}/\hat\sigma_{T_{LC}}+J_c.}
#' With some regularity conditions, under the null hypothesis
#' \eqn{H_{0c}: \mathbf{\Sigma}_1 = \mathbf{\Sigma}_2},
#' the test statistic \eqn{T_{PE}} converges in distribution to
#' a standard normal distribution as \eqn{n_1, n_2, p \rightarrow \infty}.
#' The asymptotic \eqn{p}-value is obtained by
#' \deqn{p\text{-value}=1-\Phi(T_{PE}),}
#' where \eqn{\Phi(\cdot)} is the cdf of the standard normal distribution.
#' @import stats
#' @usage
#' covtest.pe.comp(dataX,dataY,delta=NULL)
#' @param dataX an \eqn{n_1} by \eqn{p} data matrix
#' @param dataY an \eqn{n_2} by \eqn{p} data matrix
#' @param delta a scalar; the thresholding value used in the construction of
#' the PE component. If not specified, the function uses a default value
#' \eqn{\delta_{cov}=4\log(\log (n_1+n_2))\log p}.
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
#' covtest.pe.comp(X,Y)

covtest.pe.comp <- function(dataX,dataY,delta=NULL)
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
    delta = 4*log(log(n1+n2))*log(p)
  }

  XTX = t(X) %*% X
  X2TX2 = t(X**2) %*% (X**2)
  YTY = t(Y) %*% Y
  Y2TY2 = t(Y**2) %*% (Y**2)

  # A1ij
  Aij1 = ( XTX^2 - X2TX2)/(n1*(n1-1))
  # B1ij
  Bij1 = ( YTY^2 -Y2TY2)/(n2*(n2-1))
  # C1ij
  Cij1 = (XTX * YTY)/(n1*n2)

  # A2ij & B2ij
  xbar = apply(X, 2, mean)
  ybar = apply(Y, 2, mean)
  xbarmat = matrix(xbar, p, n1, byrow = F)
  ybarmat = matrix(ybar, p, n2, byrow = F)
  Aij2 = 2/(n1*(n1-1)*(n1-2)) *( n1*n1*(XTX)*(xbar %*% t(xbar))-
                                   n1*(t(X)*xbarmat) %*% (X**2)-
                                   n1* t(X**2) %*% (X*t(xbarmat))-
                                   XTX*XTX + 2*X2TX2 )
  Bij2 = 2/(n2*(n2-1)*(n2-2)) *( n2*n2*(YTY)*(ybar %*% t(ybar))-
                                   n2*(t(Y)*ybarmat) %*% (Y**2)-
                                   n2* t(Y**2) %*% (Y*t(ybarmat))-
                                   YTY*YTY + 2*Y2TY2 )
  # C23ij & C32ij
  Cij2 = -((n1*n1*xbar %*% t(xbar) - XTX) * (YTY))/(n1*n2*(n1-1))
  Cij3 = -((n2*n2*ybar %*% t(ybar) - YTY) * (XTX))/(n1*n2*(n2-1))

  # A3ij & B3ij
  x2bar = apply(X^2, 2, mean)
  y2bar = apply(Y^2, 2, mean)
  Aij3 = (n1^4*(xbar %*% t(xbar))^2 +
            -n1*n1* XTX * (xbar %*% t(xbar)) +
            -n1*n1*n1*((xbar*xbar) %*% t(x2bar)) +
            -n1*n1*(xbar %*% t(xbar))*XTX +
            -n1*n1*n1*( x2bar %*% t(xbar*xbar) ) +
            -n1*n1* (xbar %*% t(xbar)) * XTX +
            2*n1*t(X**2) %*% (X*t(xbarmat)) +
            2*n1*(t(X)*xbarmat) %*% (X**2) +
            n1*n1* x2bar %*% t(x2bar) +
            XTX * XTX +
            -n1*n1* XTX * (xbar %*% t(xbar)) +
            2*n1*(t(X)*xbarmat) %*% (X**2) +
            2*n1*(t(X**2) %*% (X*t(xbarmat))) +
            XTX*XTX - 6*X2TX2)/(n1*(n1-1)*(n1-2)*(n1-3))

  Bij3 = (n2^4*(ybar %*% t(ybar))^2 +
            -n2*n2* YTY * (ybar %*% t(ybar)) +
            -n2*n2*n2*((ybar*ybar) %*% t(y2bar)) +
            -n2*n2*(ybar %*% t(ybar))*YTY +
            -n2*n2*n2*( y2bar %*% t(ybar*ybar) ) +
            -n2*n2* (ybar %*% t(ybar)) * YTY +
            2*n2*t(Y**2) %*% (Y*t(ybarmat)) +
            2*n2*(t(Y)*ybarmat) %*% (Y**2) +
            n2*n2* y2bar %*% t(y2bar) +
            YTY * YTY +
            -n2*n2* YTY * (ybar %*% t(ybar)) +
            2*n2*(t(Y)*ybarmat) %*% (Y**2) +
            2*n2*(t(Y**2) %*% (Y*t(ybarmat))) +
            YTY*YTY - 6*Y2TY2)/(n2*(n2-1)*(n2-2)*(n2-3))

  # C4ij
  Cij4 = (n1*n1*xbar %*% t(xbar) - XTX) * (n2*n2*ybar %*% t(ybar) - YTY)/(n1*n2*(n1-1)*(n2-1))

  Anij=Aij1-Aij2+Aij3; Bnij=Bij1-Bij2+Bij3
  Cnij=Cij1+Cij2+Cij3+Cij4
  Tij = Anij + Bnij - 2*Cnij

  covx = cov(X)
  covy = cov(Y)
  varx = diag(covx)
  vary = diag(covy)
  varxmat = varx %*% t(varx)
  varymat = vary %*% t(vary)
  varmat = 2*((1/n1)*covx^2 + (1/n2)*covy^2+ (1/n1)*varxmat + (1/n2)*varymat)^2
  sdmat = sqrt(varmat)

  stat_normalized = Tij/sdmat
  delta = 4*log(log(n1+n2))*log(p)
  idx = sqrt(2)*stat_normalized + 1 > delta
  J0 = sqrt(p)*sum(stat_normalized[idx])

  stat_cov_LC = covtest.lc(X,Y)$stat
  stat_cov_PE = stat_cov_LC + J0
  pval_cov_PE = 1-pnorm(stat_cov_PE)

  return(list(stat = stat_cov_PE,
              pval = pval_cov_PE))
}
