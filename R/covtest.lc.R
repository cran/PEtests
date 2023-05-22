#' Two-sample high-dimensional covariance test (Li and Chen, 2012)
#' @description
#' This function implements the two-sample \eqn{l_2}-norm-based high-dimensional covariance test
#' proposed by Li and Chen (2012).
#' Suppose \eqn{\{\mathbf{X}_1, \ldots, \mathbf{X}_{n_1}\}} are i.i.d.
#' copies of \eqn{\mathbf{X}}, and \eqn{\{\mathbf{Y}_1, \ldots, \mathbf{Y}_{n_2}\}}
#' are i.i.d. copies of \eqn{\mathbf{Y}}. The test statistic \eqn{T_{LC}} is
#' defined as
#' \deqn{T_{LC} = A_{n_1}+B_{n_2}-2C_{n_1,n_2},}
#' where \eqn{A_{n_1}}, \eqn{B_{n_2}}, and \eqn{C_{n_1,n_2}} are unbiased estimators for
#' \eqn{\mathrm{tr}(\mathbf{\Sigma}^2_1)}, \eqn{\mathrm{tr}(\mathbf{\Sigma}^2_2)},
#' and  \eqn{\mathrm{tr}(\mathbf{\Sigma}_1\mathbf{\Sigma}_2)}, respectively.
#' Under the null hypothesis \eqn{H_{0c}: \mathbf{\Sigma}_1 = \mathbf{\Sigma}_2},
#' the leading variance of \eqn{T_{LC}} is
#'  \eqn{\sigma^2_{T_{LC}} = 4(\frac{1}{n_1}+\frac{1}{n_2})^2 \rm{tr}^2(\mathbf{\Sigma}^2)},
#'  which can be consistently estimated by \eqn{\hat\sigma^2_{LC}}.
#'  The explicit formulas of \eqn{A_{n_1}}, \eqn{B_{n_2}}, \eqn{C_{n_1,n_2}}
#' and \eqn{\hat\sigma^2_{T_{LC}}} can be found in
#' Equations (2.1), (2.2) and Theorem 1 of Li and Chen (2012).
#' With some regularity conditions, under the null hypothesis \eqn{H_{0c}: \mathbf{\Sigma}_1 = \mathbf{\Sigma}_2},
#' the test statistic \eqn{T_{LC}} converges in distribution to a standard normal distribution
#' as \eqn{n_1, n_2, p \rightarrow \infty}.
#' The asymptotic  \eqn{p}-value is obtained by
#' \deqn{p_{LC} = 1-\Phi(T_{LC}/\hat\sigma_{T_{LC}}),}
#' where \eqn{\Phi(\cdot)} is the cdf of the standard normal distribution.
#' @import stats
#' @usage
#' covtest.lc(dataX,dataY)
#' @param dataX an \eqn{n_1} by \eqn{p} data matrix
#' @param dataY an \eqn{n_2} by \eqn{p} data matrix
#' @return
#' `stat` the value of test statistic
#'
#' `pval` the p-value for the test.
#' @export
#' @references
#' Li, J. and Chen, S. X. (2012). Two sample tests for high-dimensional
#' covariance matrices. \emph{The Annals of Statistics}, 40(2):908–940.
#' @examples
#' n1 = 100; n2 = 100; pp = 500
#' set.seed(1)
#' X = matrix(rnorm(n1*pp), nrow=n1, ncol=pp)
#' Y = matrix(rnorm(n2*pp), nrow=n2, ncol=pp)
#' covtest.lc(X,Y)

covtest.lc <- function(dataX,dataY)
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
  XXT=X%*%t(X);YYT=Y%*%t(Y);XYT=X%*%t(Y);
  A.1=(sum((XXT)^2)-tr((XXT)^2))/(n1*(n1-1));
  B.1=(sum((YYT)^2)-tr((YYT)^2))/(n2*(n2-1));
  A.2=A_2(X);A.3=A_3(X);B.2=A_2(Y);B.3=A_3(Y)
  An1=A.1-A.2+A.3;Bn2=B.1-B.2+B.3;
  C.1=(sum((XYT)^2))/(n1*n2);
  C.2=C_23(X,Y);C.3=C_23(Y,X);C.4=C_4(X,Y)
  Cn=C.1+C.2+C.3+C.4
  T.value=An1+Bn2-2*Cn;
  sd.est=(2*An1)/n2+(2*Bn2)/n1;
  stat_cov_LC = T.value/sd.est;
  pval_cov_LC = 1-pnorm(stat_cov_LC)
  return(list(stat = stat_cov_LC,
              pval = pval_cov_LC))
}

tr <- function(A)
{
  return(sum(diag(A)))
}

center <-function(X)
{
  n1 = dim(X)[1]
  p = dim(X)[2]
  X.c <- X-matrix(rep(apply(X,2,mean),n1),nrow=n1,ncol=p,byrow=T)
  return(X.c)
}

A_2 <- function(X)
{
  n=nrow(X);p=ncol(X)
  XXT = X %*% t(X)
  xbar = apply(X, 2, mean)
  Xvv = matrix(diag(XXT), n, n, byrow = T)
  Xv_sumk = matrix(X %*% xbar * n, n, n, byrow = T)
  mat1 = XXT
  mat2 = Xv_sumk - Xvv - XXT
  est = sum(mat1*mat2)-sum(diag(mat1*mat2))
  ans = 2*est/(n*(n-1)*(n-2))
  return(ans)
}

A_3 <-function(X)
{
  n=nrow(X);p=ncol(X); XXT=X%*%t(X)
  xbar = apply(X, 2, mean)
  X_all = matrix(sum(XXT) - sum(diag(XXT)), n, n)
  X_uvvec = -2*X%*%xbar*n + 2* diag(XXT)
  Xuu = matrix(X_uvvec, n, n, byrow = F)
  Xvv = matrix(X_uvvec, n, n, byrow = T)
  mat1 = XXT
  mat2 = X_all + Xuu + Xvv + 2* XXT
  est = sum(mat1*mat2)-sum(diag(mat1*mat2))
  ans = est/(n*(n-1)*(n-2)*(n-3))
  return(ans)
}

C_23 <-function(X,Y)
{
  n1=nrow(X);n2=nrow(Y);p=ncol(X)
  xbar = apply(X, 2, mean)
  XYT = X%*%t(Y)
  Yv_sumk = matrix(Y %*% xbar * n1, n1, n2, byrow = T)
  mat1 = XYT
  mat2 = Yv_sumk - XYT
  est = sum(mat1*mat2)
  ans = -est/(n1*n2*(n1-1))
  return(ans)
}

C_4<-function(X,Y)
{
  n1=nrow(X);n2=nrow(Y);p=ncol(X);
  xbar = apply(X, 2, mean); ybar = apply(Y, 2, mean);
  XYT=X%*%t(Y)
  XY_sumkl = matrix(sum(XYT), n1, n2)
  Xk_suml = matrix(X %*% ybar* n2, n1, n2, byrow = F)
  Yl_sumk = matrix(Y %*% xbar* n1, n1, n2, byrow = T)
  mat1 = XYT
  mat2 = XY_sumkl - Xk_suml - Yl_sumk + XYT
  est = sum(mat1*mat2)
  ans = est/(n1*n2*(n1-1)*(n2-1))
  return(ans)
}
