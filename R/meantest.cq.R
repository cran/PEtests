#' Two-sample high-dimensional mean test (Chen and Qin, 2010)
#' @description This function implements the two-sample \eqn{l_2}-norm-based high-dimensional
#' mean test proposed by Chen and Qin (2010).
#' Suppose \eqn{\{\mathbf{X}_1, \ldots, \mathbf{X}_{n_1}\}} are i.i.d.
#' copies of \eqn{\mathbf{X}}, and \eqn{\{\mathbf{Y}_1, \ldots, \mathbf{Y}_{n_2}\}}
#' are i.i.d. copies of \eqn{\mathbf{Y}}.
#' The test statistic \eqn{M_{CQ}} is defined as
#' \deqn{M_{CQ} = \frac{1}{n_1(n_1-1)}\sum_{u\neq v}^{n_1} \mathbf{X}_{u}'\mathbf{X}_{v}
#' +\frac{1}{n_2(n_2-1)}\sum_{u\neq v}^{n_2} \mathbf{Y}_{u}'\mathbf{Y}_{v}
#' -\frac{2}{n_1n_2}\sum_u^{n_1}\sum_v^{n_2} \mathbf{X}_{u}'\mathbf{Y}_{v}.}
#' Under the null hypothesis \eqn{H_{0m}: \boldsymbol{\mu}_1 = \boldsymbol{\mu}_2},
#' the leading variance of \eqn{M_{CQ}} is
#' \eqn{\sigma^2_{M_{CQ}}=\frac{2}{n_1(n_1-1)}\text{tr}(\mathbf{\Sigma}_1^2)+
#' \frac{2}{n_2(n_2-1)}\text{tr}(\mathbf{\Sigma}_2^2)+
#' \frac{4}{n_1n_2}\text{tr}(\mathbf{\Sigma}_1\mathbf{\Sigma}_2)},
#' which can be consistently estimated by \eqn{\widehat\sigma^2_{M_{CQ}}=
#' \frac{2}{n_1(n_1-1)}\widehat{\text{tr}(\mathbf{\Sigma}_1^2)}+
#' \frac{2}{n_2(n_2-1)}\widehat{\text{tr}(\mathbf{\Sigma}_2^2)}+
#' \frac{4}{n_1n_2}\widehat{\text{tr}(\mathbf{\Sigma}_1\mathbf{\Sigma}_2)}.}
#' The explicit formulas of \eqn{\widehat{\text{tr}(\mathbf{\Sigma}_1^2)}},
#' \eqn{\widehat{\text{tr}(\mathbf{\Sigma}_2^2)}}, and
#' \eqn{\widehat{\text{tr}(\mathbf{\Sigma}_1\mathbf{\Sigma}_2)}}
#' can be found in Section 3 of Chen and Qin (2010).
#' With some regularity conditions, under the null hypothesis
#' \eqn{H_{0m}: \boldsymbol{\mu}_1 = \boldsymbol{\mu}_2},
#' the test statistic \eqn{M_{CQ}} converges in distribution to a standard normal distribution
#' as \eqn{n_1, n_2, p \rightarrow \infty}.
#' The asymptotic \eqn{p}-value is obtained by
#' \deqn{p_{CQ} = 1-\Phi(M_{CQ}/\hat\sigma_{M_{CQ}}),}
#' where \eqn{\Phi(\cdot)} is the cdf of the standard normal distribution.
#' @import stats
#' @usage
#' meantest.cq(dataX,dataY)
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
#' @examples
#' n1 = 100; n2 = 100; pp = 500
#' set.seed(1)
#' X = matrix(rnorm(n1*pp), nrow=n1, ncol=pp)
#' Y = matrix(rnorm(n2*pp), nrow=n2, ncol=pp)
#' meantest.cq(X,Y)

meantest.cq <- function(dataX,dataY)
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
  n1=nrow(X);n2=nrow(Y);p=ncol(X)
  XXT=X%*%t(X);YYT=Y%*%t(Y);XYT=X%*%t(Y);
  A.1=(sum(XXT)-tr(XXT))/(n1*(n1-1));
  B.1=(sum(YYT)-tr(YYT))/(n2*(n2-1));
  C.1=(sum(XYT))/(n1*n2);
  T.value = A.1 + B.1 - 2*C.1
  tr1 = hat_tr_sigma(X)
  tr2 = hat_tr_sigma(Y)
  tr3 = hat_tr_inter(X,Y)
  sd.est = sqrt(2*tr1/(n1*(n1-1))+2*tr2/(n2*(n2-1))+4*tr3/(n1*n2))
  stat_mean_CQ = T.value/sd.est
  pval_mean_CQ = 1-pnorm(stat_mean_CQ)
  return(list(stat = stat_mean_CQ,
              pval = pval_mean_CQ))
}

hat_tr_sigma <- function(sample)
{
  X=sample; n=nrow(X); p=ncol(X)
  XXT=X%*%t(X)
  diagvec = diag(XXT)
  Xuu = matrix(diagvec, n, n, byrow = F)
  Xvv = matrix(diagvec, n, n, byrow = T)
  xbar = apply(X, 2, mean)
  Xub = matrix(X%*%xbar, n, n, byrow = F)
  Xvb = matrix(X%*%xbar, n, n, byrow = T)
  mat1 = (n-1)/(n-2)*XXT + 1/(n-2)*Xuu - n/(n-2)*Xub
  mat2 = (n-1)/(n-2)*XXT + 1/(n-2)*Xvv - n/(n-2)*Xvb
  est = sum(mat1*mat2)-sum(diag(mat1*mat2))
  ans = est/(n*(n-1))
  return(ans)
}

hat_tr_inter <- function(sample1,sample2)
{
  X = sample1; Y = sample2
  n1=nrow(X);n2=nrow(Y);p=ncol(Y);
  xbar = apply(X, 2, mean); ybar = apply(Y, 2, mean)
  XYT = X%*%t(Y)
  XYu = matrix(X%*%ybar, n1, n2, byrow = F)
  YXv = matrix(Y%*%xbar, n1, n2, byrow = T)
  mat1 = n2/(n2-1)*XYT - n2/(n2-1)*XYu
  mat2 = n1/(n1-1)*XYT - n1/(n1-1)*YXv
  est = sum(mat1*mat2)
  ans = est/(n1*n2)
  return(ans)
}
