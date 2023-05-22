#' Two-sample high-dimensional covariance test (Cai, Liu and Xia, 2013)
#' @description
#' This function implements the two-sample \eqn{l_\infty}-norm-based
#' high-dimensional covariance test proposed in Cai, Liu and Xia (2013).
#' Suppose \eqn{\{\mathbf{X}_1, \ldots, \mathbf{X}_{n_1}\}} are i.i.d.
#' copies of \eqn{\mathbf{X}}, and \eqn{\{\mathbf{Y}_1, \ldots, \mathbf{Y}_{n_2}\}}
#' are i.i.d. copies of \eqn{\mathbf{Y}}. The test statistic is defined as
#' \deqn{T_{CLX} = \max_{1\leq i,j \leq p} \frac{(\hat\sigma_{ij1}-\hat\sigma_{ij2})^2}
#' {\hat\theta_{ij1}/n_1+\hat\theta_{ij2}/n_2},}
#' where \eqn{\hat\sigma_{ij1}} and \eqn{\hat\sigma_{ij2}} are the sample covariances,
#' and \eqn{\hat\theta_{ij1}/n_1+\hat\theta_{ij2}/n_2} estimates the variance of
#' \eqn{\hat{\sigma}_{ij1}-\hat{\sigma}_{ij2}}.
#' The explicit formulas of \eqn{\hat\sigma_{ij1}}, \eqn{\hat\sigma_{ij2}},
#' \eqn{\hat\theta_{ij1}} and \eqn{\hat\theta_{ij2}} can be found
#' in Section 2 of Cai, Liu and Xia (2013).
#' With some regularity conditions, under the null hypothesis \eqn{H_{0c}: \mathbf{\Sigma}_1 = \mathbf{\Sigma}_2},
#' the test statistic \eqn{T_{CLX}-4\log p+\log\log p} converges in distribution to
#' a Gumbel distribution \eqn{G_{cov}(x) = \exp(-\frac{1}{\sqrt{8\pi}}\exp(-\frac{x}{2}))}
#' as \eqn{n_1, n_2, p \rightarrow \infty}.
#' The asymptotic \eqn{p}-value is obtained by
#' \deqn{p_{CLX} = 1-G_{cov}(T_{CLX}-4\log p+\log\log p).}
#' @import stats
#' @usage
#' covtest.clx(dataX,dataY)
#' @param dataX an \eqn{n_1} by \eqn{p} data matrix
#' @param dataY an \eqn{n_2} by \eqn{p} data matrix
#' @return
#' `stat` the value of test statistic
#'
#' `pval` the p-value for the test.
#' @export
#' @references
#' Cai, T. T., Liu, W., and Xia, Y. (2013). Two-sample covariance matrix testing
#' and support recovery in high-dimensional and sparse settings.
#' \emph{Journal of the American Statistical Association}, 108(501):265–277.
#' @examples
#' n1 = 100; n2 = 100; pp = 500
#' set.seed(1)
#' X = matrix(rnorm(n1*pp), nrow=n1, ncol=pp)
#' Y = matrix(rnorm(n2*pp), nrow=n2, ncol=pp)
#' covtest.clx(X,Y)

covtest.clx <- function(dataX,dataY)
{
  X=dataX;Y=dataY
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

  Sigma.hat.1=((n1-1)/n1)*cov(X)
  Sigma.hat.2=((n2-1)/n2)*cov(Y)
  X.c=center(X);Y.c=center(Y);
  Theta.hat.1=matrix(NA,p,p)
  Theta.hat.1=(1/n1)*t(X.c^2)%*%(X.c^2)-(2/n1)*(Sigma.hat.1)*(t(X.c)%*%(X.c))+Sigma.hat.1^2
  Theta.hat.2=matrix(NA,p,p)
  Theta.hat.2=(1/n2)*t(Y.c^2)%*%(Y.c^2)-(2/n2)*(Sigma.hat.2)*(t(Y.c)%*%(Y.c))+Sigma.hat.2^2
  M.value=max(((Sigma.hat.1- Sigma.hat.2)^2)/(Theta.hat.1/n1+Theta.hat.2/n2))
  stat_cov_CLX = M.value-4*log(p)+log(log(p))
  pval_cov_CLX = 1-F.extreme.cov(stat_cov_CLX)
  return(list(stat = stat_cov_CLX,
              pval = pval_cov_CLX))
}


F.extreme.cov <- function(x)
{
  return(exp(-exp(-x/2)/(sqrt(8*pi))))
}




