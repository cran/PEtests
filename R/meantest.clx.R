#' Two-sample high-dimensional mean test (Cai, Liu and Xia, 2014)
#' @description
#' This function implements the two-sample \eqn{l_\infty}-norm-based
#' high-dimensional mean test proposed in Cai, Liu and Xia (2014).
#' Suppose \eqn{\{\mathbf{X}_1, \ldots, \mathbf{X}_{n_1}\}} are i.i.d.
#' copies of \eqn{\mathbf{X}}, and \eqn{\{\mathbf{Y}_1, \ldots, \mathbf{Y}_{n_2}\}}
#' are i.i.d. copies of \eqn{\mathbf{Y}}.
#' The test statistic is defined as
#' \deqn{
#' M_{CLX}=\frac{n_1n_2}{n_1+n_2}\max_{1\leq j\leq p}
#'  \frac{(\bar{X_j}-\bar{Y_j})^2}
#'  {\frac{1}{n_1+n_2} [\sum_{u=1}^{n_1} (X_{uj}-\bar{X_j})^2+\sum_{v=1}^{n_2} (Y_{vj}-\bar{Y_j})^2]  }
#'  }
#' With some regularity conditions, under the null hypothesis \eqn{H_{0c}: \mathbf{\Sigma}_1 = \mathbf{\Sigma}_2},
#' the test statistic \eqn{M_{CLX}-2\log p+\log\log p} converges in distribution to
#' a Gumbel distribution \eqn{G_{mean}(x) = \exp(-\frac{1}{\sqrt{\pi}}\exp(-\frac{x}{2}))}
#' as \eqn{n_1, n_2, p \rightarrow \infty}.
#' The asymptotic \eqn{p}-value is obtained by
#' \deqn{p_{CLX} = 1-G_{mean}(M_{CLX}-2\log p+\log\log p).}
#' @import stats
#' @usage
#' meantest.clx(dataX,dataY)
#' @param dataX an \eqn{n_1} by \eqn{p} data matrix
#' @param dataY an \eqn{n_2} by \eqn{p} data matrix
#' @return
#' `stat` the value of test statistic
#'
#' `pval` the p-value for the test.
#' @export
#' @references
#' Cai, T. T., Liu, W., and Xia, Y. (2014). Two-sample test of high dimensional
#'  means under dependence. \emph{Journal of the Royal Statistical Society:
#'  Series B: Statistical Methodology}, 76(2):349–372.
#' @examples
#' n1 = 100; n2 = 100; pp = 500
#' set.seed(1)
#' X = matrix(rnorm(n1*pp), nrow=n1, ncol=pp)
#' Y = matrix(rnorm(n2*pp), nrow=n2, ncol=pp)
#' meantest.clx(X,Y)

meantest.clx <- function(dataX,dataY)
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
  xbar = apply(X, 2, mean)
  ybar = apply(Y, 2, mean)
  deltaVec = xbar-ybar
  pooledCov = ((n1-1)*cov(X)+(n2-1)*cov(Y))/(n1+n2)
  pooledVar = diag(pooledCov)
  M.value = n1*n2/(n1+n2)*max(deltaVec*deltaVec/pooledVar)
  stat_mean_CLX = M.value-2*log(p)+log(log(p))
  pval_mean_CLX = 1-F.extreme.mean(stat_mean_CLX)
  return(list(stat = stat_mean_CLX,
              pval = pval_mean_CLX))
}

F.extreme.mean <- function(x)
{
  return(exp(-exp(-x/2)/(sqrt(pi))))
}



