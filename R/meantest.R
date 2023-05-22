#' Two-sample mean tests for high-dimensional data
#' @description
#' This function implements five two-sample mean tests on high-dimensional
#' mean vectors.
#' Let \eqn{\mathbf{X} \in \mathbb{R}^p} and \eqn{\mathbf{Y} \in \mathbb{R}^p}
#' be two \eqn{p}-dimensional populations with mean vectors
#' \eqn{(\boldsymbol{\mu}_1, \boldsymbol{\mu}_2)} and covariance matrices
#' \eqn{(\mathbf{\Sigma}_1, \mathbf{\Sigma}_2)}, respectively.
#' The problem of interest is to test the equality of the two
#' mean vectors of the two populations:
#' \deqn{H_{0m}: \boldsymbol{\mu}_1 = \boldsymbol{\mu}_2.}
#' Suppose \eqn{\{\mathbf{X}_1, \ldots, \mathbf{X}_{n_1}\}} are i.i.d.
#' copies of \eqn{\mathbf{X}}, and \eqn{\{\mathbf{Y}_1, \ldots, \mathbf{Y}_{n_2}\}}
#' are i.i.d. copies of \eqn{\mathbf{Y}}. We denote
#' \code{dataX=}\eqn{(\mathbf{X}_1, \ldots, \mathbf{X}_{n_1})^\top\in\mathbb{R}^{n_1\times p}}
#' and \code{dataY=}\eqn{(\mathbf{Y}_1, \ldots, \mathbf{Y}_{n_2})^\top\in\mathbb{R}^{n_2\times p}}.
#' @import stats
#' @usage
#' meantest(dataX,dataY,method='pe.comp',delta=NULL)
#' @param dataX an \eqn{n_1} by \eqn{p} data matrix
#' @param dataY an \eqn{n_2} by \eqn{p} data matrix
#' @param method the method type (default = \code{'pe.comp'});
#'  chosen from
#' * \code{'clx'}: the \eqn{l_\infty}-norm-based mean test, proposed in Cai et al. (2014); \cr
#'                 see \code{\link{meantest.clx}} for details.
#' * \code{'cq'}: the \eqn{l_2}-norm-based mean test, proposed in Chen and Qin (2010); \cr
#'                 see \code{\link{meantest.cq}} for details.
#' * \code{'pe.cauchy'}: the PE mean test via Cauchy combination; \cr
#'                    see \code{\link{meantest.pe.cauchy}} for details.
#' * \code{'pe.comp'}: the PE mean test via the construction of PE components; \cr
#'                    see \code{\link{meantest.pe.comp}} for details.
#' * \code{'pe.fisher'}: the PE mean test via Fisher's combination; \cr
#'                    see \code{\link{meantest.pe.fisher}} for details.
#' @param delta This is needed only in \code{method='pe.comp'};
#' see \code{\link{meantest.pe.comp}} for details.
#' The default is NULL.
#' @return
#' `method` the method type
#'
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
#'
#' Yu, X., Li, D., Xue, L., and Li, R. (2022). Power-enhanced simultaneous test
#' of high-dimensional mean vectors and covariance matrices with application
#' to gene-set testing. \emph{Journal of the American Statistical Association},
#' (in press):1–14.
#' @examples
#' n1 = 100; n2 = 100; pp = 500
#' set.seed(1)
#' X = matrix(rnorm(n1*pp), nrow=n1, ncol=pp)
#' Y = matrix(rnorm(n2*pp), nrow=n2, ncol=pp)
#' meantest(X,Y)

meantest <- function(dataX,dataY, method="pe.comp", delta=NULL)
{
  X=dataX; Y=dataY
  n1=nrow(X);n2=nrow(Y);p=ncol(X);p2=ncol(Y)
  #if(p!= p2)
  #{
  #  stop(" The data dimensions ncol(X) and ncol(Y) do not match!")
  #}

  #if(p <= 30)
  #{
  #  warning(paste0("These methods are designed for high-dimensional data!
  #                 The data dimension (p=", p, ") is small in the input data,
  #                 which may results in an inflated Type-I error rate."))
  #}

  if(method != "clx" & method != "cq" & method != "pe.cauchy" &
     method != "pe.comp" & method != "pe.fisher"){
    stop("'method' must be one of 'clx', 'cq', 'pe.cauchy',
          'pe.comp', or 'pe.fisiher'!")
  }

  if(method == "pe.comp")
  {
    res <- meantest.pe.comp(X,Y, delta)
  }

  if(method == "clx")
  {
    res <- meantest.clx(X,Y)
  }

  if(method == "cq")
  {
    res <- meantest.cq(X,Y)
  }

  if(method == "pe.cauchy")
  {
    res <- meantest.pe.cauchy(X,Y)
  }

  if(method == "pe.fisher")
  {
    res <- meantest.pe.fisher(X,Y)
  }

  output <- c(list(method = method), res)
  return(output)


}
