#' Two-sample covariance tests for high-dimensional data
#' @description
#' This function implements five two-sample covariance tests on high-dimensional
#' covariance matrices.
#' Let \eqn{\mathbf{X} \in \mathbb{R}^p} and \eqn{\mathbf{Y} \in \mathbb{R}^p}
#' be two \eqn{p}-dimensional populations with mean vectors
#' \eqn{(\boldsymbol{\mu}_1, \boldsymbol{\mu}_2)} and covariance matrices
#' \eqn{(\mathbf{\Sigma}_1, \mathbf{\Sigma}_2)}, respectively.
#' The problem of interest is to test the equality of the two
#' covariance matrices:
#' \deqn{H_{0c}: \mathbf{\Sigma}_1 = \mathbf{\Sigma}_2. }
#' Suppose \eqn{\{\mathbf{X}_1, \ldots, \mathbf{X}_{n_1}\}} are i.i.d.
#' copies of \eqn{\mathbf{X}}, and \eqn{\{\mathbf{Y}_1, \ldots, \mathbf{Y}_{n_2}\}}
#' are i.i.d. copies of \eqn{\mathbf{Y}}. We denote
#' \code{dataX=}\eqn{(\mathbf{X}_1, \ldots, \mathbf{X}_{n_1})^\top\in\mathbb{R}^{n_1\times p}}
#' and \code{dataY=}\eqn{(\mathbf{Y}_1, \ldots, \mathbf{Y}_{n_2})^\top\in\mathbb{R}^{n_2\times p}}.
#' @import stats
#' @usage
#' covtest(dataX,dataY,method='pe.comp',delta=NULL)
#' @param dataX an \eqn{n_1} by \eqn{p} data matrix
#' @param dataY an \eqn{n_2} by \eqn{p} data matrix
#' @param method the method type (default = \code{'pe.comp'});
#'  chosen from
#' * \code{'clx'}: the \eqn{l_\infty}-norm-based covariance test, proposed in Cai et al. (2013); \cr
#'                 see \code{\link{covtest.clx}} for details.
#' * \code{'lc'}: the \eqn{l_2}-norm-based covariance test, proposed in Li and Chen (2012); \cr
#'                 see \code{\link{covtest.lc}} for details.
#' * \code{'pe.cauchy'}: the PE covariance test via Cauchy combination; \cr
#'                    see \code{\link{covtest.pe.cauchy}} for details.
#' * \code{'pe.comp'}: the PE covariance test via the construction of PE components; \cr
#'                    see \code{\link{covtest.pe.comp}} for details.
#' * \code{'pe.fisher'}: the PE covariance test via Fisher's combination; \cr
#'                    see \code{\link{covtest.pe.fisher}} for details.
#' @param delta This is needed only in \code{method='pe.comp'};
#' see \code{\link{covtest.pe.comp}} for details.
#' The default is NULL.
#' @return
#' `method` the method type
#'
#' `stat` the value of test statistic
#'
#' `pval` the p-value for the test.
#' @export
#' @references
#'
#' Cai, T. T., Liu, W., and Xia, Y. (2013). Two-sample covariance matrix testing
#' and support recovery in high-dimensional and sparse settings.
#' \emph{Journal of the American Statistical Association}, 108(501):265–277.
#'
#' Li, J. and Chen, S. X. (2012). Two sample tests for high-dimensional
#' covariance matrices. \emph{The Annals of Statistics}, 40(2):908–940.
#'
#' Yu, X., Li, D., and Xue, L. (2022). Fisher’s combined probability test
#' for high-dimensional covariance matrices. \emph{Journal of the American
#' Statistical Association}, (in press):1–14.
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
#' covtest(X,Y)


covtest <- function(dataX,dataY,method="pe.comp",delta=NULL)
{
  n1=nrow(dataX);n2=nrow(dataY); p=ncol(dataX); p2=ncol(dataY)
  #if(p!= p2)
  #{
  #  stop(" The data dimensions ncol(dataX) and ncol(dataY) do not match!")
  #}

  #if(p <= 30)
  #{
  #  warning(paste0("These methods are designed for high-dimensional data!
  #                 The data dimension (p=", p ") is small in the input data,
  #                 which may results in an inflated Type-I error rate.")
  #}

  if(method != "clx" & method != "lc" & method != "pe.cauchy" &
     method != "pe.comp" & method != "pe.fisher"){
    stop("'method' must be one of 'clx', 'lc', 'pe.cauchy',
          'pe.comp', or 'pe.fisiher'!")
  }

  if(method == "pe.comp")
  {
    res <- covtest.pe.comp(dataX,dataY,delta)
  }

  if(method == "clx")
  {
    res <- covtest.clx(dataX,dataY)
  }

  if(method == "lc")
  {
    res <- covtest.lc(dataX,dataY)
  }

  if(method == "pe.cauchy")
  {
    res <- covtest.pe.cauchy(dataX,dataY)
  }

  if(method == "pe.fisher")
  {
    res <- covtest.pe.fisher(dataX,dataY)
  }

  output <- c(list(method = method), res)
  return(output)
}

