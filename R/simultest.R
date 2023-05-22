#' Two-sample simultaneous tests on high-dimensional mean vectors and
#' covariance matrices
#' @description
#' This function implements six two-sample simultaneous tests
#' on high-dimensional mean vectors and covariance matrices.
#' Let \eqn{\mathbf{X} \in \mathbb{R}^p} and \eqn{\mathbf{Y} \in \mathbb{R}^p}
#' be two \eqn{p}-dimensional populations with mean vectors
#' \eqn{(\boldsymbol{\mu}_1, \boldsymbol{\mu}_2)} and covariance matrices
#' \eqn{(\mathbf{\Sigma}_1, \mathbf{\Sigma}_2)}, respectively.
#' The problem of interest is the simultaneous inference on the equality of
#' mean vectors and covariance matrices of the two populations:
#' \deqn{H_0: \boldsymbol{\mu}_1 = \boldsymbol{\mu}_2 \ \text{ and }
#' \ \mathbf{\Sigma}_1 = \mathbf{\Sigma}_2. }
#' Suppose \eqn{\{\mathbf{X}_1, \ldots, \mathbf{X}_{n_1}\}} are i.i.d.
#' copies of \eqn{\mathbf{X}}, and \eqn{\{\mathbf{Y}_1, \ldots, \mathbf{Y}_{n_2}\}}
#' are i.i.d. copies of \eqn{\mathbf{Y}}. We denote
#' \code{dataX=}\eqn{(\mathbf{X}_1, \ldots, \mathbf{X}_{n_1})^\top\in\mathbb{R}^{n_1\times p}}
#' and \code{dataY=}\eqn{(\mathbf{Y}_1, \ldots, \mathbf{Y}_{n_2})^\top\in\mathbb{R}^{n_2\times p}}.
#' @import stats
#' @usage
#' simultest(dataX, dataY, method='pe.fisher', delta_mean=NULL, delta_cov=NULL)
#' @param dataX an \eqn{n_1} by \eqn{p} data matrix
#' @param dataY an \eqn{n_2} by \eqn{p} data matrix
#' @param method the method type (default = \code{'pe.fisher'}); chosen from
#' * \code{'cauchy'}: the simultaneous test via Cauchy combination;\cr
#'                    see \code{\link{simultest.cauchy}} for details.
#' * \code{'chisq'}: the simultaneous test via chi-squared approximation;\cr
#'                    see \code{\link{simultest.chisq}} for details.
#' * \code{'fisher'}: the simultaneous test via Fisher's combination;\cr
#'                    see \code{\link{simultest.fisher}} for details.
#' * \code{'pe.cauchy'}: the PE simultaneous test via Cauchy combination;\cr
#'                    see \code{\link{simultest.pe.cauchy}} for details.
#' * \code{'pe.chisq'}: the PE simultaneous test via chi-squared approximation;\cr
#'                    see \code{\link{simultest.pe.chisq}} for details.
#' * \code{'pe.fisher'}: the PE simultaneous test via Fisher's combination;\cr
#'                    see \code{\link{simultest.pe.fisher}} for details.
#'
#' @param delta_mean the thresholding value used in the construction of
#' the PE component for the mean test statistic. It is needed only in PE methods such as
#' \code{method='pe.cauchy'}, \code{method='pe.chisq'}, and
#' \code{method='pe.fisher'}; see \code{\link{simultest.pe.cauchy}}, \cr
#' \code{\link{simultest.pe.chisq}},
#' and \code{\link{simultest.pe.fisher}}
#' for details. The default is NULL.
#' @param delta_cov the thresholding value used in the construction of
#' the PE component for the covariance test statistic. It is needed only in PE methods such as
#' \code{method='pe.cauchy'}, \code{method='pe.chisq'}, and
#' \code{method='pe.fisher'}; see {\code{\link{simultest.pe.cauchy}}}, \cr
#' \code{\link{simultest.pe.chisq}},
#' and \code{\link{simultest.pe.fisher}}
#' for details. The default is NULL.
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
#' simultest(X,Y)

simultest <- function(dataX, dataY, method = "pe.fisher",
                     delta_mean=NULL, delta_cov=NULL)
{
  X=dataX; Y=dataY;
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

  if(method != "cauchy" & method != "chisq" & method != "fisher" &
     method != "pe.cauchy" & method != "pe.chisq" & method != "pe.fisher"){
    stop("'method' must be one of 'cauchy', 'chisq', 'fisher',
          'pe.cauchy', 'pe.chisq', or 'pe.fisher'!")
  }

  if(method == "pe.fisher")
  {
    res <- simultest.pe.fisher(X,Y,delta_mean,delta_cov)
  }

  if(method == "cauchy")
  {
    res <- simultest.cauchy(X,Y)
  }

  if(method == "chisq")
  {
    res <- simultest.chisq(X,Y)
  }

  if(method == "fisher")
  {
    res <- simultest.fisher(X,Y)
  }

  if(method == "pe.cauchy")
  {
    res <- simultest.pe.cauchy(X,Y,delta_mean,delta_cov)
  }

  if(method == "pe.chisq")
  {
    res <- simultest.pe.chisq(X,Y,delta_mean,delta_cov)
  }

  output <- c(list(method = method), res)
  return(output)
}



