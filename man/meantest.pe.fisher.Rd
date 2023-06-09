% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meantest.pe.fisher.R
\name{meantest.pe.fisher}
\alias{meantest.pe.fisher}
\title{Two-sample PE mean test for high-dimensional data via Fisher's combination}
\usage{
meantest.pe.fisher(dataX,dataY)
}
\arguments{
\item{dataX}{an \eqn{n_1} by \eqn{p} data matrix}

\item{dataY}{an \eqn{n_2} by \eqn{p} data matrix}
}
\value{
\code{stat} the value of test statistic

\code{pval} the p-value for the test.
}
\description{
This function implements the two-sample PE covariance test via
Fisher's combination.
Suppose \eqn{\{\mathbf{X}_1, \ldots, \mathbf{X}_{n_1}\}} are i.i.d.
copies of \eqn{\mathbf{X}}, and \eqn{\{\mathbf{Y}_1, \ldots, \mathbf{Y}_{n_2}\}}
are i.i.d. copies of \eqn{\mathbf{Y}}.
Let \eqn{p_{CQ}} and \eqn{p_{CLX}} denote the \eqn{p}-values associated with
the \eqn{l_2}-norm-based covariance test (see \code{\link{meantest.cq}} for details)
and the \eqn{l_\infty}-norm-based covariance test
(see \code{\link{meantest.clx}} for details), respectively.
The PE covariance test via Fisher's combination is defined as
\deqn{M_{Fisher} = -2\log(p_{CQ})-2\log(p_{CLX}).}
It has been proved that with some regularity conditions, under the null hypothesis
\eqn{H_{0m}: \boldsymbol{\mu}_1 = \boldsymbol{\mu}_2,}
the two tests are asymptotically independent as \eqn{n_1, n_2, p\rightarrow \infty},
and therefore \eqn{M_{Fisher}} asymptotically converges in distribution to a \eqn{\chi_4^2} distribution.
The asymptotic \eqn{p}-value is obtained by
\deqn{p\text{-value} = 1-F_{\chi_4^2}(M_{Fisher}),}
where \eqn{F_{\chi_4^2}(\cdot)} is the cdf of the \eqn{\chi_4^2} distribution.
}
\examples{
n1 = 100; n2 = 100; pp = 500
set.seed(1)
X = matrix(rnorm(n1*pp), nrow=n1, ncol=pp)
Y = matrix(rnorm(n2*pp), nrow=n2, ncol=pp)
meantest.pe.fisher(X,Y)
}
\references{
Chen, S. X. and Qin, Y. L. (2010). A two-sample test for high-dimensional data
with applications to gene-set testing.
\emph{Annals of Statistics}, 38(2):808–835.

Cai, T. T., Liu, W., and Xia, Y. (2014). Two-sample test of high dimensional
means under dependence. \emph{Journal of the Royal Statistical Society:
Series B: Statistical Methodology}, 76(2):349–372.
}
