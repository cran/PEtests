% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meantest.pe.comp.R
\name{meantest.pe.comp}
\alias{meantest.pe.comp}
\title{Two-sample PE mean test for high-dimensional data via PE component}
\usage{
meantest.pe.comp(dataX,dataY,delta=NULL)
}
\arguments{
\item{dataX}{an \eqn{n_1} by \eqn{p} data matrix}

\item{dataY}{an \eqn{n_2} by \eqn{p} data matrix}

\item{delta}{a scalar; the thresholding value used in the construction of
the PE component. If not specified, the function uses a default value
\eqn{\delta_{mean}=2\log(\log (n_1+n_2))\log p}.}
}
\value{
\code{stat} the value of test statistic

\code{pval} the p-value for the test.
}
\description{
This function implements the two-sample PE mean via the
construction of the PE component. Let \eqn{M_{CQ}/\hat\sigma_{M_{CQ}}}
denote the \eqn{l_2}-norm-based mean test statistic
(see \code{\link{meantest.cq}} for details).
The PE component is constructed by
\deqn{J_m = \sqrt{p}\sum_{i=1}^p M_i\widehat\nu^{-1/2}_i
\mathcal{I}\{ \sqrt{2}M_i\widehat\nu^{-1/2}_i + 1 >  \delta_{mean} \}, }
where \eqn{\delta_{mean}} is a threshold for the screening procedure,
recommended to take the value of \eqn{\delta_{mean}=2\log(\log (n_1+n_2))\log p}.
The explicit forms of \eqn{M_{i}} and \eqn{\widehat\nu_{j}}
can be found in Section 3.1 of Yu et al. (2022).
The PE covariance test statistic is defined as
\deqn{M_{PE}=M_{CQ}/\hat\sigma_{M_{CQ}}+J_m.}
With some regularity conditions, under the null hypothesis
\eqn{H_{0m}:  \boldsymbol{\mu}_1 = \boldsymbol{\mu}_2},
the test statistic \eqn{M_{PE}} converges in distribution to
a standard normal distribution as \eqn{n_1, n_2, p \rightarrow \infty}.
The asymptotic  \eqn{p}-value is obtained by
\deqn{p\text{-value}= 1-\Phi(M_{PE}),}
where \eqn{\Phi(\cdot)} is the cdf of the standard normal distribution.
}
\examples{
n1 = 100; n2 = 100; pp = 500
set.seed(1)
X = matrix(rnorm(n1*pp), nrow=n1, ncol=pp)
Y = matrix(rnorm(n2*pp), nrow=n2, ncol=pp)
meantest.pe.comp(X,Y)
}
\references{
Yu, X., Li, D., Xue, L., and Li, R. (2022). Power-enhanced simultaneous test
of high-dimensional mean vectors and covariance matrices with application
to gene-set testing. \emph{Journal of the American Statistical Association},
(in press):1–14.
}
