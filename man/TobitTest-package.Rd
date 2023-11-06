\name{TobitTest-package}
\alias{TobitTest-package}
\alias{TobitTest}
\docType{package}
\title{
  Linear Hypothesis Testing for High Dimensional Tobit Models
}
\description{
  Compute partial penalized estimators and partial penalized Wald, score, and likelihood ratio test statistics for high dimensional Tobit models.
}
\details{
  Compute partial penalized estimators for high dimensional Tobit models with \code{\link{tobitADMM}}.
  Select the penalty parameter \eqn{\lambda} with \code{\link{pick_lambda_tobit}}.
  Compute the partial penalized Wald, score, and likelihood ratio test statistics with \code{\link{compute_wald}}, \code{\link{compute_score}}, and \code{\link{compute_lrt}}.
}
\author{
Tate Jacobson

Maintainer: Tate Jacobson <tate.jacobson@oregonstate.edu>
}
\references{
  Jacobson, T. and Zou, H. (\emph{accepted}) ``Linear Hypothesis Testing for High Dimensional Tobit Models,'' \emph{Statistica Sinica}.
}
\keyword{ censored regression, high dimensional inference }
