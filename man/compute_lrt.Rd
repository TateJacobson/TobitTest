\name{compute_lrt}
\alias{compute_lrt}
\title{Compute the partial penalized likelihood ratio test statistic}

\description{
Computes the likelihood ratio test statistic for a partial penalized Tobit model.
}

\usage{
compute_lrt(tADMM.full, tADMM.null, x, y, ic = c("max", "GIC", "BIC"))
}

\arguments{
\item{tADMM.full}{a fitted \code{tobitADMM} object for the full model. Must have been fit with \code{full.model = TRUE}. }
\item{tADMM.null}{a fitted \code{tobitADMM} object for the reduced model. Must have been fit with \code{full.model = FALSE}. }
\item{x}{numeric predictor matrix with \code{n} rows and \code{p} columns, where each row corresponds to an observation and each column corresponds to a predictor}
\item{y}{numeric response vector}
\item{ic}{a string specifying which information criterion penalty coefficient to use. Default is \code{max}, the maximum of the GIC and BIC penalty coefficients.}
}
    
\value{
\item{TL}{the computed partial penalized likelihood ratio test statistic}
\item{selected_model}{the model selected using \code{pick_lambda_tobit}}
\item{tADMM.full}{the user-provided \code{tADMM.full} argument}
\item{tADMM.null}{the user-provided \code{tADMM.null} argument}
}

\details{
    This function first selects the penalty parameter \eqn{\lambda} for the user-provided \code{tADMM.full} and \code{tADMM.null} objects using \code{pick_lambda_tobit}.
    It then computes the partial penalized likelihood ratio test statistic for the selected model.
    See Jacobson and Zou (\emph{accepted}) for details on how \eqn{\lambda} is selected and the partial penalized likelihood ratio test statistic is computed.
}

\author{
Tate Jacobson
}

\references{
Jacobson, T. and Zou, H. (\emph{accepted}) ``Linear Hypothesis Testing for High Dimensional Tobit Models,'' \emph{Statistica Sinica}.
}

\examples{
n = 100
p = 10

x = vapply(1:p, function(i) rnorm(n, 0, sd = 1), FUN.VALUE = numeric(n))

b0 = 1
beta = c(2, -2, rep(0, p-2))
epsilon = rnorm(n, 0, 1)

left = 1

ystar = rep(b0,n) + x\%*\%beta + epsilon
y = pmax(ystar, left)

M = c(1,2)
t = c(0)
C = matrix(c(1,1), nrow = 1, ncol = 2)

tadmm.null = tobitADMM(x = x, y = y, left = left, C = C, M = M, t = t, full.model = FALSE)
tadmm.full = tobitADMM(x = x, y = y, left = left, C = C, M = M, t = t, full.model = TRUE)
lrt = compute_lrt(tADMM.full = tadmm.full, tADMM.null = tadmm.null, x = x, y = y)
}