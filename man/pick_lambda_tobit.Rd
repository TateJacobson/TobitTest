\name{pick_lambda_tobit}
\alias{pick_lambda_tobit}
\title{Tune the penalty parameter for a partial penalized Tobit model}

\description{
Selects the \eqn{\lambda} penalty parameter for a partial penalized Tobit model using a generalized information criterion (GIC).
}

\usage{
pick_lambda_tobit(tADMM, x, y, ic = c("max", "GIC", "BIC"))
}

\arguments{
\item{tADMM}{a fitted tobitADMM object}
\item{x}{numeric predictor matrix with \code{n} rows and \code{p} columns, where each row corresponds to an observation and each column corresponds to a predictor}
\item{y}{numeric response vector}
\item{ic}{a string specifying which information criterion penalty coefficient to use to select \eqn{\lambda}. Default is \code{max}, the maximum of the GIC and BIC penalty coefficients.}
}

\value{
\item{gic}{the information criterion values computed along the \code{lambda} path from tADMM}
\item{cn}{the information criterion penalty coefficient value}
\item{neg.loglik}{the negative log likelihood computed along the \code{lambda} path from tADMM}
\item{beta.zero.norm}{the number of nonzero estimated coefficients for each model along the \code{lambda} path from tADMM}
\item{ind}{the index of the selected model within the solution path from tADMM}
\item{lambda}{the \eqn{\lambda} penalty value selected using the information criterion}
\item{b0}{the estimated intercept value from the selected model}
\item{beta}{the estimated coefficient vector from the selected model}
\item{sigma}{the estimated sigma value from the selected model}
}

\details{
See Jacobson and Zou (\emph{accepted}) for details on how \eqn{\lambda} is selected.
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

tadmm = tobitADMM(x = x, y = y, left = left, C = C, M = M, t = t)
selected_tadmm = pick_lambda_tobit(tADMM = tadmm, x = x, y = y)
}
