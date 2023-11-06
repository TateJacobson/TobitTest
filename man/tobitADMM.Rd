\name{tobitADMM}
\alias{tobitADMM}
\title{Fit the partial penalized Tobit model}

\description{
  Fits the partial penalized Tobit model with a SCAD penalty using an ADMM algorithm.
}

\usage{
tobitADMM(x, y, left, C, M, t, a = 3.7, 
            nlambda = 100, lambda.factor = 0.01, lambda = NULL, 
            rho = 1, abs_tol = 1e-4, rel_tol = 1e-2, nr_tol = 1e-6, theta_tol = 1e-6, 
            full.model = FALSE, maxit = 1000)
}

\arguments{
  \item{x}{numeric predictor matrix with \code{n} rows and \code{p} columns, where each row corresponds to an observation and each column corresponds to a predictor}
  \item{y}{numeric response vector}
  \item{left}{left-censoring value for the response. Default is 0.}
  \item{C}{a numeric matrix, the constraint matrix under the linear null hypothesis you want to test}
  \item{M}{a numeric vector, the indices of the coefficients you want to test. Must be ordered from low to high.}
  \item{t}{a numeric vector, the vector defining the linear null hypothesis you want to test}
  \item{a}{a SCAD penalty hyperparameter. Default is 3.7.}
  \item{nlambda}{if \code{lambda} is not provided, this is the number of \eqn{\lambda}{\lambda} values used for the solution path. Default is 100.}
  \item{lambda.factor}{if \code{lambda} is not provided, this is the factor multiplied by \eqn{\lambda_{\max}}{\lambdamax}, the largest \eqn{\lambda}{\lambda} value in the automatically generated sequence for the solution path, to get \eqn{\lambda_{\min}}{\lambdamin}, the smallest \eqn{\lambda}{\lambda} value for the solution path. Default is 0.01.}
  \item{lambda}{a sequence of \eqn{\lambda} penalty parameter values for the SCAD penalty}
  \item{rho}{tuning parameter for ADMM algorithm. Default is 1.}
  \item{abs_tol}{absolute tolerance for declaring ADMM convergence. Default is 1e-4.}
  \item{rel_tol}{relative tolerance for declaring ADMM convergence. Default is 1e-2.}
  \item{nr_tol}{tolerance for declaring convergence for Newton Raphson update steps within the ADMM algorithm. Default is 1e-6.}
  \item{theta_tol}{tolerance for declaring convergence for \eqn{\theta}{\theta} updates within the ADMM algorithm. Default is 1e-6.}
  \item{full.model}{boolean, if TRUE, then \code{tobitADMM} fits the full model rather than the reduced model. Default is FALSE.}
  \item{maxit}{maximum number of ADMM iterations. Default is 1000.}
}

\value{
  An object with S3 class \dQuote{\code{tobitADMM}}.
   \item{b0}{vector of estimated intercept values of length \code{length(lambda)}}
   \item{beta}{matrix of estimated coefficient values with \code{p} rows and \code{length(lambda)} columns}
   \item{sigma}{vector of estimated sigma values of length \code{length(lambda)}}
   \item{lambda}{user-provided or automatically generated \code{lambda} path}
   \item{left}{user-provided \code{left} argument}
   \item{C}{user-provided \code{C} argument}
   \item{M}{user-provided \code{M} argument}
   \item{t}{user-provided \code{t} argument}
   \item{full.model}{user-provided \code{full.model} argument}
}

\details{
  The \code{tobitADMM} function fits the partial penalized Tobit model with a SCAD penalty using an alternating direction method of multipliers (ADMM) algorithm.
  The function can be used to fit both the full model and the reduced model, which satisfies the linear constraints of the null hypothesis: \eqn{H_0: C \beta_M = t}. See Jacobson and Zou (\emph{accepted}) for details.
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
}
