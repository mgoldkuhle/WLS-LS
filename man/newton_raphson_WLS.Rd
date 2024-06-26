% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_fit.R
\name{newton_raphson_WLS}
\alias{newton_raphson_WLS}
\title{Fisher Scoring/Newton-Raphson for the \code{LocationScaleRegressionWLS} model class}
\usage{
newton_raphson_WLS(
  model,
  abstol = 1e-04,
  maxit = 1000,
  verbose = FALSE,
  careful = FALSE
)
}
\arguments{
\item{model}{A \code{\link{LocationScaleRegressionWLS}} object.}

\item{abstol}{The absolute convergence tolerance. The algorithm stops if the
absolute value of the gradient drops below this value.}

\item{maxit}{The maximum number of iterations.}

\item{verbose}{Whether to print the progress of the algorithm.}

\item{diff}{The stepsize for the differential quotient.}
}
\value{
The updated model, invisibly.
}
\description{
This function optimizes the log-likelihood of the given location-scale
regression model by Newton-Raphson or Fisher Scoring. The Fisher Information
(derivative of the gradient) is approximated by the differential quotient.
The function has a side effect on the \code{model} object.
}
\examples{
y <- rnorm(30)
model <- LocationScaleRegressionWLS$new(y ~ 1)
newton_raphson_WLS(model)

}
