% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_cjs.R
\name{fit_s4t_cjs_ml}
\alias{fit_s4t_cjs_ml}
\title{Fit space-for-time mark-recapture model in likelihood framework}
\usage{
fit_s4t_cjs_ml(
  p_formula,
  theta_formula,
  ageclass_formula,
  cov_p = NULL,
  cov_theta = NULL,
  groups = NULL,
  s4t_ch,
  ndeps = 0.001,
  lmm = 5,
  maxit = 500,
  fixed_age = FALSE
)
}
\arguments{
\item{p_formula}{an object of class "formula" for the formula for detection probabilities.}

\item{theta_formula}{an object of class "formula" for the formula for transition probabilities.}

\item{ageclass_formula}{an object of class "formula" for the ageclass sub-model}

\item{cov_p}{a \code{data.frame} or \code{list} of \code{data.frame}'s containing the covariates for p. See details.}

\item{cov_theta}{a \code{data.frame} or \code{list} of \code{data.frame}'s containing the covariates for theta. See details.}

\item{groups}{a \code{character} vector containing the names of the covariates that comprise the groups.}

\item{s4t_ch}{a \code{s4t_ch} object}

\item{ndeps}{a \code{numeric} value .....}

\item{lmm}{an \code{integer} of ....}

\item{maxit}{an \code{integer} of the max....}

\item{fixed_age}{a \code{logical} object that determines whether the ageclass model will be run
as a separate model (TRUE) or whether it is estimated along with the CJS model (FALSE).}
}
\value{
a \code{s4t_cjs} object.
}
\description{
Uses optim to fit the age-specific space-for-time model.
}
\examples{
# don't run

}
