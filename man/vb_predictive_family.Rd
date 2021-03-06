% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vb_predictive_family.R
\name{vb_predictive_family}
\alias{vb_predictive_family}
\title{Bayesian linear mixed model with multiple random effects for family data}
\usage{
vb_predictive_family(fit, epsilon = 1e-3, maf_beta = 0.5, maf_u = 0.5)
}
\arguments{
\item{fit}{\code{vb} fit from \code{vb_fit_family}}

\item{epsilon}{a threshold for the increase of the variance}

\item{maf_beta}{a filtering threshold for the weight of parameter beta}

\item{maf_u}{a filtering threshold for the weight of parameter U}
}
\description{
Bayesian linear mixed model with multiple random effects for family data
}
\examples{
fit <- vb_fit_family(...)
vb_predictive_family(fit, epsilon = 1e-3, maf_beta = 0.5, maf_u = 0.5)
}
