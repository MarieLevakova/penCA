% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pen_svd.R
\name{pen.svd}
\alias{pen.svd}
\title{Fit VEC model with penalty imposed on SVD of Pi}
\usage{
pen.svd(
  Y,
  X,
  nrank,
  ic.type = c("AIC", "BIC", "HQ", "CV"),
  control = list(),
  n.cv = 20
)
}
\arguments{
\item{Y}{Matrix of differences of the process}

\item{X}{Matrix of lagged values of the process}

\item{nrank}{cointegration rank}

\item{control}{parameters specifying the iterative algorithm}

\item{n.cv}{number of repetitions of the crossvalidation procedure}

\item{ic}{Criterion for choosing optimal penalty value, the options are
\code{"CV"} (crossvalidated error), \code{"AIC"} (Akaike information criterion),
\code{"BIC"} (Bayes information criterion) and \code{"HQ"} (Hannan-Quinn criterion).}
}
\value{
A list containing the following components:
}
\description{
This function loads data matrix X and fits a VEC model,
where penalty is imposed on singular values and vectors
in the singular vector decomposition of Pi
}
\section{Slots}{

\describe{
\item{\code{Y}}{matrix of differences of the process}

\item{\code{X}}{matrix of lagged values of the process}

\item{\code{U.path}}{the regularization path of the matrix of left single vectors}

\item{\code{V.path}}{the regularization path of the matrix of right single vectors}

\item{\code{D.path}}{the regularization path of the singular values}

\item{\code{ic.path}}{the regularization path of the objective function}

\item{\code{U}}{the final left singular vectors}

\item{\code{V}}{the final right singular vectors}

\item{\code{D}}{the final singular values}

\item{\code{rank}}{the final cointegration rank}

\item{\code{PI}}{the final estimate of the matrix Pi}
}}

