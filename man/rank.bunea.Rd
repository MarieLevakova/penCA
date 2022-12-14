% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rank.R
\name{rank.bunea}
\alias{rank.bunea}
\title{Determine the cointegration rank by the method of Bunea et al.}
\usage{
rank.bunea(X, figure = FALSE)
}
\arguments{
\item{X}{Matrix of lagged values of the process}

\item{figure}{logical value indicating if an illustrative figure should be plotted}
}
\value{
A list containing the following components:
}
\description{
This function loads data matrix X and determines the cointegration rank by
the method proposed by Bunea et al. using a rank penalty
}
\section{Slots}{

\describe{
\item{\code{r}}{cointegration rank}

\item{\code{lambdas}}{eigenvalues}

\item{\code{mu}}{threshold of statistically significant eigenvalues}
}}

