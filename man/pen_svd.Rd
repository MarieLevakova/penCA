\name{pen_svd}
\alias{pen_svd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reduced-rank regression with a sparse singular value decomposition
}
\description{
Fitting of a cointegrated VAR model with a sparse singular value decomposition using the iterative exclusive extraction algorithm.
}
\usage{
pen.svd(Y, X, nrank, ic.type = c("AIC", "BIC", "HQ"),
        control = list(), n.cv = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  matrix with differenced time series
  }
  \item{X}{
  matrix with lagged time series
  }
  \item{nrank}{
  cointegration rank
  }
  \item{ic.type}{
  Information criterion to use in selection of the tuning parameter
  }
  \item{control}{
  other parameters to control the fitting procedure (see Details)
  }
  \item{n.cv}{
  Number of repetitions of the cross-validation procedure
}
}
\details{
The model fitting can be controled through argument \code{control}.
The available elements include
\itemize{
   \item{\code{maxit}}: maximum number of iterations.
   \item{\code{epsilon}}: convergence tolerance.
   \item{\code{innerMaxit}}: maximum number of iterations for inner steps.
   \item{\code{innerEpsilon}}: convergence tolerance for inner steps.
   \item{\code{nlambda}}: number of tuning parameters.
   \item{\code{adaptive}}: if \code{TRUE}, use adaptive penalization.
   \item{\code{gamma0}}: power parameter for constructing adaptive weights.
   \item{\code{minLambda}}: multiplicate factor to determine the minimum lambda.
   \item{\code{niter.eea}}: the number of iterations in the iterative exclusive
extraction algorithm.
   \item{\code{df.tol}}: tolerance.
}
}
\value{
A list consisting of
  \item{Y}{Response matrix}
  \item{X}{Predictor matrix}
  \item{U.path}{solution path of U}
  \item{V.path}{solution path of V}
  \item{D.path}{solution path of D}
  \item{ic.path}{solution path of the information criterion}
  \item{U}{estimated left singular matrix}
  \item{V}{estimated right singular matrix}
  \item{D}{estimated singular values}
  \item{rank}{estimated cointegration rank}
  \item{PI}{estimated matrix \eqn{\Pi}}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
