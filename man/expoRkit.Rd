\docType{package}
\name{expoRkit}
\alias{expoRkit}
\alias{expoRkit-package}
\title{Expokit in R}
\description{
  R-interface to the Fortran package Expokit for matrix
  exponentiation.
}
\details{
  \tabular{ll}{ Package: \tab expoRkit \cr Type: \tab
  Package \cr Version: \tab 0.9 \cr Date: \tab
  2012-10-08\cr License: \tab GPL (>= 2)\cr LazyLoad: \tab
  yes \cr }

  Expokit is an efficient Fortran implementation for
  computing the matrix exponential, or rather, its action
  on a vector, for large sparse matrices. This can also be
  understood as computing the solution of a system of
  linear ordinary first order differential equations. This
  package provides an R-interface to some of the Fortran
  subroutines from Expokit.
}
\author{
  Niels Richard Hansen \email{Niels.R.Hansen@math.ku.dk}
}
\references{
  Sidje, R. B. (1998) Expokit. Software Package for
  Computing Matrix Exponentials. ACM Trans. Math. Softw.
  24(1), 130-156.
}
\seealso{
  \code{\link{expv}}
}
\keyword{package}

