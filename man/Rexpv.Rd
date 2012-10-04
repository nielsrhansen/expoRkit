\name{Rexpv}
\alias{Rexpv}
\title{Expokit __EXPV and __PHIV wrapper.}
\usage{
  Rexpv(a, ia, ja, n, v, t = 1, storage = "CCS", u = NULL,
    anorm = max(abs(a)), Markov = FALSE, m = 30L,
    tol = 1e-10, itrace = 0L)
}
\arguments{
  \item{a}{\code{numeric} or \code{complex} non-zero
  entries in the \eqn{x}-matrix.}

  \item{ia}{\code{integer} index/pointer. Precise meaning
  depends on storage format.}

  \item{ja}{\code{integer} index/pointer. Precise meaning
  depends on storage format.}

  \item{n}{dimension of the (square) matrix.}

  \item{v}{\code{numeric} or \code{complex} vector.}

  \item{t}{time. Default \code{1}.}

  \item{storage}{\code{character}, one of \code{'CCS'}
  (Compressed Column Storage), \code{'CRS'} (Compressed Row
  Storage) or \code{'COO'} (COOrdinate list). Default
  \code{'CCS'}.}

  \item{u}{\code{numeric} or \code{complex} vector. Default
  \code{NULL}.}

  \item{anorm}{A norm of the matrix. Default is the
  sup-norm.}

  \item{Markov}{\code{logical}, if \code{TRUE} the
  (transposed) matrix is taken to be an intensity matrix
  and steps are taken to ensure that the computed result is
  a probability vector. Default \code{FALSE}.}

  \item{m}{\code{integer}, the maximum size for the Krylov
  basis.}

  \item{tol}{\code{numeric}. Default \code{1e-10}.}

  \item{itrace}{\code{integer}, 0 (default) means no trace
  information from Expokit, and 1 means trace information
  printed from Expokit.}
}
\value{
  The solution, \eqn{w}, of the ODE as a \code{numeric} or
  \code{complex} vector of length \eqn{n}.
}
\description{
  R wrapper of the Expokit Fortran subroutines __EXPV and
  __PHIV for sparse matrix exponentiation. In general,
  these routines compute the solution at time point \eqn{t}
  of the ODE \deqn{w'(t) = x w(t) + u} with initial
  condition \eqn{w(0) = v}.
}
\details{
  The \code{Rexpv} function is the low level wrapper of the
  Fortran subroutines in the Expokit package. It is not
  intended to be used directly but rather via the
  \code{\link{expv}} methods. In the call the correct
  storage format in terms of the vectors \code{a},
  \code{ia} and \code{ja} has to be specified via the
  \code{storage} argument. For CCS, \code{ia} contains
  1-based row numbers of non-zero elements and \code{ja}
  contains 1-based pointers to the initial element for each
  column. For CRS, \code{ia} contains 1-based column
  numbers of non-zero elements and \code{ja} are 1-based
  pointers to the initial element for each row. For COO,
  \code{ia} and \code{ja} contain the 1-based column and
  row numbers, respectively, for the non-zero elements.
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
  \code{\link[Matrix]{expm}}, \code{\link[expm]{expm}},
  \code{\link{expv}}
}

