\name{expv}
\alias{expv}
\alias{expv,CsparseMatrix,vector-method}
\alias{expv,TsparseMatrix,vector-method}
\alias{expv,dgCMatrix,vector-method}
\alias{expv,matrix,vector-method}
\alias{expv,matrix.coo,vector-method}
\alias{expv,matrix.csc,vector-method}
\alias{expv,matrix.csr,vector-method}
\alias{expv-methods}
\title{Matrix exponentiation and more.}
\usage{
  \S4method{expv}{matrix,vector}(x, v, t = 1, u = NULL,
    Markov = FALSE, transpose = Markov, ...)

  \S4method{expv}{CsparseMatrix,vector}(x, v, t = 1,
    u = NULL, Markov = FALSE, transpose = Markov, ...)

  \S4method{expv}{dgCMatrix,vector}(x, v, t = 1, u = NULL,
    Markov = FALSE, transpose = Markov, ...)

  \S4method{expv}{TsparseMatrix,vector}(x, v, t = 1,
    u = NULL, Markov = FALSE, transpose = Markov, ...)

  \S4method{expv}{matrix.csc,vector}(x, v, t = 1, u = NULL,
    Markov = FALSE, transpose = Markov, ...)

  \S4method{expv}{matrix.csr,vector}(x, v, t = 1, u = NULL,
    Markov = FALSE, transpose = Markov, ...)

  \S4method{expv}{matrix.coo,vector}(x, v, t = 1, u = NULL,
    Markov = FALSE, transpose = Markov, ...)
}
\arguments{
  \item{x}{a matrix.}

  \item{v}{a \code{numeric} vector. The initial value.}

  \item{t}{a \code{numeric} vector of time points at which
  the solution is computed.}

  \item{u}{a \code{numeric} vector. Default \code{NULL}.}

  \item{Markov}{\code{logical}. If \code{TRUE} the matrix
  is taken to be an rate matrix and steps are taken to
  ensure that the computed result is a probability vector.
  Default \code{FALSE}.}

  \item{transpose}{\code{logical}. If \code{TRUE} transpose
  the matrix before the solution is computed. Default
  equals \code{Markov}.}

  \item{...}{other arguments passed to
  \code{\link{Rexpv}}.}
}
\value{
  An \eqn{n} by \eqn{k} matrix with \eqn{k} the length of
  \code{t} and \eqn{n} the dimension of the matrix
  \code{x}. The \eqn{i}'th column contains the solution of
  the ODE at time point \eqn{i}.
}
\description{
  Methods for computing the solution of the ODE \deqn{w'(t)
  = x w(t) + u} with initial condition \eqn{w(0) = v} at
  one or more time points.
}
\details{
  Analytically the solution is given as \deqn{w(t) =
  \exp(tx)v + t \phi(tx)u}{w(t) = exp(tx)v + tphi(tx)u}
  with \eqn{\phi(z) = (\exp(z)-1)/z}{phi(z) =
  (exp(z)-1)/z}. For large matrices \eqn{x} the computation
  of the full matrices \eqn{\exp(tx)}{exp(tx)} and
  \eqn{\phi(tx)}{phi(tx)} is slow. An alternative is to
  compute the solution \eqn{w} directly for a given initial
  value \eqn{v} using Krylov subspace methods. This is, in
  particular, efficient for large sparse matrices.

  Note that if \eqn{Q} is a rate matrix for a homogeneous
  continuous time Markov process (non-negative
  off-diagonals and row-sums 0) and \eqn{v} is a
  probability vector (the initial distribution at time
  point 0), then the distribution at time point \eqn{t}
  solves \deqn{w'(t) = Q^T w(t).} In this case we want to
  take \eqn{x} to be \eqn{Q^T} to get the desired solution.

  The solution is computed using the Fortran package
  Expokit. Methods are available for matrix classes
  implemented in the Matrix package as well as the SparseM
  package. The implementation avoids the computation of the
  full matrix exponential of \eqn{tx} and the approach is
  advantageous when we want to compute \eqn{w(t)} for one
  or a few initial values \eqn{v}. The full matrix
  exponential should \emph{not} be computed this way by
  looping over \eqn{n} different initial values.

  Though there is a method implemented for ordinary (dense)
  matrices, such a matrix is simply coerced into a
  \code{CsparseMatrix} before the solution is computed. It
  is recommended that large sparse matrices are stored and
  handled as such, e.g. using the classes and methods from
  the Matrix package. Dense intermediates should be
  avoided.

  The \code{x} matrix is allowed to be a dense complex
  matrix in which case \code{v} and \code{u} are also
  allowed to be complex.
}
\examples{
### Small 3 by 3 example.
x <- matrix(c(-1, 0, 1,
              0, -2, -3,
              1, 0, 0),
            3, 3)
v <- c(1, 1, 1)

require(Matrix)
require(SparseM)

w <- cbind(padm(x) \%*\% v,
           expv(x, v),
           expv(Matrix(x, sparse = TRUE), v),
           expv(as.matrix.coo(x), v),
           expv(as.matrix.csr(x), v),
           expv(as.matrix.csc(x), v)
           )

stopifnot(all.equal(w[, 1], w[, 2]),
          all.equal(w[, 1], w[, 3]),
          all.equal(w[, 1], w[, 4]),
          all.equal(w[, 1], w[, 5]),
          all.equal(w[, 1], w[, 6]))

u <- c(2, 0, 1)
ex <- padm(x)
w <- cbind(ex \%*\% v + (ex - diag(1, 3)) \%*\% solve(x, u),
           expv(x, v, u = u),
           expv(Matrix(x, sparse = TRUE), v, u = u),
           expv(as.matrix.coo(x), v, u = u),
           expv(as.matrix.csr(x), v, u = u),
           expv(as.matrix.csc(x), v, u = u)
           )

stopifnot(all.equal(w[, 1], w[, 2]),
          all.equal(w[, 1], w[, 3]),
          all.equal(w[, 1], w[, 4]),
          all.equal(w[, 1], w[, 5]),
          all.equal(w[, 1], w[, 6]))

############################################################
### Linear birth-death Markov process with immigration
############################################################

alpha <- 2.1  ## Death rate per individual
beta <- 2     ## Birth rate per individual
delta <- 20   ## Immigration rate

n <- 500L     ## state space {0, ..., n-1}
i <- seq(1, n)
rates <- c(alpha * i[-1], ## subdiagonal
           -(c(0, alpha * i[-1]) +
             c(delta, beta * i[-c(1,n)] + delta, 0)), ## diagonal
           c(delta, beta * i[-c(1, n)] + delta)) ## superdiagonal
j <- c(i[-n], i, i[-1])
i <- c(i[-1], i, i[-n])

## Sparse rate matrix constructed without dense intermediate
Q <- sparseMatrix(i = i, j = j, x = rates, dims = c(n, n))

## Evolution of uniform initial distribution
p0 <- rep(1, n)/n
time <- seq(0, 10, 0.2)
Pt <- expv(Q, p0, t = time, Markov = TRUE)

\dontrun{matplot(1:n, Pt, type = "l")
image(time, 0:(n-1), -t(Pt), col = terrain.colors(100))
}
}
\author{
  Niels Richard Hansen \email{Niels.R.Hansen@math.ku.dk}
}
\seealso{
  \code{\link{Rexpv}}, \code{\link[Matrix]{expm}},
  \code{\link[expm]{expm}}
}

