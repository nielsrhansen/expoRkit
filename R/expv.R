### R interface for calling Expokit Fortran routines for computing the
### matrix exponential. This file contains S4-methods for different matrix
### classes from the Matrix package and SparseM package together with a
### method for ordinary matrices. 
###
###
###     Copyright (C) 2012 Niels Richard Hansen.
###
### This program is free software; you can redistribute it and/or modify it
### under the terms of the GNU General Public License as published by the
### Free Software Foundation; either version 2, or (at your option) any
### later version.
###
### These functions are distributed in the hope that they will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, a copy is available at
### http://www.r-project.org/Licenses/

setGeneric("expv", function(x, v, ...) standardGeneric("expv"))

setMethod("expv", signature("matrix", "vector"),
          function(x, v, t = 1.0, u = NULL, ...) {
            if (is.numeric(x)) {
            ## For a numeric 'matrix' argument coerce the argument into a
            ## sparse matrix in CCS format.
            w <- callGeneric(Matrix(x, sparse = TRUE), v = v, t = t, u = u, ...)
          } else if (is.complex(x)) {
            ## For a complex 'matrix' the coercion to a sparse matrix
            ## is done manually as the Matrix package does not support
            ## complex matrices (yet).
            pattern <- x != 0
            sparsePattern <- as(pattern, 'Matrix')
            w <- matrix(0, length(v), length(t))
            if (length(t) > 1)
              t <- c(t[1], diff(t))
            for(k in seq_along(t)) {
              w[, k] <- Rexpv(a = x[pattern],
                              ia = sparsePattern@i + 1,
                              ja = sparsePattern@p + 1,
                              n = nrow(x),
                              v = v,
                              u = u,
                              t = t[k],
                              Markov = FALSE,
                              storage = 'CCS',
                              ...)
              v <- w[, k]
            }
          } else {
            stop("Argument 'x' must be a 'numeric' or 'complex' matrix.")
          }
            w
          }
          )

setMethod("expv", signature("CsparseMatrix", "vector"),
          function(x, v, t = 1.0, u = NULL, Markov = FALSE,
                   transpose = Markov, ...) {
            if (transpose)
              x <- t(x)
            w <- matrix(0, length(v), length(t))
            if (length(t) > 1)
              t <- c(t[1], diff(t))
            for(k in seq_along(t)) {
              w[, k] <- Rexpv(a = x@x,
                              ia = x@i + 1,
                              ja = x@p + 1,
                              n = nrow(x),
                              v = v,
                              u = u,
                              t = t[k],
                              Markov = Markov,
                              storage = 'CCS',
                              ...)
              v <- w[, k]
            }
            w
          }
          )

setMethod("expv", signature("TsparseMatrix", "vector"),
          function(x, v, t = 1.0, u = NULL, Markov = FALSE,
                   transpose = Markov, ...) {
            if (transpose)
              x <- t(x)
            w <- matrix(0, length(v), length(t))
            if (length(t) > 1)
              t <- c(t[1], diff(t))
            for(k in seq_along(t)) {
              w[, k] <- Rexpv(a = x@x,
                              ia = x@i + 1,
                              ja = x@j + 1,
                              n = nrow(x),
                              v = v,
                              u = u,
                              t = t[k],
                              Markov = Markov,
                              storage = 'COO',
                              ...)
              v <- w[, k]
            }
            w
          }
          )

setMethod("expv", signature("matrix.csc", "vector"),
          function(x, v, t = 1.0, u = NULL, Markov = FALSE,
                   transpose = Markov, ...) {
            if (transpose)
              x <- t(x)
            w <- matrix(0, length(v), length(t))
            if (length(t) > 1)
              t <- c(t[1], diff(t))
            for(k in seq_along(t)) {
              w[, k] <- Rexpv(a = x@ra,
                              ia = x@ja,
                              ja = x@ia,
                              n = nrow(x),
                              v = v,
                              u = u,
                              t = t[k],
                              Markov = Markov,
                              storage = 'CCS',
                              ...)
              v <- w[, k]
            }
            w
          }
          )

setMethod("expv", signature("matrix.csr", "vector"),
          function(x, v, t = 1.0, u = NULL, Markov = FALSE,
                   transpose = Markov, ...) {
            if (transpose)
              x <- t(x)
            w <- matrix(0, length(v), length(t))
            if (length(t) > 1)
              t <- c(t[1], diff(t))
            for(k in seq_along(t)) {
              w[, k] <- Rexpv(a = x@ra,
                              ia = x@ia,
                              ja = x@ja,
                              n = nrow(x),
                              v = v,
                              u = u,
                              t = t[k],
                              Markov = Markov,
                              storage = 'CRS',
                              ...)
              v <- w[, k]
            }
            w
          }
          )

setMethod("expv", signature("matrix.coo", "vector"),
          function(x, v, t = 1.0, u = NULL, Markov = FALSE,
                   transpose = Markov, ...) {
            if (transpose)
              x <- t(x)
            w <- matrix(0, length(v), length(t))
            if (length(t) > 1)
              t <- c(t[1], diff(t))
            for(k in seq_along(t)) {
              w[, k] <- Rexpv(a = x@ra,
                              ia = x@ia,
                              ja = x@ja,
                              n = nrow(x),
                              v = v,
                              u = u,
                              t = t[k],
                              Markov = Markov,
                              storage = 'COO'
                              ...)
              v <- w[, k]
            }
            w
          }
          )




            
