\name{SparseArray-combine}

\alias{SparseArray-combine}
\alias{SparseArray_combine}

\alias{arbind,COO_SparseArray-method}
\alias{acbind,COO_SparseArray-method}
\alias{arbind,SVT_SparseArray-method}
\alias{acbind,SVT_SparseArray-method}

\alias{rbind,SparseArray-method}
\alias{cbind,SparseArray-method}
\alias{bindROWS,SparseArray-method}

\title{Combine SparseArray objects by rows or columns}

\description{
  Like ordinary matrices in base R, \link{SparseMatrix} derivatives can be
  combined by rows or columns, with \code{rbind()} or \code{cbind()}.

  The \pkg{SparseArray} package extends these operations to \link{SparseArray}
  derivatives with an arbitrary number of dimensions, as long as the objects
  to combine have the same number of dimensions and their dimensions are
  compatible with the binding operation. In such case, \code{rbind()} or
  \code{arbind()} can be used to combine the objects along their first
  dimension, and \code{cbind()} or \code{acbind()} can be used to combine
  them along their second dimension.
}

\seealso{
  \itemize{
    \item \code{\link[base]{cbind}} in base R.

    \item \code{\link[S4Arrays]{arbind}} in the \pkg{S4Arrays} package.

    \item \link{SparseArray} objects.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## Combining SparseMatrix objects:

m1 <- matrix(1:15, nrow=3, ncol=5,
             dimnames=list(NULL, paste0("M1y", 1:5)))
m2 <- matrix(101:135, nrow=7, ncol=5,
             dimnames=list(paste0("M2x", 1:7), paste0("M2y", 1:5)))
sm1 <- SparseArray(m1)
sm2 <- SparseArray(m2)
sm1
sm2

rbind(sm1, sm2)

## Combining SparseArray objects with 3 dimensions:

a1 <- array(1:60, c(3, 5, 4),
            dimnames=list(NULL, paste0("A1y", 1:5), NULL))
a2 <- array(101:240, c(7, 5, 4),
            dimnames=list(paste0("A2x", 1:7), paste0("A2y", 1:5), NULL))
sa1 <- SparseArray(a1)
sa2 <- SparseArray(a2)
sa1
sa2

rbind(sa1, sa2)  # same as arbind(sa1, sa2)

## Sanity checks:
sm3 <- rbind(sm1, sm2)
m3 <- rbind(m1, m2)
stopifnot(identical(as.array(sm3), m3), identical(sm3, SparseArray(m3)))
sa3 <- rbind(sa1, sa2)
a3 <- arbind(a1, a2)
stopifnot(identical(as.array(sa3), a3), identical(sa3, SparseArray(a3)))
}
\keyword{methods}
\keyword{array}
