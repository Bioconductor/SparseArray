\name{SparseArray-abind}

\alias{SparseArray-abind}
\alias{SparseArray_abind}
\alias{SparseArray-combine}
\alias{SparseArray_combine}
\alias{abind,SparseArray-method}

\title{Bind SparseArray objects with abind()}

\description{
  Like ordinary matrices and arrays in base R, \link{SparseMatrix} derivatives
  can be combined by rows or columns, with \code{rbind()} or \code{cbind()},
  and multidimensional \link{SparseArray} derivatives can be bound along
  any dimension with \code{\link[S4Arrays]{abind}()}.

  Note that \code{\link[S4Arrays]{arbind}()} can also be used to combine
  multidimensional array-like objects along their first dimension, and
  \code{\link[S4Arrays]{acbind}()} can be used to combine them along their
  second dimension.
}

\seealso{
  \itemize{
    \item \code{\link[base]{cbind}} in base R.

    \item \code{\link[S4Arrays]{abind}} in the \pkg{S4Arrays} package.

    \item \link{SparseArray} objects.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## Binding SparseMatrix objects with rbind() or cbind()
## ---------------------------------------------------------------------

m1 <- matrix(1:8, ncol=4,
             dimnames=list(NULL, paste0("M1y", 1:4)))
m2 <- matrix(101:120, ncol=4,
             dimnames=list(paste0("M2x", 1:5), NULL))
sm1 <- SparseArray(m1)
sm2 <- SparseArray(m2)

rbind(sm1, sm2)

m3 <- matrix(0.5 - 1:28, ncol=4,
             dimnames=list(paste0("M3x", 1:7), paste0("M3y", 1:4)))
sm3 <- SparseArray(m3)

rbind(sm2, sm3, sm1)

cbind(sm3, rbind(sm2, sm1))

## ---------------------------------------------------------------------
## Binding 3D SparseArray objects with abind()
## ---------------------------------------------------------------------

a1 <- array(1:105, dim=c(5, 7, 3),
            dimnames=list(NULL, paste0("A1y", 1:7), NULL))
a2 <- array(1001:1105, dim=c(5, 7, 3),
            dimnames=list(paste0("A2x", 1:5), paste0("A2y", 1:7), NULL))
sa1 <- SparseArray(a1)
sa2 <- SparseArray(a2)

abind(sa1, sa2)               # same as 'abind(sa1, sa2, along=3)'
abind(sa1, sa2, rev.along=0)  # same as 'abind(sa1, sa2, along=4)'

a3 <- array(1:60, dim=c(3, 5, 4),
            dimnames=list(NULL, paste0("A3y", 1:5), NULL))
a4 <- array(101:240, dim=c(7, 5, 4),
            dimnames=list(paste0("A4x", 1:7), paste0("A4y", 1:5), NULL))
sa3 <- SparseArray(a3)
sa4 <- SparseArray(a4)

arbind(sa3, sa4)  # same as 'abind(sa3, sa4, along=1)'

## ---------------------------------------------------------------------
## Sanity checks
## ---------------------------------------------------------------------

stopifnot(
    identical(as.matrix(rbind(sm1, sm2)),
              rbind(m1, m2)),
    identical(as.matrix(rbind(sm2, sm3, sm1)),
              rbind(m2, m3, m1)),
    identical(as.matrix(cbind(sm3, rbind(sm2, sm1))),
              cbind(m3, rbind(m2, m1))),
    identical(as.array(abind(sa1, sa2)),
              abind(a1, a2)),
    identical(as.array(abind(sa1, sa2, rev.along=0)),
              abind(a1, a2, rev.along=0)),
    identical(as.array(arbind(sa3, sa4)),
              arbind(a3, a4))
)
}
\keyword{array}
\keyword{methods}
\keyword{manip}
