\name{SparseArray-abind}

\alias{SparseArray-abind}
\alias{SparseArray_abind}
\alias{SparseArray-combine}
\alias{SparseArray_combine}

\alias{abind,SparseArray-method}

\alias{rbind,SparseArray-method}
\alias{cbind,SparseArray-method}
\alias{bindROWS,SparseArray-method}

\title{Combine multidimensional SparseArray objects}

\description{
  Like ordinary matrices and arrays in base R, \link{SparseMatrix} derivatives
  can be combined by rows or columns, with \code{rbind()} or \code{cbind()},
  and multidimensional \link{SparseArray} derivatives can be bound along
  any dimension with \code{abind()}.

  Note that \code{arbind()} can also be used to combine the objects along
  their first dimension, and \code{acbind()} can be used to combine them
  along their second dimension.
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
## COMBINING SparseMatrix OBJECTS
## ---------------------------------------------------------------------

m1a <- matrix(1:15, nrow=3, ncol=5,
              dimnames=list(NULL, paste0("M1y", 1:5)))
m1b <- matrix(101:135, nrow=7, ncol=5,
              dimnames=list(paste0("M2x", 1:7), paste0("M2y", 1:5)))
sm1a <- SparseArray(m1a)
sm1b <- SparseArray(m1b)

rbind(sm1a, sm1b)

## ---------------------------------------------------------------------
## COMBINING SparseArray OBJECTS WITH 3 DIMENSIONS
## ---------------------------------------------------------------------

a2a <- array(1:105, dim=c(5, 7, 3),
             dimnames=list(NULL, paste0("A1y", 1:7), NULL))
a2b <- array(1001:1105, dim=c(5, 7, 3),
             dimnames=list(paste0("A2x", 1:5), paste0("A2y", 1:7), NULL))
sa2a <- SparseArray(a2a)
sa2b <- SparseArray(a2b)

abind(sa2a, sa2b)               # same as 'abind(sa2a, sa2b, along=3)'
abind(sa2a, sa2b, rev.along=0)  # same as 'abind(sa2a, sa2b, along=4)'

a3a <- array(1:60, dim=c(3, 5, 4),
             dimnames=list(NULL, paste0("A1y", 1:5), NULL))
a3b <- array(101:240, dim=c(7, 5, 4),
             dimnames=list(paste0("A2x", 1:7), paste0("A2y", 1:5), NULL))
sa3a <- SparseArray(a3a)
sa3b <- SparseArray(a3b)

arbind(sa3a, sa3b)  # same as 'abind(sa3a, sa3b, along=1)'

## ---------------------------------------------------------------------
## Sanity checks
## ---------------------------------------------------------------------

sm1 <- rbind(sm1a, sm1b)
m1 <- rbind(m1a, m1b)
stopifnot(identical(as.array(sm1), m1), identical(sm1, SparseArray(m1)))

sa2 <- abind(sa2a, sa2b)
stopifnot(identical(sa2, abind(sa2a, sa2b, along=3)))
a2 <- abind(a2a, a2b, along=3)
stopifnot(identical(as.array(sa2), a2), identical(sa2, SparseArray(a2)))

sa2 <- abind(sa2a, sa2b, rev.along=0)
stopifnot(identical(sa2, abind(sa2a, sa2b, along=4)))
a2 <- abind(a2a, a2b, along=4)
stopifnot(identical(as.array(sa2), a2), identical(sa2, SparseArray(a2)))

sa3 <- arbind(sa3a, sa3b)
a3 <- arbind(a3a, a3b)
stopifnot(identical(as.array(sa3), a3), identical(sa3, SparseArray(a3)))
}
\keyword{array}
\keyword{methods}
\keyword{manip}
