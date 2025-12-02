\name{NaArray-subassignment}

\alias{NaArray-subassignment}
\alias{NaArray_subassignment}

\title{NaArray subassignment}

\description{
  Like ordinary arrays in base R, \link{NaArray} objects support
  subassignment via the \code{[<-} operator.
}

\seealso{
  \itemize{
    \item \code{\link[base]{[<-}} in base R.

    \item \link{NaArray} objects.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
a <- array(NA, dim=5:3)
a[c(1:2, 8, 10, 15:17, 20, 24, 40, 56:60)] <- 0.01 * (1:15)
naa <- NaArray(a)
naa

naa[ , 1:2, 3] <- - naa[ , 1:2 ,1]

naa[5:3, c(4,2,4), 1:2] <- c(0, NA, 0)

## Sanity checks:
a[ , 1:2, 3] <- -a[ , 1:2 ,1]
a[5:3, c(4,2,4), 1:2] <- c(0, NA, 0)
stopifnot(identical(as.array(naa), a), identical(naa, NaArray(a)))
}
\keyword{array}
\keyword{methods}
