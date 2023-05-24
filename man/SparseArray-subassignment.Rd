\name{SparseArray-subassignment}

\alias{SparseArray-subassignment}
\alias{SparseArray_subassignment}

\alias{[<-,SVT_SparseArray,ANY,ANY,ANY-method}

\title{SparseArray subassignment}

\description{
  Like ordinary arrays in base R, \link{SparseArray} derivatives support
  subassignment via the \code{[<-} operator.
}

\seealso{
  \itemize{
    \item \code{\link[base]{[<-}} in base R.

    \item \link{SparseArray} objects.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
a <- array(0L, dim=5:3)
a[c(1:2, 8, 10, 15:17, 20, 24, 40, 56:60)] <- (1:15)*10L
svt <- SparseArray(a)
svt

svt[5:3, c(4,2,4), 2:3] <- -99L

## Sanity checks:
a[5:3, c(4,2,4), 2:3] <- -99L
stopifnot(identical(as.array(svt), a), identical(svt, SparseArray(a)))
}
\keyword{array}
\keyword{methods}
