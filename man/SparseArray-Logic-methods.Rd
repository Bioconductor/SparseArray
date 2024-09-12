\name{SparseArray-Logic-methods}

\alias{SparseArray-Logic-methods}
\alias{SparseArray_Logic-methods}
\alias{SparseArray-Logic}
\alias{SparseArray_Logic}

\alias{!,SparseArray-method}

\alias{Logic,SVT_SparseArray,vector-method}
\alias{Logic,vector,SVT_SparseArray-method}
\alias{Logic,SVT_SparseArray,SVT_SparseArray-method}
\alias{Logic,SVT_SparseArray,array-method}
\alias{Logic,array,SVT_SparseArray-method}

\title{'Logic' operations on SparseArray objects}

\description{
  \link{SparseArray} derivatives support operations from the \code{Logic}
  group, with some restrictions.
  See \code{?\link[methods]{S4groupGeneric}} in the \pkg{methods} package
  for more information about the \code{Logic} group generic.

  IMPORTANT NOTES:
  \itemize{
    \item Only \link{SVT_SparseArray} objects are supported at the moment.
          Support for \link{COO_SparseArray} objects might be added in the
          future.
    \item In base R, \code{Logic} operations (i.e. \code{&} and \code{|})
          support input of \code{type()} \code{"logical"}, \code{"integer"},
          \code{"double"}, or \code{"complex"}. However, the corresponding
          methods for \link{SVT_SparseArray} objects only support objects
          of \code{type()} \code{"logical"} for now.
  }
}

\value{
  A \link{SparseArray} derivative of the same dimensions as the input
  object(s).
}

\seealso{
  \itemize{
    \item \code{\link[methods]{S4groupGeneric}} in the \pkg{methods} package.

    \item \link{SparseArray} objects.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
svt1 <- svt2 <- SVT_SparseArray(dim=c(15, 6))
svt1[cbind(1:15, 2)] <- c(TRUE, FALSE, NA)
svt1[cbind(1:15, 5)] <- c(TRUE, NA, NA, FALSE, TRUE)
svt2[c(2, 6, 12:17, 22:33, 55, 59:62, 90)] <- c(TRUE, NA)

svt1 & svt2
svt1 | svt2

## Sanity check:
m1 <- as.matrix(svt1)
m2 <- as.matrix(svt2)
stopifnot(
  identical(as.matrix(svt1 & svt2), m1 & m2),
  identical(as.matrix(svt1 | svt2), m1 | m2)
)
}
\keyword{array}
\keyword{methods}
\keyword{algebra}
