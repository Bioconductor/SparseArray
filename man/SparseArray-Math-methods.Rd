\name{SparseArray-Math-methods}

\alias{SparseArray-Math-methods}
\alias{SparseArray_Math_methods}
\alias{SparseArray-Math}
\alias{SparseArray_Math}
\alias{SparseArray-math-methods}
\alias{SparseArray_math_methods}
\alias{SparseArray-math}
\alias{SparseArray_math}
\alias{SparseArray-Math2-methods}
\alias{SparseArray_Math2_methods}
\alias{SparseArray-Math2}
\alias{SparseArray_Math2}
\alias{SparseArray-math2-methods}
\alias{SparseArray_math2_methods}
\alias{SparseArray-math2}
\alias{SparseArray_math2}

\alias{Math,SVT_SparseArray-method}
\alias{Math2,SVT_SparseArray-method}
\alias{round,SVT_SparseArray-method}
\alias{signif,SVT_SparseArray-method}

\title{'Math' and 'Math2' methods for SparseArray objects}

\description{
  \link{SparseArray} derivatives support a \emph{subset} of operations
  from the \code{Math} and \code{Math2} groups.
  See \code{?\link[methods]{S4groupGeneric}} in the \pkg{methods} package
  for more information about the \code{Math} and \code{Math2} group generics.

  IMPORTANT NOTES:
  \itemize{
    \item Only operations from these groups that preserve sparsity are
          supported. For example, \code{sqrt()}, \code{trunc()},
          \code{log1p()}, and \code{sin()} are supported, but \code{cumsum()},
          \code{log()}, \code{cos()}, or \code{gamma()} are not.
    \item Only \link{SVT_SparseArray} objects are supported at the moment.
          Support for \link{COO_SparseArray} objects might be added in the
          future.
    \item \code{Math} and \code{Math2} operations only support
          \link{SVT_SparseArray} objects of \code{type()} \code{"double"}
          at the moment.
  }
}

\value{
  A \link{SparseArray} derivative of the same dimensions as the input
  object.
}

\seealso{
  \itemize{
    \item \code{\link[methods]{S4groupGeneric}} in the \pkg{methods} package.

    \item \link{SparseArray} objects.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
m <- matrix(0, nrow=15, ncol=6)
m[c(2, 6, 12:17, 22:33, 55, 59:62, 90)] <-
    c(runif(22)*1e4, Inf, -Inf, NA, NaN)
svt <- SparseArray(m)

svt2 <- trunc(sqrt(svt))
svt2

## Sanity check:
m2 <- suppressWarnings(trunc(sqrt(m)))
stopifnot(identical(as.matrix(svt2), m2))
}
\keyword{array}
\keyword{methods}
\keyword{arith}
