\name{NaArray-Math-methods}

\alias{NaArray-Math-methods}
\alias{NaArray_Math_methods}
\alias{NaArray-Math}
\alias{NaArray_Math}
\alias{NaArray-math-methods}
\alias{NaArray_math_methods}
\alias{NaArray-math}
\alias{NaArray_math}
\alias{NaArray-Math2-methods}
\alias{NaArray_Math2_methods}
\alias{NaArray-Math2}
\alias{NaArray_Math2}
\alias{NaArray-math2-methods}
\alias{NaArray_math2_methods}
\alias{NaArray-math2}
\alias{NaArray_math2}

\alias{Math,NaArray-method}
\alias{Math2,NaArray-method}
\alias{round,NaArray-method}
\alias{signif,NaArray-method}

\title{'Math' and 'Math2' methods for NaArray objects}

\description{
  EXPERIMENTAL!!!

  \link{NaArray} objects support all operations from the \code{Math}
  and \code{Math2} groups with a few exceptions.
  See \code{?\link[methods]{S4groupGeneric}} in the \pkg{methods} package
  for more information about the \code{Math} and \code{Math2} group generics.

  Note that \code{Math} and \code{Math2} operations only support
  \link{NaArray} objects of \code{type()} \code{"double"} at the moment.
  \link{NaArray} objects of \code{type()} \code{"integer"} are not
  supported yet.
}

\value{
  A \link{NaArray} derivative of the same dimensions as the input
  object.
}

\seealso{
  \itemize{
    \item \code{\link[methods]{S4groupGeneric}} in the \pkg{methods} package.

    \item \link{NaArray} objects.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
nam <- NaArray(dim=c(15, 6))
nam[c(2, 6, 12:17, 22:33, 55, 59:62, 90)] <-
               c(runif(22)*1e4, Inf, -Inf, NA, NaN)

log(nam)
exp(nam)
cos(nam)
gamma(nam)

## Sanity checks:
m <- as.matrix(nam)
stopifnot(
  suppressWarnings(identical(as.matrix(log(nam)), log(m))),
  identical(as.matrix(exp(nam)), exp(m)),
  suppressWarnings(identical(as.matrix(cos(nam)), cos(m))),
  suppressWarnings(identical(as.matrix(gamma(nam)), gamma(m)))
)
}
\keyword{array}
\keyword{methods}
\keyword{arith}
