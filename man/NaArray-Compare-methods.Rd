\name{NaArray-Compare-methods}

\alias{NaArray-Compare-methods}
\alias{NaArray_Compare-methods}
\alias{NaArray-Compare}
\alias{NaArray_Compare}

\alias{Compare,NaArray,vector-method}
\alias{Compare,vector,NaArray-method}

\title{'Compare' operations on NaArray objects}

\description{
  EXPERIMENTAL!!!

  \link{NaArray} objects support operations from the \code{Compare} group.
  See \code{?\link[methods]{S4groupGeneric}} in the \pkg{methods} package
  for more information about the \code{Compare} group generic.
}

\details{
  COMING SOON...
}

\value{
  An \link{NaArray} object of the same dimensions as the input object(s).
}

\seealso{
  \itemize{
    \item \code{\link[methods]{S4groupGeneric}} in the \pkg{methods} package.

    \item \link{NaArray} objects.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
nam <- NaArray(dim=c(15, 6), type="double")
nam[c(2, 6, 12:17, 22:33, 55, 59:62, 90)] <- runif(26)
nam

nam >= 0.2
nam != 0

## Sanity check:
m <- as.matrix(nam)
stopifnot(
  identical(as.matrix(nam >= 0.2), m >= 0.2),
  identical(as.matrix(nam != 0), m != 0)
)
}
\keyword{array}
\keyword{methods}
\keyword{algebra}
