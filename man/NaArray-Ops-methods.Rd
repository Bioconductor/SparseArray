\name{NaArray-Ops-methods}

\alias{NaArray-Ops-methods}
\alias{NaArray_Ops_methods}
\alias{NaArray-Ops}
\alias{NaArray_Ops}
\alias{NaArray-Compare}
\alias{NaArray_Compare}

\alias{Compare,NaArray,vector-method}
\alias{Compare,vector,NaArray-method}

\title{'Ops' methods for NaArray objects}

\description{
  EXPERIMENTAL!!!

  \link{NaArray} derivatives support operations from the \code{Arith},
  \code{Compare}, and \code{Logic} groups, with some restrictions.
  All together, these groups are referred to as the \code{Ops} group.
  See \code{?\link[methods]{S4groupGeneric}} in the \pkg{methods} package
  for more information about the \code{Ops} group generic.

  WORK IN PROGRESS: Only operations from the \code{Compare} group are
  supported on \link{NaArray} objects at the moment. Also they're only
  supported between an \link{NaArray} object \code{naa} and a single
  value \code{y}: \preformatted{    naa op y
    y op naa}
}

\details{
  COMING SOON...
}

\value{
  An \link{NaArray} derivative of the same dimensions as the input
  object(s).
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
\keyword{arith}
