\name{NaArray-Compare-methods}

\alias{NaArray-Compare-methods}
\alias{NaArray_Compare-methods}
\alias{NaArray-Compare}
\alias{NaArray_Compare}
\alias{NaArray-compare-methods}
\alias{NaArray_compare-methods}
\alias{NaArray-compare}
\alias{NaArray_compare}

\alias{Compare,NaArray,vector-method}
\alias{Compare,vector,NaArray-method}
\alias{Compare,NaArray,NaArray-method}
\alias{Compare,NaArray,array-method}
\alias{Compare,array,NaArray-method}
\alias{Compare,NaArray,SVT_SparseArray-method}
\alias{Compare,SVT_SparseArray,NaArray-method}

\title{'Compare' operations on NaArray objects}

\description{
  EXPERIMENTAL!!!

  \link{NaArray} objects support all operations from the \code{Compare} group.
  See \code{?\link[methods]{S4groupGeneric}} in the \pkg{methods} package
  for more information about the \code{Compare} group generic.
}

\details{
  Three forms of 'Compare' operations involving \link{NaArray} objects are
  supported:
  \enumerate{
    \item Between an \link{NaArray} object \code{naa} and a single
          value \code{y}: \preformatted{    naa op y
    y op naa}

    \item Between two \link{NaArray} objects \code{naa1} and \code{naa2}
          of same dimensions (a.k.a. \emph{conformable arrays}):
          \preformatted{    naa1 op naa2}

    \item Between an \link{NaArray} object \code{naa} and an
          \link{SVT_SparseArray} object \code{svt} of same dimensions
          (a.k.a. \emph{conformable arrays}):
          \preformatted{    naa op svt
    svt op naa}
  }
}

\value{
  An \link{NaArray} object of \code{type()} \code{"logical"} and same
  dimensions as the input object(s).
}

\seealso{
  \itemize{
    \item \code{\link[methods]{S4groupGeneric}} in the \pkg{methods} package.

    \item \link{NaArray} objects.

    \item \link{SVT_SparseArray} objects.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
nam1 <- NaArray(dim=c(15, 6), type="double")
nam1[c(2, 6, 12:17, 22:33, 55, 59:62, 90)] <- runif(26)
nam1

nam1 >= 0.2
nam1 != 0

nam2 <- NaArray(dim=dim(nam1), type="integer")
nam2[cbind(1:15, 2)] <- 100:114
nam2[cbind(1:15, 5)] <- -(114:100)
nam2

nam1 < nam2

## Sanity checks:
m1 <- as.matrix(nam1)
m2 <- as.matrix(nam2)
stopifnot(
  identical(as.matrix(nam1 >= 0.2), m1 >= 0.2),
  identical(as.matrix(nam1 != 0), m1 != 0),
  identical(as.matrix(nam1 < nam2), m1 < m2)
)
}
\keyword{array}
\keyword{methods}
\keyword{algebra}
