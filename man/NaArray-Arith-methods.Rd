\name{NaArray-Arith-methods}

\alias{NaArray-Arith-methods}
\alias{NaArray_Arith-methods}
\alias{NaArray-Arith}
\alias{NaArray_Arith}

\alias{+,NaArray,missing-method}
\alias{-,NaArray,missing-method}

\alias{Arith,NaArray,vector-method}
\alias{Arith,vector,NaArray-method}
\alias{Arith,NaArray,NaArray-method}
\alias{Arith,NaArray,array-method}
\alias{Arith,array,NaArray-method}
\alias{Arith,NaArray,SVT_SparseArray-method}
\alias{Arith,SVT_SparseArray,NaArray-method}

\title{'Arith' operations on NaArray objects}

\description{
  EXPERIMENTAL!!!

  \link{NaArray} objects support all operations from the \code{Arith} group.
  See \code{?\link[methods]{S4groupGeneric}} in the \pkg{methods} package
  for more information about the \code{Arith} group generic.

  Note that \link{NaArray} of \code{type()} \code{"complex"} don't support
  \code{Arith} operations at the moment.
}

\details{
  Three forms of 'Arith' operations involving \link{NaArray} objects are
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
  An \link{NaArray} object of the same dimensions as the input object(s).
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
nam1 <- NaArray(dim=c(15, 6), type="integer")
nam1[cbind(1:15, 2)] <- 100:114
nam1[cbind(1:15, 5)] <- -(114:100)
nam1

nam1 * -0.01

nam1 * 10   # result is of type "double"
nam1 * 10L  # result is of type "integer"

nam1 / 10L

nam1 ^ 3.5

nam1 \%\% 5L
nam1 \%/\% 5L

nam2 <- NaArray(dim=dim(nam1), type="double")
nam2[c(2, 6, 12:17, 22:33, 55, 59:62, 90)] <- runif(26)
nam2

nam2 + nam1
nam2 - nam1
nam2 * nam1
nam2 / nam1
nam2 ^ nam1
nam2 \%\% nam1
nam2 \%/\% nam1

## Sanity checks:
m1 <- as.matrix(nam1)
m2 <- as.matrix(nam2)
stopifnot(
  identical(as.matrix(nam1 * -0.01), m1 * -0.01),
  identical(as.matrix(nam1 * 10), m1 * 10),
  identical(as.matrix(nam1 * 10L), m1 * 10L),
  identical(as.matrix(nam1 / 10L), m1 / 10L),
  identical(as.matrix(nam1 ^ 3.5), m1 ^ 3.5),
  identical(as.matrix(nam1 \%\% 5L), m1 \%\% 5L),
  identical(as.matrix(nam1 \%/\% 5L), m1 \%/\% 5L),
  identical(as.matrix(nam2 + nam1), m2 + m1),
  identical(as.matrix(nam2 - nam1), m2 - m1),
  identical(as.matrix(nam2 * nam1), m2 * m1),
  identical(as.matrix(nam2 / nam1), m2 / m1),
  identical(as.matrix(nam2 ^ nam1), m2 ^ m1),
  identical(as.matrix(nam2 \%\% nam1), m2 \%\% m1),
  identical(as.matrix(nam2 \%/\% nam1), m2 \%/\% m1)
)
}
\keyword{array}
\keyword{methods}
\keyword{algebra}
\keyword{arith}
