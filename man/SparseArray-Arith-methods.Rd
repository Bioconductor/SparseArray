\name{SparseArray-Arith-methods}

\alias{SparseArray-Arith-methods}
\alias{SparseArray_Arith-methods}
\alias{SparseArray-Arith}
\alias{SparseArray_Arith}

\alias{+,SparseArray,missing-method}
\alias{-,SparseArray,missing-method}

\alias{Arith,SVT_SparseArray,vector-method}
\alias{Arith,vector,SVT_SparseArray-method}
\alias{Arith,SVT_SparseArray,SVT_SparseArray-method}
\alias{Arith,SVT_SparseArray,array-method}
\alias{Arith,array,SVT_SparseArray-method}

\title{'Arith' operations on SparseArray objects}

\description{
  \link{SparseArray} derivatives support operations from the \code{Arith}
  group, with some restrictions.
  See \code{?\link[methods]{S4groupGeneric}} in the \pkg{methods} package
  for more information about the \code{Arith} group generic.

  IMPORTANT NOTES:
  \itemize{
    \item Only \link{SVT_SparseArray} objects are supported at the moment.
          Support for \link{COO_SparseArray} objects might be added in the
          future.
    \item \link{SVT_SparseArray} of \code{type()} \code{"complex"} don't
          support \code{Arith} operations at the moment.
  }
}

\details{
  Two forms of 'Arith' operations are supported:
  \enumerate{
    \item Between an \link{SVT_SparseArray} object \code{svt} and a single
          value \code{y}: \preformatted{    svt op y
    y op svt}

          The \code{Arith} operations that support this form are: \code{*},
          \code{/}, \code{^}, \code{\%\%},\code{\%/\%}.
          Note that, except for \code{*} (for which both \code{svt * y}
          and \code{y * svt} are supported), single value \code{y} must
          be on the right e.g. \code{svt ^ 3}.

    \item Between two \link{SVT_SparseArray} objects \code{svt1} and
          \code{svt2} of same dimensions (a.k.a. \emph{conformable arrays}):
          \preformatted{    svt1 op svt2}

          The \code{Arith} operations that support this form are: \code{+},
          \code{-}, \code{*}.
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
## ---------------------------------------------------------------------
## Basic examples
## ---------------------------------------------------------------------

svt1 <- SVT_SparseArray(dim=c(15, 6), type="integer")
svt1[cbind(1:15, 2)] <- 100:114
svt1[cbind(1:15, 5)] <- -(114:100)
svt1

svt1 * -0.01

svt1 * 10   # result is of type "double"
svt1 * 10L  # result is of type "integer"

svt1 / 10L

svt1 ^ 3.5

svt1 \%\% 5L
svt1 \%/\% 5L

svt2 <- SVT_SparseArray(dim=dim(svt1), type="double")
svt2[c(2, 6, 12:17, 22:33, 55, 59:62, 90)] <- runif(26)
svt2

svt1 + svt2
svt1 - svt2
svt1 * svt2

svt2 * (0.1 * svt1 - svt2 ^ 2) + svt1 / sum(svt2)

## Sanity checks:
m1 <- as.matrix(svt1)
m2 <- as.matrix(svt2)
stopifnot(
  identical(as.matrix(svt1 * -0.01), m1 * -0.01),
  identical(as.matrix(svt1 * 10), m1 * 10),
  identical(as.matrix(svt1 * 10L), m1 * 10L),
  identical(as.matrix(svt1 / 10L), m1 / 10L),
  identical(as.matrix(svt1 ^ 3.5), m1 ^ 3.5),
  identical(as.matrix(svt1 \%\% 5L), m1 \%\% 5L),
  identical(as.matrix(svt1 \%/\% 5L), m1 \%/\% 5L),
  identical(as.matrix(svt1 + svt2), m1 + m2),
  identical(as.matrix(svt1 - svt2), m1 - m2),
  identical(as.matrix(svt1 * svt2), m1 * m2),
  all.equal(as.matrix(svt2 * (0.1 * svt1 - svt2 ^ 2) + svt1 / sum(svt2)),
            m2 * (0.1 * m1 - m2 ^ 2) + m1 / sum(m2))
)

## ---------------------------------------------------------------------
## An example combining operations from the 'Arith', 'Compare',
## and 'Logic' groups
## ---------------------------------------------------------------------

m3 <- matrix(0L, nrow=15, ncol=6)
m3[c(2, 6, 12:17, 22:33, 55, 59:62, 90)] <- 101:126
svt3 <- SparseArray(m3)

## Can be 5x or 10x faster than with a dgCMatrix object on a big
## SVT_SparseMatrix object!
svt4 <- (svt3^1.5 + svt3) \%\% 100 - 0.2 * svt3 > 0
svt4

## Sanity check:
m4 <- (m3^1.5 + m3) \%\% 100 - 0.2 * m3 > 0
stopifnot(identical(as.matrix(svt4), m4))
}
\keyword{array}
\keyword{methods}
\keyword{algebra}
\keyword{arith}
