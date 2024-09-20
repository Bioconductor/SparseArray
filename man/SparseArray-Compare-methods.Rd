\name{SparseArray-Compare-methods}

\alias{SparseArray-Compare-methods}
\alias{SparseArray_Compare-methods}
\alias{SparseArray-Compare}
\alias{SparseArray_Compare}
\alias{SparseArray-compare-methods}
\alias{SparseArray_compare-methods}
\alias{SparseArray-compare}
\alias{SparseArray_compare}

\alias{Compare,SVT_SparseArray,vector-method}
\alias{Compare,vector,SVT_SparseArray-method}
\alias{Compare,SVT_SparseArray,SVT_SparseArray-method}
\alias{Compare,SVT_SparseArray,array-method}
\alias{Compare,array,SVT_SparseArray-method}

\title{'Compare' operations on SparseArray objects}

\description{
  \link{SparseArray} derivatives support operations from the \code{Compare}
  group, with some restrictions.
  See \code{?\link[methods]{S4groupGeneric}} in the \pkg{methods} package
  for more information about the \code{Compare} group generic.

  IMPORTANT NOTE: Only \link{SVT_SparseArray} objects are supported at the
  moment. Support for \link{COO_SparseArray} objects might be added in the
  future.
}

\details{
  Two forms of 'Compare' operations are supported:
  \enumerate{
    \item Between an \link{SVT_SparseArray} object \code{svt} and a single
          value \code{y}: \preformatted{    svt op y
    y op svt}

          All operations from the \code{Compare} group support this form,
          with single value \code{y} either on the left or the right.
          However, there are some operation-dependent restrictions on the
          value of \code{y}.

    \item Between two \link{SVT_SparseArray} objects \code{svt1} and
          \code{svt2} of same dimensions (a.k.a. \emph{conformable arrays}):
          \preformatted{    svt1 op svt2}

          The \code{Compare} operations that support this form are:
          \code{!=}, \code{<}, \code{>}.
  }
}

\value{
  A \link{SparseArray} derivative of \code{type()} \code{"logical"} and
  same dimensions as the input object(s).
}

\seealso{
  \itemize{
    \item \code{\link[methods]{S4groupGeneric}} in the \pkg{methods} package.

    \item \link{SparseArray} objects.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
svt1 <- SVT_SparseArray(dim=c(15, 6), type="double")
svt1[c(2, 6, 12:17, 22:33, 55, 59:62, 90)] <- runif(26)
svt1

svt1 >= 0.2
svt1 != 0

## Sanity checks:
m1 <- as.matrix(svt1)
stopifnot(
  identical(as.matrix(svt1 >= 0.2), m1 >= 0.2),
  identical(as.matrix(svt1 != 0), m1 != 0)
)
}
\keyword{array}
\keyword{methods}
\keyword{algebra}
