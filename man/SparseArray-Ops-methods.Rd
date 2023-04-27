\name{SparseArray-Ops-methods}

\alias{SparseArray-Ops-methods}
\alias{SparseArray_Ops_methods}
\alias{SparseArray-Ops}
\alias{SparseArray_Ops}
\alias{SparseArray-Arith}
\alias{SparseArray_Arith}
\alias{SparseArray-Compare}
\alias{SparseArray_Compare}
\alias{SparseArray-Logic}
\alias{SparseArray_Logic}

\alias{+,SparseArray,missing-method}
\alias{-,SparseArray,missing-method}

\alias{Arith,SVT_SparseArray,vector-method}
\alias{Arith,vector,SVT_SparseArray-method}
\alias{Arith,SVT_SparseArray,SVT_SparseArray-method}
\alias{Arith,SVT_SparseArray,array-method}
\alias{Arith,array,SVT_SparseArray-method}

\alias{Compare,SVT_SparseArray,vector-method}
\alias{Compare,vector,SVT_SparseArray-method}
\alias{Compare,SVT_SparseArray,SVT_SparseArray-method}
\alias{Compare,SVT_SparseArray,array-method}
\alias{Compare,array,SVT_SparseArray-method}

\alias{Logic,SVT_SparseArray,vector-method}
\alias{Logic,vector,SVT_SparseArray-method}
\alias{Logic,SVT_SparseArray,SVT_SparseArray-method}
\alias{Logic,SVT_SparseArray,array-method}
\alias{Logic,array,SVT_SparseArray-method}

\title{'Ops' methods for SparseArray objects}

\description{
  \link{SparseArray} derivatives support operations from the \code{Arith},
  \code{Compare}, and \code{Logic} groups, with some restrictions.
  All together, these groups are referred to as the \code{Ops} group. See
  \code{?\link[methods]{S4groupGeneric}} in the \pkg{methods} package for
  more information.

  IMPORTANT NOTES:
  \itemize{
    \item Only \link{SVT_SparseArray} objects are supported at the moment.
          \link{COO_SparseArray} might be added in the future.
    \item \code{Arith} operations don't support \link{SVT_SparseArray}
          objects of \code{type()} \code{"complex"} at the moment.
    \item Operations from the \code{Logic} groups are not ready yet but
          will be soon.
  }
}

\details{
  Two forms of operations are supported:
  \enumerate{
    \item Between an \link{SVT_SparseArray} object \code{svt} and a single
          value \code{y}: \preformatted{    svt op y
    y op svt}

          The operations from the \code{Arith} group that support this form
          are: \code{*}, \code{/}, \code{^}, \code{\%\%},\code{\%/\%}.
          Note that, except for \code{*} (for which both \code{svt * y}
          and \code{y * svt} are supported), single value \code{y} must
          be on the right e.g. \code{svt ^ 3}.

          All operations from the \code{Compare} group support this form,
          with single value \code{y} either on the left or the right.
          However, there are some operation-dependent restrictions on the
          value of \code{y}.

    \item Between two \link{SVT_SparseArray} objects \code{svt1} and
          \code{svt2} of same dimensions (a.k.a. \emph{conformable arrays}):
          \preformatted{    svt1 op svt2}

          The operations from the \code{Arith} group that support this form
          are: \code{+}, \code{-}, \code{*}.

          The operations from the \code{Compare} group that support this form
          are: \code{!=}, \code{<}, \code{>}.
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

    \item \link[base]{array} objects in base R.
  }
}

\examples{
m <- matrix(0L, nrow=15, ncol=6)
m[c(2, 6, 12:17, 22:33, 55, 59:62, 90)] <- 101:126
svt <- SparseArray(m)

## Can be 5x or 10x faster than with a dgCMatrix object on a big
## SVT_SparseMatrix object!
svt2 <- (svt^1.5 + svt) %% 100 - 0.2 * svt > 0

## Sanity check:
m2 <- (m^1.5 + m) %% 100 - 0.2 * m > 0
stopifnot(identical(as.matrix(svt2), m2))
}
\keyword{array}
\keyword{methods}
\keyword{algebra}
\keyword{arith}
