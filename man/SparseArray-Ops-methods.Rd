\name{SparseArray-Ops-methods}

\alias{SparseArray_Ops_methods}
\alias{SparseArray_Ops}

\alias{+,SparseArray,missing-method}
\alias{-,SparseArray,missing-method}

\alias{Arith,SVT_SparseArray,vector-method}
\alias{Arith,vector,SVT_SparseArray-method}
\alias{Arith,SVT_SparseArray,SVT_SparseArray-method}
\alias{Arith,SVT_SparseArray,array-method}
\alias{Arith,array,SVT_SparseArray-method}

\alias{Compare,SVT_SparseArray,SVT_SparseArray-method}
\alias{Compare,SVT_SparseArray,array-method}
\alias{Compare,array,SVT_SparseArray-method}

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
    \item Operations from the \code{Compare} and \code{Logic} groups are not
          ready yet but will be soon.
    \item Operations involving complex numbers are not supported yet.
  }
}

\details{
  Two forms of operations are supported:
  \enumerate{
    \item Between an \link{SVT_SparseArray} object \code{svt} and a single
          number \code{x}: \preformatted{    svt op x}
          The operations from the \code{Arith} group that support this form
          are: \code{*}, \code{/}, \code{^}, \code{\%\%},\code{\%/\%}.
          Note that, except for \code{*} (for which both \code{svt * x}
          and \code{x * svt} are supported), the single number must be on
          the right e.g. \code{svt ^ 3}.

    \item Between two \link{SVT_SparseArray} objects \code{svt1} and
          \code{svt2} of same dimensions (a.k.a. \emph{conformable arrays}):
          \preformatted{    svt1 op svt2}
          The operations from the \code{Arith} group that support this form
          are: \code{+}, \code{-}, \code{*}.
  }
}

\value{
  An \link{SparseArray} derivative of the same dimensions as the input
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
svt2 <- (svt ^ 1.5 + svt) \%\% 50 - 10 * svt

## Sanity check:
stopifnot(identical(as.matrix(svt2), (m ^ 1.5 + m) \%\% 50 - 10 * m))
}
\keyword{array}
\keyword{methods}
\keyword{algebra}
\keyword{arith}
