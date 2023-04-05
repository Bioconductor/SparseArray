\name{SparseMatrix-mult}

\alias{SparseMatrix-mult}
\alias{SparseMatrix_mult}

\alias{crossprod,SVT_SparseMatrix,missing-method}
\alias{crossprod,SVT_SparseMatrix,SVT_SparseMatrix-method}
\alias{crossprod,SVT_SparseMatrix,ANY-method}
\alias{crossprod,ANY,SVT_SparseMatrix-method}

\alias{tcrossprod,SVT_SparseMatrix,missing-method}
\alias{tcrossprod,SVT_SparseMatrix,SVT_SparseMatrix-method}
\alias{tcrossprod,SVT_SparseMatrix,ANY-method}
\alias{tcrossprod,ANY,SVT_SparseMatrix-method}

\alias{\%*\%,SVT_SparseMatrix,SVT_SparseMatrix-method}
\alias{\%*\%,SVT_SparseMatrix,ANY-method}
\alias{\%*\%,ANY,SVT_SparseMatrix-method}

\title{SparseMatrix multiplication and cross-product}

\description{
  Like ordinary matrices in base R, \link{SparseMatrix} derivatives can
  be multiplied with the \code{\%*\%} operator. They also support
  \code{\link[base]{crossprod}()} and \code{\link[base]{tcrossprod}()}.

  Note that \code{\%*\%}, \code{crossprod()} and \code{tcrossprod()}
  always return an \emph{ordinary} matrix of \code{type()} \code{"double"}
  on \link{SparseMatrix} derivatives.
}

\seealso{
  \itemize{
    \item \code{\link[base]{\%*\%}} and \code{\link[base]{crossprod}} in base R.

    \item \link{SparseMatrix} objects.

    \item \code{S4Arrays::\link[S4Arrays]{type}} in the \pkg{S4Arrays}
          package to get the type of the elements of an array-like object.

    \item \link[base]{matrix} objects in base R.
  }
}

\examples{
m1 <- matrix(0L, nrow=15, ncol=6)
m1[c(2, 6,12:17, 22:33, 55, 59:62, 90)] <- 101:126
svt1 <- as(m1, "SVT_SparseMatrix")

set.seed(333)
svt2 <- poissonSparseMatrix(nrow=6, ncol=7, density=0.2)

svt1 \%*\% svt2
m1   \%*\% svt2

## Unary crossprod() and tcrossprod():
crossprod(svt1)               # same as t(svt1) \%*\% svt1
tcrossprod(svt1)              # same as svt1 \%*\% t(svt1)

## Binary crossprod() and tcrossprod():
crossprod(svt1[1:6, ], svt2)  # same as t(svt1[1:6, ]) \%*\% svt2
tcrossprod(svt1, t(svt2))     # same as svt1 \%*\% svt2

## Sanity checks:
m12 <- m1 \%*\% as.matrix(svt2)
stopifnot(
  identical(svt1 \%*\% svt2, m12),
  identical(m1 \%*\% svt2, m12),
  identical(crossprod(svt1), t(svt1) \%*\% svt1),
  identical(tcrossprod(svt1), svt1 \%*\% t(svt1)),
  identical(crossprod(svt1[1:6, ], svt2), t(svt1[1:6, ]) \%*\% svt2),
  identical(tcrossprod(svt1, t(svt2)), m12)
)
}
\keyword{array}
\keyword{methods}
