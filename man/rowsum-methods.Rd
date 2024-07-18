\name{rowsum-methods}

\alias{rowsum-methods}
\alias{rowsum_methods}

\alias{rowsum}
\alias{rowsum.SparseMatrix}
\alias{rowsum,SparseMatrix-method}
\alias{rowsum.dgCMatrix}
\alias{rowsum,dgCMatrix-method}
\alias{colsum}
\alias{colsum,SparseMatrix-method}
\alias{colsum,dgCMatrix-method}

\title{rowsum() methods for sparse matrices}

\description{
  The \pkg{SparseArray} package provides memory-efficient
  \code{\link[S4Arrays]{rowsum}()} and \code{\link[S4Arrays]{colsum}()}
  methods for \link{SparseMatrix} and \linkS4class{dgCMatrix} objects.
}

\usage{
\S4method{rowsum}{SparseMatrix}(x, group, reorder=TRUE, ...)

\S4method{rowsum}{dgCMatrix}(x, group, reorder=TRUE, ...)

\S4method{colsum}{SparseMatrix}(x, group, reorder=TRUE, ...)

\S4method{colsum}{dgCMatrix}(x, group, reorder=TRUE, ...)
}

\arguments{
  \item{x}{
    A \link{SparseMatrix} or \linkS4class{dgCMatrix} object.
  }
  \item{group, reorder}{
    See \code{?base::\link[base]{rowsum}} for a description of
    these arguments.
  }
  \item{...}{
    Like the default S3 \code{rowsum()} method defined in the \pkg{base}
    package, the methods documented in this man page support additional
    argument \code{na.rm}, set to \code{FALSE} by default.
    If \code{TRUE}, missing values (\code{NA} or \code{NaN}) are omitted
    from the calculations.
  }
}

\value{
  An \emph{ordinary} matrix, like the default \code{rowsum()} method.
  See \code{?base::\link[base]{rowsum}} for how the matrix returned
  by the default \code{rowsum()} method is obtained.
}

\seealso{
  \itemize{
    \item \code{\link[base]{rowsum}} in base R.

    \item \code{S4Arrays::\link[S4Arrays]{rowsum}} in the \pkg{S4Arrays}
          package for the \code{rowsum()} and \code{colsum()} S4 generic
          functions.

    \item \link{SparseMatrix} objects.

    \item \linkS4class{dgCMatrix} objects implemented in the \pkg{Matrix}
          package.
  }
}

\examples{
svt0 <- randomSparseMatrix(7e5, 100, density=0.15)
dgcm0 <- as(svt0, "dgCMatrix")
m0 <- as.matrix(svt0)

group <- sample(10, nrow(m0), replace=TRUE)

## Calling rowsum() on the sparse representations is usually faster
## than on the dense representation:
rs1 <- rowsum(m0, group)
rs2 <- rowsum(svt0, group)   # about 3x faster
rs3 <- rowsum(dgcm0, group)  # also about 3x faster

## Sanity checks:
stopifnot(identical(rs1, rs2))
stopifnot(identical(rs1, rs3))
}
\keyword{array}
\keyword{methods}
\keyword{algebra}
\keyword{arith}
