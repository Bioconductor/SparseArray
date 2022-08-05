\name{randomSparseArray}

\alias{randomSparseArray}
\alias{randomSparseMatrix}
\alias{poissonSparseArray}
\alias{poissonSparseMatrix}

\title{Random SparseArray object}

\description{
  \code{randomSparseArray()} and \code{poissonSparseArray()} can be used
  to generate a random \link{SparseArray} object efficiently.
}

\usage{
randomSparseArray(dim, density=0.05)
poissonSparseArray(dim, lambda=-log(0.95), density=NA)

## Convenience wrappers for the 2D case:
randomSparseMatrix(nrow, ncol, density=0.05)
poissonSparseMatrix(nrow, ncol, lambda=-log(0.95), density=NA)
}

\arguments{
  \item{dim}{
    The dimensions (specified as an integer vector) of the \link{SparseArray}
    object to generate.
  }
  \item{density}{
    The desired density (specified as a number >= 0 and <= 1) of the
    \link{SparseArray} object to generate, that is, the ratio between its
    number of nonzero elements and its total number of elements.
    This is \code{nzcount(x)/length(x)} or \code{1 - sparsity(x)}.

    Note that for \code{poissonSparseArray()} and \code{poissonSparseMatrix()}
    \code{density} must be < 1 and the \emph{actual} density of the returned
    object won't be exactly as requested but will typically be very close.
  }
  \item{lambda}{
    The mean of the Poisson distribution. Passed internally to the calls
    to \code{\link[stats]{rpois}()}.

    Only one of \code{lambda} and \code{density} can be specified.

    When \code{density} is requested, \code{rpois()} is called internally
    with \code{lambda} set to \code{-log(1 - density)}. This is expected
    to generate Poisson data with the requested density.

    Finally note that the default value for \code{lambda} corresponds to
    a requested density of 0.05.
  }
  \item{nrow, ncol}{
    Number of rows and columns of the \link{SparseMatrix} object to generate.
  }
}

\details{
  \code{randomSparseArray()} mimics the \code{\link[Matrix]{rsparsematrix}()}
  function from the \pkg{Matrix} package but returns a \link{SparseArray}
  object instead of a dgCMatrix object.

  \code{poissonSparseArray()} populates a \link{SparseArray} object with
  Poisson data i.e. it's equivalent to:
  \preformatted{    a <- array(rpois(prod(dim), lambda), dim)
    as(a, "SparseArray")}
  but is faster and more memory efficient because intermediate dense array
  \code{a} is never generated.
}

\value{
  A \link{SparseArray} derivative (of class \link{SVT_SparseArray} or
  \link{SVT_SparseMatrix}) with the requested dimensions and density.

  The type of the returned object is \code{"double"} for
  \code{randomSparseArray()} and \code{randomSparseMatrix()},
  and \code{"integer"} for \code{poissonSparseArray()} and
  \code{poissonSparseMatrix()}.
}

\note{
  Unlike with \code{Matrix::\link[Matrix]{rsparsematrix}()} there's no
  limit on the number of nonzero elements that can be contained in the
  returned \link{SparseArray} object.

  For example \code{Matrix::rsparsematrix(3e5, 2e4, density=0.5)} will fail
  with an error but \code{randomSparseMatrix(3e5, 2e4, density=0.5)} should
  work (even though it will take some time and the memory footprint of the
  resulting object will be about 18 Gb).
}

\seealso{
  \itemize{
    \item The \code{Matrix::\link[Matrix]{rsparsematrix}} function in
          the \pkg{Matrix} package.

    \item The \code{stats::\link[stats]{rpois}} function in the
          \pkg{stats} package.

    \item \link{SVT_SparseArray} objects.
  }
}

\examples{
## ---------------------------------------------------------------------
## randomSparseArray() / randomSparseMatrix()
## ---------------------------------------------------------------------
set.seed(123)
dgcm1 <- rsparsematrix(2500, 950, density=0.1)
set.seed(123)
svt1 <- randomSparseMatrix(2500, 950, density=0.1)
svt1
type(svt1)  # "double"

stopifnot(identical(as(svt1, "dgCMatrix"), dgcm1))

## ---------------------------------------------------------------------
## poissonSparseArray() / poissonSparseMatrix()
## ---------------------------------------------------------------------
svt2 <- poissonSparseMatrix(2500, 950, density=0.1)
svt2
type(svt2)  # "integer"
1 - sparsity(svt2)  # very close to the requested density

set.seed(123)
svt3 <- poissonSparseArray(c(600, 1700, 80), lambda=0.01)
set.seed(123)
a3 <- array(rpois(length(svt3), lambda=0.01), dim(svt3))
stopifnot(identical(svt3, SparseArray(a3)))

## The memory footprint of 'svt3' is 10x smaller than that of 'a3':
object.size(svt3)
object.size(a3)
as.double(object.size(a3) / object.size(svt3))
}
\keyword{utilities}
