\name{SparseArray-dim-tuning}

\alias{SparseArray-dim-tuning}
\alias{SparseArray_dim_tuning}
\alias{dim-tuning}
\alias{dim_tuning}
\alias{SparseArray-ineffective-dims}
\alias{SparseArray_ineffective_dims}
\alias{ineffective-dims}
\alias{ineffective_dims}

\alias{drop}
\alias{dim<-}

\title{Add/drop ineffective dims to/from a SparseArray object}

\description{
  The \emph{ineffective dimensions} of an array-like object are its
  dimensions that have an extent of 1.

  Drop all \emph{ineffective dimensions} from \link{SparseArray} object
  \code{x} with \code{drop(x)}.

  Add and/or drop arbitrary \emph{ineffective dimensions} to/from
  \link{SparseArray} object \code{x} with the \code{dim()} setter.
}

\details{
  The \emph{ineffective dimensions} of an array-like object are its
  dimensions that have an extent of 1.
  For example, for a 1 x 1 x 15 x 1 x 6 array, the \emph{ineffective
  dimensions} are its 1st, 2nd, and 4th dimensions.

  Note that \emph{ineffective dimensions} can be dropped or added
  from/to an array-like object \code{x} without changing its length
  (\code{prod(dim(x))}) or altering its content.

  \describe{
    \item{\code{drop(x)}:}{
      Drop all \emph{ineffective dimensions} from \link{SparseArray}
      derivative \code{x}. If \code{x} has at most one effective dimension,
      then the result is returned as an ordinary vector. Otherwise it's
      returned as a \link{SparseArray} derivative.
    }
    \item{\code{dim(x) <- value}:}{
      Add and/or drop arbitrary \emph{ineffective dimensions} to/from
      \link{SparseArray} derivative \code{x}. \code{value} must be a
      vector of dimensions compatible with \code{dim(x)}, that is, it
      must preserve all the effective dimenions in their original order.
    }
  }
}

\seealso{
  \itemize{
    \item \code{\link[base]{drop}()} in base R.

    \item The \code{\link[base]{dim}()} getter and setter in base R.

    \item \link{SparseArray} objects.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
## An array with ineffective 1st and 4th dimensions:
a <- array(0L, dim=c(1, 1, 5, 4, 1, 3))
dimnames(a) <- list(NULL, NULL, letters[1:5], NULL, NULL, LETTERS[1:3])
a[c(1:2, 8, 10, 15:17, 20, 24, 40, 56:60)] <- (1:15)*10L
svt <- SparseArray(a)
dim(svt)

## Drop ineffective dims:
dim(drop(svt))  # the 1st, 2nd, and 5th dimensions were dropped

## Drop some ineffective dims and adds new ones:
svt2 <- svt
dim(svt2) <- c(1, 5, 4, 1, 1, 3, 1)
dim(svt2)

## Sanity check:
stopifnot(identical(as.array(drop(svt)), drop(a)))
a2 <- `dim<-`(a, c(1, 5, 4, 1, 1, 3, 1))
dimnames(a2)[c(2, 6)] <- dimnames(a)[c(3, 6)]
stopifnot(identical(as.array(svt2), a2))
}
\keyword{array}
\keyword{methods}
