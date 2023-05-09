\name{SparseArray-ineffective-dims}

\alias{SparseArray-ineffective-dims}
\alias{SparseArray_ineffective_dims}

\alias{drop}
\alias{drop,SVT_SparseArray-method}

\title{Drop/add ineffective dims from/to a SparseArray object}

\description{
  Drop/add \emph{ineffective} dimensions from/to a \link{SparseArray}
  derivative.

  NOTE: Only dropping ineffective dims is supported at the moment.
  Adding ineffective dims will be supported soon.
}

\usage{
\S4method{drop}{SVT_SparseArray}(x)
}

\arguments{
  \item{x}{
    An \link{SVT_SparseArray} object.
  }
}

\details{
  The \emph{ineffective} dimensions of an array-like object are its
  dimensions that have an extent of 1.
  For example, for a 1 x 1 x 15 x 1 x 6 array, the \emph{ineffective}
  dimensions are its 1st, 2nd, and 4th dimensions.

  Note that \emph{ineffective} dimensions can be dropped or added
  from/to an array-like object \code{x} without changing its length
  (\code{prod(dim(x))}) or altering its content.
}

\value{
  An \link{SVT_SparseArray} object with the same content as \code{x},
  but with its \emph{ineffective} dimensions dropped.
}

\seealso{
  \itemize{
    \item \code{\link[base]{drop}} in base R.

    \item \link{SparseArray} objects.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
## An array with ineffective 1st and 4th dimensions:
a <- array(0L, dim=c(1, 5, 4, 1, 3))
dimnames(a) <- list(NULL, letters[1:5], NULL, NULL, LETTERS[1:3])
a[c(1:2, 8, 10, 15:17, 20, 24, 40, 56:60)] <- (1:15)*10L
svt <- SparseArray(a)

dim(svt)

dim(drop(svt))  # the 1st and 4th dimensions were dropped

## Sanity check:
stopifnot(identical(as.array(drop(svt)), drop(a)))
}
\keyword{array}
\keyword{methods}
