\name{NaArray-misc-methods}

\alias{NaArray-misc-methods}
\alias{NaArray_misc_methods}
\alias{NaArray-misc}
\alias{NaArray_misc}

\alias{is.na,NaArray-method}
\alias{is.nan,NaArray-method}
\alias{is.infinite,NaArray-method}

\title{Miscellaneous operations on a NaArray object}

\description{
  This man page documents various base array operations that are
  supported by \link{NaArray} derivatives, and that didn't belong
  to any of the groups of operations documented in the other man pages
  of the \pkg{SparseArray} package.
}

\usage{
# --- unary isometric array transformations ---

\S4method{is.nan}{NaArray}(x)

\S4method{is.infinite}{NaArray}(x)

# --- N-ary isometric array transformations ---

# COMING SOON...
}

\arguments{
  \item{x}{
    An \link{NaArray} object.
  }
}

\details{
  More operations will be added in the future.
}

\value{
  \code{is.nan()} and \code{is.infinite()} return a \link{SparseArray}
  object of \code{type()} \code{"logical"} and same dimensions as the
  input object.
}

\seealso{
  \itemize{
    \item \code{base::\link[base]{is.nan}} and
          \code{base::\link[base]{is.infinite}} in base R.

    \item \link{NaArray} objects.

    \item \link{SparseArray} objects.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
a <- array(c(NA, 2.77, NaN, Inf, NA, -Inf), dim=5:3)
naa <- NaArray(a)  # NaArray object
naa

is.nan(naa)            # SparseArray object of type "logical"
is.infinite(naa)       # SparseArray object of type "logical"

## Sanity checks:
res <- is.nan(naa)
stopifnot(is(res, "SparseArray"), type(res) == "logical",
          identical(as.array(res), is.nan(a)))
res <- is.infinite(naa)
stopifnot(is(res, "SparseArray"), type(res) == "logical",
          identical(as.array(res), is.infinite(a)))
}
\keyword{array}
\keyword{methods}
