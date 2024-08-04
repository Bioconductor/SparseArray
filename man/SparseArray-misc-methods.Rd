\name{SparseArray-misc-methods}

\alias{SparseArray-misc-methods}
\alias{SparseArray_misc_methods}
\alias{SparseArray-misc}
\alias{SparseArray_misc}

\alias{is.na}
\alias{is.na,COO_SparseArray-method}
\alias{is.na,SVT_SparseArray-method}
\alias{is.nan}
\alias{is.nan,COO_SparseArray-method}
\alias{is.nan,SVT_SparseArray-method}
\alias{is.infinite}
\alias{is.infinite,COO_SparseArray-method}
\alias{is.infinite,SVT_SparseArray-method}
\alias{tolower}
\alias{tolower,COO_SparseArray-method}
\alias{toupper}
\alias{toupper,COO_SparseArray-method}
\alias{nchar}
\alias{nchar,COO_SparseArray-method}

\alias{pmin}
\alias{pmin,SparseArray-method}
\alias{pmax}
\alias{pmax,SparseArray-method}

\title{Miscellaneous operations on a SparseArray object}

\description{
  This man page documents various base array operations that are
  supported by \link{SparseArray} derivatives, and that didn't belong
  to any of the groups of operations documented in the other man pages
  of the \pkg{SparseArray} package.
}

\usage{
# --- unary isometric array transformations ---

\S4method{is.na}{COO_SparseArray}(x)
\S4method{is.na}{SVT_SparseArray}(x)

\S4method{is.nan}{COO_SparseArray}(x)
\S4method{is.nan}{SVT_SparseArray}(x)

\S4method{is.infinite}{COO_SparseArray}(x)
\S4method{is.infinite}{SVT_SparseArray}(x)

\S4method{tolower}{COO_SparseArray}(x)

\S4method{toupper}{COO_SparseArray}(x)

\S4method{nchar}{COO_SparseArray}(x, type="chars", allowNA=FALSE, keepNA=NA)

# --- N-ary isometric array transformations ---

\S4method{pmin}{SparseArray}(..., na.rm=FALSE)
\S4method{pmax}{SparseArray}(..., na.rm=FALSE)
}

\arguments{
  \item{x}{
    A \link{SparseArray} derivative.
  }
  \item{type, allowNA, keepNA}{
    See \code{?base::\link[base]{nchar}} for a description of
    these arguments.
  }
  \item{...}{
    \link{SparseArray} derivatives.
  }
  \item{na.rm}{
    See \code{?base::\link[base]{pmin}} for a description of
    this argument.
  }
}

\details{
  More operations will be supported in the future.
}

\value{
  See man pages of the corresponding base functions (e.g.
  \code{?base::\link[base]{is.na}}, \code{?base::\link[base]{nchar}},
  \code{?base::\link[base]{pmin}}, etc...) for the value returned by
  these methods.

  Note that, like the base functions, the methods documented in this man page
  are \emph{endomorphisms} i.e. they return an array-like object of the same
  class as the input.
}

\seealso{
  \itemize{
    \item \code{base::\link[base]{is.na}} and
          \code{base::\link[base]{is.infinite}} in base R.

    \item \code{base::\link[base]{tolower}} in base R.

    \item \code{base::\link[base]{nchar}} in base R.

    \item \code{base::\link[base]{pmin}} in base R.

    \item \link{SparseArray} objects.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
a <- array(c(0, 2.77, NA, 0, NaN, -Inf), dim=5:3)
svt <- SparseArray(a)  # SVT_SparseArray object
class(svt)

is.na(svt)             # SVT_SparseArray object of type "logical"
is.nan(svt)            # SVT_SparseArray object of type "logical"
is.infinite(svt)       # SVT_SparseArray object of type "logical"

svt1 <- poissonSparseMatrix(500, 20, density=0.2)
svt2 <- poissonSparseMatrix(500, 20, density=0.25) * 0.77
pmin(svt1, svt2)
pmax(svt1, svt2)

## Sanity checks:
res <- is.na(svt)
stopifnot(is(res, "SVT_SparseArray"), type(res) == "logical",
          identical(as.array(res), is.na(a)))
res <- is.nan(svt)
stopifnot(is(res, "SVT_SparseArray"), type(res) == "logical",
          identical(as.array(res), is.nan(a)))
res <- is.infinite(svt)
stopifnot(is(res, "SVT_SparseArray"), type(res) == "logical",
          identical(as.array(res), is.infinite(a)))
res <- pmin(svt1, svt2)
stopifnot(is(res, "SVT_SparseArray"),
          identical(as.array(res), pmin(as.array(svt1), as.array(svt2))))
res <- pmax(svt1, svt2)
stopifnot(is(res, "SVT_SparseArray"),
          identical(as.array(res), pmax(as.array(svt1), as.array(svt2))))
}
\keyword{array}
\keyword{methods}
