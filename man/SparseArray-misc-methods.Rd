\name{SparseArray-misc-methods}

\alias{SparseArray-misc-methods}
\alias{SparseArray_misc_methods}
\alias{SparseArray-misc}
\alias{SparseArray_misc}

\alias{is.na,COO_SparseArray-method}
\alias{is.infinite,COO_SparseArray-method}
\alias{is.nan,COO_SparseArray-method}
\alias{tolower,COO_SparseArray-method}
\alias{toupper,COO_SparseArray-method}
\alias{nchar,COO_SparseArray-method}

\alias{which,COO_SparseArray-method}

\title{Miscellaneous operations on a SparseArray object}

\description{
  This man page documents various base array operations that are
  supported by \link{SparseArray} derivatives, and that didn't belong
  to any of the groups of operations documented in the other man pages
  of the \pkg{SparseArray} package.

  Note that only \link{COO_SparseArray} objects support these operations
  at the moment.
}

\usage{
\S4method{is.na}{COO_SparseArray}(x)

\S4method{is.infinite}{COO_SparseArray}(x)

\S4method{is.nan}{COO_SparseArray}(x)

\S4method{tolower}{COO_SparseArray}(x)

\S4method{toupper}{COO_SparseArray}(x)

\S4method{nchar}{COO_SparseArray}(x, type="chars", allowNA=FALSE, keepNA=NA)

\S4method{which}{COO_SparseArray}(x, arr.ind=FALSE, useNames=TRUE)
}

\arguments{
  \item{x}{
    An \link{COO_SparseArray} object.
  }
  \item{type, allowNA, keepNA}{
    See \code{?base::\link[base]{nchar}} for a description of
    these arguments.
  }
  \item{arr.ind, useNames}{
    See \code{?base::\link[base]{which}} for a description of
    these arguments.
  }
}

\details{
  More operations will be added in the future. For example
  \link{SVT_SparseArray} objects also need to support \code{which()},
  \code{is.na()}, \code{is.infinite()}, \code{is.nan()}.
}

\value{
  See man pages for the corresponding default methods in the
  \pkg{base} package (e.g. \code{?base::\link[base]{which}},
  \code{?base::\link[base]{is.na}}, etc...) for the value returned
  by these methods.
}

\seealso{
  \itemize{
    \item \code{base::\link[base]{is.na}} and
          \code{base::\link[base]{is.infinite}} in base R.

    \item \code{base::\link[base]{tolower}} in base R.

    \item \code{base::\link[base]{nchar}} in base R.

    \item \code{base::\link[base]{which}} in base R.

    \item \link{SparseArray} objects.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
a <- array(FALSE, dim=5:3)
nzidx <- c(1:2, 8, 10, 15:17, 20, 24, 40, 56:60)
a[nzidx] <- TRUE
coo <- as(a, "COO_SparseArray")

which(coo)

stopifnot(identical(which(coo), as.integer(nzidx)))
}
\keyword{array}
\keyword{methods}
