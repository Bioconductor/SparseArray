\name{is_nonzero}

\alias{is_nonzero}
\alias{is_nonzero,ANY-method}
\alias{is_nonzero,sparseMatrix-method}

\alias{nzcount}
\alias{nzcount,ANY-method}
\alias{nzcount,CsparseMatrix-method}
\alias{nzcount,RsparseMatrix-method}

\alias{nzwhich}
\alias{nzwhich,ANY-method}
\alias{nzwhich,CsparseMatrix-method}
\alias{nzwhich,RsparseMatrix-method}

\alias{nzvals}
\alias{nzvals,ANY-method}
\alias{nzvals,dgCMatrix-method}
\alias{nzvals,lgCMatrix-method}
\alias{nzvals,ngCMatrix-method}
\alias{nzvals,ngRMatrix-method}

\alias{nzvals<-}
\alias{nzvals<-,ANY-method}

\alias{sparsity}

\title{is_nonzero() and the nz*() functions}

\description{
  A set of functions for direct manipulation of the nonzero elements
  of an array-like object.
}

\usage{
is_nonzero(x)

nzcount(x)
nzwhich(x, arr.ind=FALSE)
nzvals(x)
nzvals(x) <- value

sparsity(x)
}

\arguments{
  \item{x}{
    Typically (but not necessarily) an array-like object that is sparse,
    like a \link{SparseArray} derivative, or a dg[C|R]Matrix or lg[C|R]Matrix
    object from the \pkg{Matrix} package.

    However, \code{x} can also be an ordinary matrix or array, or any
    matrix-like or array-like object.
  }
  \item{arr.ind}{
    If \code{arr.ind=FALSE} (the default), the indices of the nonzero
    array elements are returned in a numeric vector (a.k.a. \emph{L-index}).
    Otherwise, they're returned in an ordinary matrix (a.k.a. \emph{M-index}).

    See \code{?\link[S4Arrays]{Lindex}} in the \pkg{S4Arrays} package for
    more information about \emph{L-index} and \emph{M-index}, and how to
    convert from one to the other.

    Note that using \code{arr.ind=TRUE} won't work if \code{nzcount(x)}
    is >= \code{.Machine$integer.max} (= 2^31), because, in that case,
    the returned \emph{M-index} would need to be a matrix with more rows
    than what is supported by base R.
  }
  \item{value}{
    A vector, typically of length \code{nzcount(x)} (or 1) and type
    \code{type(x)}.
  }
}

\details{
  \code{nzcount(x)} and \code{nzwhich(x)} are equivalent to, but
  typically more efficient than, \code{sum(is_nonzero(x))} and
  \code{which(is_nonzero(x))}, respectively.

  \code{nzvals(x)} is equivalent to, but typically more efficient than,
  \code{x[nzwhich(x)]} (or \code{x[is_nonzero(x)]}).

  \code{nzvals(x) <- value} replaces the values of the nonzero array
  elements in \code{x} with the supplied values. It's equivalent to,
  but typically more efficient than, \code{x[nzwhich(x)] <- value}.

  Note that \code{nzvals(x) <- nzvals(x)} is guaranteed to be a no-op.
}

\value{
  \code{is_nonzero()}: An array-like object of \code{type()} \code{"logical"}
  and same dimensions as the input object.

  \code{nzcount()}: The number of nonzero array elements in \code{x}.

  \code{nzwhich()}: The indices of the nonzero array elements in \code{x},
  either as an \emph{L-index} (if \code{arr.ind} is \code{FALSE}) or as
  an \emph{M-index} (if \code{arr.ind} is \code{TRUE}).
  Note that the indices are returned sorted in strictly ascending order.

  \code{nzvals()}: A vector of the same \code{type()} as \code{x} and
  containing the values of the nonzero array elements in \code{x}.
  Note that the returned vector is guaranteed to be \emph{parallel}
  to \code{nzwhich(x)}.

  \code{sparsity(x)}: The ratio between the number of zero-valued
  elements in array-like object \code{x} and its total number of
  elements (\code{length(x)} or \code{prod(dim(x))}).
  More precisely, \code{sparsity(x)} is \code{1 - nzcount(x)/length(x)}.
}

\seealso{
  \itemize{
    \item \link{is_nonna} for \code{is_nonna()} and \code{nna*()} functions
          \code{nnacount()}, \code{nnawhich()}, etc...

    \item \link{SparseArray} objects.

    \item S4 classes \linkS4class{dgCMatrix}, \linkS4class{lgCMatrix},
          and \linkS4class{ngCMatrix} defined in the \pkg{Matrix} package.

    \item Ordinary \link[base]{array} objects in base R.

    \item \code{base::\link[base]{which}} in base R.
  }
}

\examples{
a <- array(rpois(120, lambda=0.3), dim=c(5, 12, 2))

is_nonzero(a)

## Get the number of nonzero array elements in 'a':
nzcount(a)

## nzwhich() returns the indices of the nonzero array elements in 'a'.
## Either as a "L-index" i.e. an integer (or numeric) vector of
## length 'nzcount(a)' containing "linear indices":
nzidx <- nzwhich(a)
length(nzidx)
head(nzidx)

## Or as an "M-index" i.e. an integer matrix with 'nzcount(a)' rows
## and one column per dimension where the rows represent "array indices"
## (a.k.a. "array coordinates"):
Mnzidx <- nzwhich(a, arr.ind=TRUE)
dim(Mnzidx)

## Each row in the matrix is an n-tuple representing the "array
## coordinates" of a nonzero element in 'a':
head(Mnzidx)
tail(Mnzidx)

## Extract the values of the nonzero array elements in 'a' and return
## them in a vector "parallel" to 'nzwhich(a)':
a_nzvals <- nzvals(a)  # equivalent to 'a[nzwhich(a)]'
length(a_nzvals)
head(a_nzvals)

nzvals(a) <- log1p(nzvals(a))
a

## Sanity checks:
stopifnot(
  identical(nzidx, which(a != 0)),
  identical(Mnzidx, which(a != 0, arr.ind=TRUE, useNames=FALSE)),
  identical(nzvals(a), a[nzidx]),
  identical(nzvals(a), a[Mnzidx]),
  identical(`nzvals<-`(a, nzvals(a)), a)
)
}

\keyword{array}
\keyword{methods}
