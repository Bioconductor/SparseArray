\name{is_nonna}

\alias{is_nonna}
\alias{is_nonna,ANY-method}

\alias{nnacount}
\alias{nnacount,ANY-method}

\alias{nnawhich}
\alias{nnawhich,ANY-method}

\alias{nnavals}
\alias{nnavals,ANY-method}

\alias{nnavals<-}
\alias{nnavals<-,ANY-method}

\title{is_nonna() and the nna*() functions}

\description{
  A set of functions for direct manipulation of the non-NA elements
  of an array-like object.
}

\usage{
is_nonna(x)

nnacount(x)
nnawhich(x, arr.ind=FALSE)
nnavals(x)
nnavals(x) <- value
}

\arguments{
  \item{x}{
    Typically (but not necessarily) an array-like object that is
    \emph{NA-sparse}, like an \link{NaArray} object.

    However, \code{x} can also be an ordinary matrix or array, or any
    matrix-like or array-like object.
  }
  \item{arr.ind}{
    If \code{arr.ind=FALSE} (the default), the indices of the non-NA
    array elements are returned in a numeric vector (a.k.a. \emph{L-index}).
    Otherwise, they're returned in an ordinary matrix (a.k.a. \emph{M-index}).

    See \code{?\link[S4Arrays]{Lindex}} in the \pkg{S4Arrays} package for
    more information about \emph{L-index} and \emph{M-index}, and how to
    convert from one to the other.

    Note that using \code{arr.ind=TRUE} won't work if \code{nnacount(x)}
    is >= \code{.Machine$integer.max} (= 2^31), because, in that case,
    the returned \emph{M-index} would need to be a matrix with more rows
    than what is supported by base R.
  }
  \item{value}{
    A vector, typically of length \code{nnacount(x)} (or 1) and type
    \code{type(x)}.
  }
}

\details{
  \code{nnacount(x)} and \code{nnawhich(x)} are equivalent to, but
  typically more efficient than, \code{sum(is_nonna(x))} and
  \code{which(is_nonna(x))}, respectively.

  \code{nnavals(x)} is equivalent to, but typically more efficient than,
  \code{x[nnawhich(x)]} (or \code{x[is_nonna(x)]}).

  \code{nnavals(x) <- value} replaces the values of the non-NA array
  elements in \code{x} with the supplied values. It's equivalent to,
  but typically more efficient than, \code{x[nnawhich(x)] <- value}.

  Note that \code{nnavals(x) <- nnavals(x)} is guaranteed to be a no-op.
}

\value{
  \code{is_nonna()}: An array-like object of \code{type()} \code{"logical"}
  and same dimensions as the input object.

  \code{nnacount()}: The number of non-NA array elements in \code{x}.

  \code{nnawhich()}: The indices of the non-NA array elements in \code{x},
  either as an \emph{L-index} (if \code{arr.ind} is \code{FALSE}) or as
  an \emph{M-index} (if \code{arr.ind} is \code{TRUE}).
  Note that the indices are returned sorted in strictly ascending order.

  \code{nnavals()}: A vector of the same \code{type()} as \code{x} and
  containing the values of the non-NA array elements in \code{x}.
  Note that the returned vector is guaranteed to be \emph{parallel}
  to \code{nnawhich(x)}.
}

\seealso{
  \itemize{
    \item \link{is_nonzero} for \code{is_nonzero()} and \code{nz*()} functions
          \code{nzcount()}, \code{nzwhich()}, etc...

    \item \link{NaArray} objects.

    \item Ordinary \link[base]{array} objects in base R.

    \item \code{base::\link[base]{which}} in base R.
  }
}

\examples{
a <- array(NA_integer_, dim=c(5, 12, 2))
a[sample(length(a), 20)] <- (-9):10

is_nonna(a)

## Get the number of non-NA array elements in 'a':
nnacount(a)

## nnawhich() returns the indices of the non-NA array elements in 'a'.
## Either as a "L-index" i.e. an integer (or numeric) vector of
## length 'nnacount(a)' containing "linear indices":
nnaidx <- nnawhich(a)
length(nnaidx)
head(nnaidx)

## Or as an "M-index" i.e. an integer matrix with 'nnacount(a)' rows
## and one column per dimension where the rows represent "array indices"
## (a.k.a. "array coordinates"):
Mnnaidx <- nnawhich(a, arr.ind=TRUE)
dim(Mnnaidx)

## Each row in the matrix is an n-tuple representing the "array
## coordinates" of a non-NA element in 'a':
head(Mnnaidx)
tail(Mnnaidx)

## Extract the values of the non-NA array elements in 'a' and return
## them in a vector "parallel" to 'nnawhich(a)':
a_nnavals <- nnavals(a)  # equivalent to 'a[nnawhich(a)]'
length(a_nnavals)
head(a_nnavals)

nnavals(a) <- 10 ^ nnavals(a)
a

## Sanity checks:
stopifnot(
  identical(nnaidx, which(!is.na(a))),
  identical(Mnnaidx, which(!is.na(a), arr.ind=TRUE, useNames=FALSE)),
  identical(nnavals(a), a[nnaidx]),
  identical(nnavals(a), a[Mnnaidx]),
  identical(`nnavals<-`(a, nnavals(a)), a)
)
}

\keyword{array}
\keyword{methods}
