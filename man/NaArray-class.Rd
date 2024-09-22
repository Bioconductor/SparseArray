\name{NaArray}
\docType{class}

\alias{class:NaArray}
\alias{NaArray-class}
\alias{NaArray}

\alias{class:NaMatrix}
\alias{NaMatrix-class}
\alias{NaMatrix}

\alias{coerce,NaArray,NaMatrix-method}
\alias{coerce,NaMatrix,NaArray-method}

\alias{dim,NaArray-method}
\alias{dimnames,NaArray-method}
\alias{dimnames<-,NaArray,ANY-method}
\alias{type,NaArray-method}
\alias{type<-,NaArray-method}

\alias{nnacount,NaArray-method}
\alias{nnawhich,NaArray-method}

\alias{as.array.NaArray}
\alias{as.array,NaArray-method}
\alias{coerce,array,NaArray-method}
\alias{coerce,matrix,NaMatrix-method}

\alias{show,NaArray-method}

\title{NaArray objects}

\description{
  EXPERIMENTAL!!!

  Like \link{SVT_SparseArray} objects but the background value is \code{NA}
  instead of zero.
}

\usage{
## Constructor function:
NaArray(x, dim=NULL, dimnames=NULL, type=NA)
}

\arguments{
  \item{x}{
    If \code{dim} is \code{NULL} (the default) then \code{x} must be an
    ordinary matrix or array, or a dgCMatrix/lgCMatrix object, or any
    matrix-like or array-like object that supports coercion to NaArray.

    If \code{dim} is provided then \code{x} can either be missing, a
    vector (atomic or list), or an array-like object. If missing, then
    an all-NA NaArray object will be constructed. Otherwise \code{x}
    will be used to fill the returned object (\code{length(x)} must be
    \code{<= prod(dim)}). Note that if \code{x} is an array-like object
    then its dimensions are ignored i.e. it's treated as a vector.
  }
  \item{dim}{
    \code{NULL} or the dimensions (supplied as an integer vector) of
    the NaArray or NaMatrix object to construct.
    If \code{NULL} (the default) then the returned object will have
    the dimensions of matrix-like or array-like object \code{x}.
  }
  \item{dimnames}{
    The \emph{dimnames} of the object to construct. Must be \code{NULL} or
    a list of length the number of dimensions. Each list element must be
    either \code{NULL} or a character vector along the corresponding dimension.
    If both \code{dim} and \code{dimnames} are \code{NULL} (the default) then
    the returned object will have the \emph{dimnames} of matrix-like or
    array-like object \code{x}.
  }
  \item{type}{
    A single string specifying the requested type of the object.

    By default the NaArray object returned by the constructor
    function will have the same \code{type()} as \code{x}. However the
    user can use the \code{type} argument to request a different type.
    Note that doing:
    \preformatted{    naa <- NaArray(x, type=type)}
    is equivalent to doing:
    \preformatted{    naa <- NaArray(x)
    type(naa) <- type}
    but the former is more convenient and will generally be more efficient.

    The supported types for NaArray objects are \code{"integer"},
    \code{"logical"}, \code{"double"}, \code{"complex"}, and
    \code{"character"}.
  }
}

\details{
  NaArray is a concrete subclass of the \link{Array} virtual class.
  This makes NaArray objects Array derivatives.

  Like with \link{SVT_SparseArray} objects, the non-NA data in an NaArray
  object is stored in a \emph{Sparse Vector Tree}.
  See \code{?\link{SVT_SparseArray}} for more information.
}

\value{
  An NaArray or NaMatrix object.
}

\seealso{
  \itemize{
    \item The \link{SVT_SparseArray} class.

    \item \link{nnafunctions} for \code{nna*()} functions \code{nnacount()},
          \code{nnawhich()}, etc...

    \item \link{NaArray_aperm} for permuting the dimensions of an
          NaArray object (e.g. transposition).

    \item \link{NaArray_subsetting} for subsetting an NaArray object.

    \item \link{NaArray_abind} for combining 2D or multidimensional
          NaArray objects.

    \item \link{NaArray_summarization} for NaArray summarization
          methods.

    \item \link{NaArray_Arith}, \link{NaArray_Compare}, and
          \link{NaArray_Logic}, for operations from the \code{Arith},
          \code{Compare}, and \code{Arith} groups on NaArray objects.

    \item \link{NaArray_Math} for operations from the \code{Math} and
          \code{Math2} groups on NaArray objects.

    \item \link{NaArray_matrixStats} for col/row summarization methods
          for NaArray objects.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## Display details of class definition & known subclasses
## ---------------------------------------------------------------------

showClass("NaArray")

## ---------------------------------------------------------------------
## The NaArray() constructor
## ---------------------------------------------------------------------

naa1 <- NaArray(dim=5:3)  # all-NA object
naa1

naa2 <- NaArray(dim=c(35000, 2e6), type="integer")  # all-NA object
naa2

## Add some non-NA values to 'naa2':
naa2[cbind( 1:99, 2:100)] <-  1L
naa2[cbind(1:100, 1:100)] <-  0L
naa2[cbind(2:100,  1:99)] <- -1L
naa2

## The dimnames can be specified at construction time, or
## added/modified later:
naa3 <- NaArray(c(NA, NA, 1L, NA, 0:7, rep(NA, 4), 12:14, NA),
                dim=4:5, dimnames=list(letters[1:4], LETTERS[1:5]))
naa3

colnames(naa3) <- LETTERS[22:26]
naa3

## Sanity checks:
stopifnot(
  is(naa1, "NaArray"),
  identical(dim(naa1), 5:3),
  identical(as.array(naa1), array(dim=5:3)),
  is(naa2, "NaMatrix"),
  all.equal(dim(naa2), c(35000, 2e6)),
  identical(nnacount(naa2), 298L),
  is(naa3, "NaMatrix"),
  identical(dim(naa3), 4:5),
  identical(nnacount(naa3), 12L)
)
}

\keyword{methods}
\keyword{classes}
