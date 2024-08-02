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

\alias{nnacount}
\alias{nnacount,NaArray-method}
\alias{nnawhich}
\alias{nnawhich,NaArray-method}

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
    COMING SOON...
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
    COMING SOON...
  }
}

\details{
  COMING SOON...
}

\value{
  COMING SOON...
}

\seealso{
  \itemize{
    \item The \link{SVT_SparseArray} class.

    \item \link{NaArray_aperm} for permuting the dimensions of an
          NaArray object (e.g. transposition).

    \item \link{NaArray_subsetting} for subsetting an NaArray object.

    \item \link{NaArray_abind} for combining 2D or multidimensional
          NaArray objects.

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

# COMING SOON...
}

\keyword{methods}
\keyword{classes}
