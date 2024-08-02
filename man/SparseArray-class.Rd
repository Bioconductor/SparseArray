\name{SparseArray}
\docType{class}

\alias{class:SparseArray}
\alias{SparseArray-class}
\alias{SparseArray}

\alias{class:SparseMatrix}
\alias{SparseMatrix-class}
\alias{SparseMatrix}

\alias{dim,SparseArray-method}
\alias{dimnames,SparseArray-method}
\alias{dimnames<-,SparseArray,ANY-method}

\alias{is_sparse,SparseArray-method}

\alias{nzcount}
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

\alias{show,SparseArray-method}

\title{SparseArray objects}

\description{
  The \pkg{SparseArray} package defines the SparseArray virtual class
  whose purpose is to be extended by other S4 classes that aim
  at representing in-memory multidimensional sparse arrays.

  It has currently two concrete subclasses, \link{COO_SparseArray}
  and \link{SVT_SparseArray}, both also defined in this package.
  Each subclass uses its own internal representation for the nonzero
  multidimensional data, the \emph{COO layout} for \link{COO_SparseArray},
  and the \emph{SVT layout} for \link{SVT_SparseArray}. The two layouts
  are described in the \link{COO_SparseArray} and \link{SVT_SparseArray}
  man pages, respectively.

  Finally, the package also defines the SparseMatrix virtual class, as
  a subclass of the SparseArray class, for the specific 2D case.
}

\usage{
## Constructor function:
SparseArray(x, type=NA)
}

\arguments{
  \item{x}{
    An ordinary matrix or array, or a dg[C|R]Matrix object, or an
    lg[C|R]Matrix object, or any matrix-like or array-like object that
    supports coercion to \link{SVT_SparseArray}.
  }
  \item{type}{
    A single string specifying the requested type of the object.

    Normally, the SparseArray object returned by the constructor
    function has the same \code{type()} as \code{x} but the user can use
    the \code{type} argument to request a different type. Note that doing:
    \preformatted{    sa <- SparseArray(x, type=type)}
    is equivalent to doing:
    \preformatted{    sa <- SparseArray(x)
    type(sa) <- type}
    but the former is more convenient and will generally be more efficient.

    Supported types are all R atomic types plus \code{"list"}.
  }
}

\details{
  The SparseArray class extends the \link[S4Arrays]{Array} virtual class
  defined in the \pkg{S4Arrays} package. Here is the full SparseArray
  sub-hierarchy as defined in the \pkg{SparseArray} package (virtual
  classes are marked with an asterisk):
  \preformatted{
: Array class :                 Array*
: hierarchy   :                   ^
                                  |
- - - - - - - - - - - - - - - - - | - - - - - - - - - - - - - - -
: SparseArray   :            SparseArray*
: sub-hierarchy :            ^    ^    ^
                             |    |    |
                 COO_SparseArray  |  SVT_SparseArray
                        ^         |         ^
- - - - - - - - - - - - | - - - - | - - - - | - - - - - - - - - -
: SparseMatrix      :   |    SparseMatrix*  |
: sub-sub-hierarchy :   |    ^         ^    |
                        |    |         |    |
                 COO_SparseMatrix    SVT_SparseMatrix}

  Any object that belongs to a class that extends SparseArray e.g. (a
  \link{SVT_SparseArray} or \link{SVT_SparseMatrix} object) is called
  a \emph{SparseArray derivative}.

  Most of the \emph{standard matrix and array API} defined in base R should
  work on SparseArray derivatives, including \code{dim()}, \code{length()},
  \code{dimnames()}, \code{`dimnames<-`()}, \code{[}, \code{drop()},
  \code{`[<-`} (subassignment), \code{t()}, \code{rbind()}, \code{cbind()},
  etc...

  SparseArray derivatives also support \code{type()}, \code{`type<-`()},
  \code{is_sparse()}, \code{nzcount()}, \code{nzwhich()}, \code{nzvals()},
  \code{`nzvals<-`()}, \code{sparsity()}, \code{arbind()}, and \code{acbind()}.

  \code{sparsity(x)} returns the ratio between the number of zero-valued
  elements in array-like object \code{x} and its total number of elements
  (\code{length(x)} or \code{prod(dim(x))}).
  More precisely, \code{sparsity(x)} is \code{1 - nzcount(x)/length(x)}.
}

\value{
  A \emph{SparseArray derivative}, that is a \link{SVT_SparseArray},
  \link{COO_SparseArray}, \link{SVT_SparseMatrix}, or
  \link{COO_SparseMatrix} object.

  The \code{type()} of the input object is preserved, except if a
  different one was requested via the \code{type} argument.

  What is considered a zero depends on the \code{type()}:
  \itemize{
    \item \code{"logical"} zero is \code{FALSE};
    \item \code{"integer"} zero is \code{0L};
    \item \code{"double"} zero is \code{0};
    \item \code{"complex"} zero is \code{0+0i};
    \item \code{"raw"} zero is \code{raw(1)};
    \item \code{"character"} zero is \code{""} (empty string);
    \item \code{"list"} zero is \code{NULL}.
  }
}

\seealso{
  \itemize{
    \item The \link{COO_SparseArray} and \link{SVT_SparseArray} classes.

    \item \link{SparseArray_aperm} for permuting the dimensions of a
          SparseArray object (e.g. transposition).

    \item \link{SparseArray_subsetting} for subsetting a SparseArray object.

    \item \link{SparseArray_subassignment} for SparseArray subassignment.

    \item \link{SparseArray_abind} for combining 2D or multidimensional
          SparseArray objects.

    \item \link{SparseArray_summarization} for SparseArray summarization
          methods.

    \item \link{SparseArray_Ops} for operations from the \code{Ops}
          group on SparseArray objects.

    \item \link{SparseArray_Math} for operations from the \code{Math} and
          \code{Math2} groups on SparseArray objects.

    \item \link{SparseArray_Complex} for operations from the \code{Complex}
          group on SparseArray objects.

    \item \link{SparseArray_misc} for miscellaneous operations on a
          SparseArray object.

    \item \link{SparseMatrix_mult} for SparseMatrix multiplication and
          cross-product.

    \item \link{matrixStats_methods} for SparseArray col/row summarization
          methods.

    \item \link{rowsum_methods} for \code{rowsum()} methods for sparse
          matrices.

    \item \code{\link{randomSparseArray}} to generate a random SparseArray
          object.

    \item \code{\link{readSparseCSV}} to read/write a sparse matrix
          from/to a CSV (comma-separated values) file.

    \item S4 classes \linkS4class{dgCMatrix}, \linkS4class{dgRMatrix},
          and \linkS4class{lgCMatrix} defined in the \pkg{Matrix} package,
          for the de facto standard for sparse matrix representations
          in the R ecosystem.

    \item \code{\link[S4Arrays]{is_sparse}} in the \pkg{S4Arrays} package.

    \item The \link[S4Arrays]{Array} class defined in the \pkg{S4Arrays}
          package.

    \item Ordinary \link[base]{array} objects in base R.

    \item \code{base::\link[base]{which}} in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## Display details of class definition & known subclasses
## ---------------------------------------------------------------------

showClass("SparseArray")

## ---------------------------------------------------------------------
## The SparseArray() constructor
## ---------------------------------------------------------------------

a <- array(rpois(9e6, lambda=0.3), dim=c(500, 3000, 6))
SparseArray(a)    # an SVT_SparseArray object

m <- matrix(rpois(9e6, lambda=0.3), ncol=500)
SparseArray(m)    # an SVT_SparseMatrix object

dgc <- sparseMatrix(i=c(4:1, 2:4, 9:12, 11:9), j=c(1:7, 1:7),
                    x=runif(14), dims=c(12, 7))
class(dgc)
SparseArray(dgc)  # an SVT_SparseMatrix object

dgr <- as(dgc, "RsparseMatrix")
class(dgr)
SparseArray(dgr)  # a COO_SparseMatrix object

## ---------------------------------------------------------------------
## nzcount(), nzwhich(), nzvals(), `nzvals<-`()
## ---------------------------------------------------------------------
x <- SparseArray(a)

## Get the number of nonzero array elements in 'x':
nzcount(x)

## nzwhich() returns the indices of the nonzero array elements in 'x'.
## Either as an integer (or numeric) vector of length 'nzcount(x)'
## containing "linear indices":
nzidx <- nzwhich(x)
length(nzidx)
head(nzidx)

## Or as an integer matrix with 'nzcount(x)' rows and one column per
## dimension where the rows represent "array indices" (a.k.a. "array
## coordinates"):
Mnzidx <- nzwhich(x, arr.ind=TRUE)
dim(Mnzidx)

## Each row in the matrix is an n-tuple representing the "array
## coordinates" of a nonzero element in 'x':
head(Mnzidx)
tail(Mnzidx)

## Extract the values of the nonzero array elements in 'x' and return
## them in a vector "parallel" to 'nzwhich(x)':
x_nzvals <- nzvals(x)  # equivalent to 'x[nzwhich(x)]'
length(x_nzvals)
head(x_nzvals)

nzvals(x) <- log1p(nzvals(x))
x

## Sanity checks:
stopifnot(identical(nzidx, which(a != 0)))
stopifnot(identical(Mnzidx, which(a != 0, arr.ind=TRUE, useNames=FALSE)))
stopifnot(identical(x_nzvals, a[nzidx]))
stopifnot(identical(x_nzvals, a[Mnzidx]))
stopifnot(identical(`nzvals<-`(x, nzvals(x)), x))
}

\keyword{methods}
\keyword{classes}
