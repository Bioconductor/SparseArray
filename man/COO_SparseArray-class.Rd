\name{COO_SparseArray-class}
\docType{class}

\alias{class:COO_SparseArray}
\alias{COO_SparseArray-class}
\alias{COO_SparseArray}

\alias{class:COO_SparseMatrix}
\alias{COO_SparseMatrix-class}
\alias{COO_SparseMatrix}
\alias{coerce,COO_SparseArray,COO_SparseMatrix-method}
\alias{coerce,COO_SparseMatrix,COO_SparseArray-method}
\alias{coerce,COO_SparseMatrix,SparseArray-method}

\alias{nzcoo}
\alias{nzcoo,COO_SparseArray-method}
\alias{nzvals}
\alias{nzvals,COO_SparseArray-method}

\alias{type,COO_SparseArray-method}
\alias{type<-,COO_SparseArray-method}

\alias{nzcount,COO_SparseArray-method}
\alias{nzwhich,COO_SparseArray-method}

\alias{as.array.COO_SparseArray}
\alias{as.array,COO_SparseArray-method}

\alias{coerce,ANY,COO_SparseArray-method}
\alias{coerce,ANY,COO_SparseMatrix-method}

\alias{coerce,COO_SparseMatrix,dgCMatrix-method}
\alias{coerce,COO_SparseMatrix,lgCMatrix-method}
\alias{coerce,COO_SparseMatrix,dgRMatrix-method}
\alias{coerce,COO_SparseMatrix,lgRMatrix-method}

\alias{coerce,dgCMatrix,COO_SparseMatrix-method}
\alias{coerce,lgCMatrix,COO_SparseMatrix-method}
\alias{coerce,dgRMatrix,COO_SparseMatrix-method}
\alias{coerce,lgRMatrix,COO_SparseMatrix-method}

\alias{coerce,Matrix,COO_SparseArray-method}

\title{COO_SparseArray objects}

\description{
  The COO_SparseArray class is a container for efficient in-memory
  representation of multidimensional sparse arrays. It uses the
  \emph{COO layout} to represent the nonzero data internally.

  A COO_SparseMatrix object is a COO_SparseArray object of 2 dimensions.

  IMPORTANT NOTE: COO_SparseArray and COO_SparseMatrix objects are now
  superseded by the new and more efficient \link{SVT_SparseArray} and
  SVT_SparseMatrix objects.
}

\usage{
## Constructor function:
COO_SparseArray(dim, nzcoo=NULL, nzvals=NULL, dimnames=NULL, check=TRUE)

## Getters (in addition to dim(), length(), and dimnames()):
nzcoo(x)
nzvals(x)
}

\arguments{
  \item{dim}{
    The dimensions (specified as an integer vector) of the
    COO_SparseArray or COO_SparseMatrix object to create.
  }
  \item{nzcoo}{
    A matrix containing the array coordinates of the nonzero values.

    This must be an integer matrix of array coordinates like
    one returned by \code{base::\link[base]{arrayInd}} or
    \code{S4Arrays::\link[S4Arrays]{Lindex2Mindex}}, that is, a matrix
    with \code{length(dim)} columns and where each row is an n-tuple
    representing the coordinates of an array element.
  }
  \item{nzvals}{
    A vector (atomic or list) of length \code{nrow(nzcoo)} containing
    the nonzero values.
  }
  \item{dimnames}{
    The \emph{dimnames} of the object to be created. Must be \code{NULL} or
    a list of length the number of dimensions. Each list element must be
    either \code{NULL} or a character vector along the corresponding dimension.
  }
  \item{check}{
    Should the object be validated upon construction?
  }
  \item{x}{
    A COO_SparseArray or COO_SparseMatrix object.
  }
}

\value{
  \itemize{
    \item For \code{COO_SparseArray()}: A COO_SparseArray or COO_SparseMatrix
          object.

    \item For \code{nzcoo()}: The matrix containing the array coordinates
          of the nonzero values.

    \item For \code{nzvals()}: The vector of nonzero values.
  }
}

\seealso{
  \itemize{
    \item The new \link{SVT_SparseArray} class for a replacement of
          of the COO_SparseArray class.

    \item The \link{SparseArray} class for the virtual parent class of
          COO_SparseArray and SVT_SparseArray.

    \item S4 classes \linkS4class{dgCMatrix} and \linkS4class{lgCMatrix}
          defined in the \pkg{Matrix} package, for the de facto standard
          of sparse matrix representations in the R ecosystem.

    \item \code{base::\link[base]{arrayInd}} in the \pkg{base} package.

    \item \code{S4Arrays::\link[S4Arrays]{Lindex2Mindex}} in the
          \pkg{S4Arrays} package for an improved (faster) version
          of \code{base::arrayInd}.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## EXAMPLE 1
## ---------------------------------------------------------------------
dim1 <- 5:3
nzcoo1 <- Lindex2Mindex(sample(60, 8), 5:3)
nzvals1 <- 11.11 * seq_len(nrow(nzcoo1))
coo1 <- COO_SparseArray(dim1, nzcoo1, nzvals1)
coo1

nzcoo(coo1)
nzvals(coo1)
type(coo1)
sparsity(coo1)

as.array(coo1)  # back to a dense representation

#as.matrix(coo1)  # error!

## ---------------------------------------------------------------------
## EXAMPLE 2
## ---------------------------------------------------------------------
m2 <- matrix(c(5:-2, rep.int(c(0L, 99L), 11)), ncol=6)
coo2 <- as(m2, "COO_SparseArray")
class(coo2)
dim(coo2)
length(coo2)
nzcoo(coo2)
nzvals(coo2)
type(coo2)
sparsity(coo2)

stopifnot(identical(as.matrix(coo2), m2))

t(coo2)
stopifnot(identical(as.matrix(t(coo2)), t(as.matrix(coo2))))

## ---------------------------------------------------------------------
## COERCION FROM/TO dg[C|R]Matrix OR lg[C|R]Matrix OBJECTS
## ---------------------------------------------------------------------
## dg[C|R]Matrix and lg[C|R]Matrix objects are defined in the Matrix
## package.

## dgCMatrix/dgRMatrix:

M2C <- as(coo2, "dgCMatrix")
stopifnot(identical(M2C, as(m2, "dgCMatrix")))

coo2C <- as(M2C, "COO_SparseArray")
## 'coo2C' is the same as 'coo2' except that 'nzvals(coo2C)' has
## type "double" instead of "integer":
stopifnot(all.equal(coo2, coo2C))
typeof(nzvals(coo2C))  # double
typeof(nzvals(coo2))   # integer

M2R <- as(coo2, "dgRMatrix")
stopifnot(identical(M2R, as(m2, "dgRMatrix")))
coo2R <- as(M2R, "COO_SparseArray")
stopifnot(all.equal(as.matrix(coo2), as.matrix(coo2R)))

## lgCMatrix/lgRMatrix:

m3 <- m2 == 99  # logical matrix
coo3 <- as(m3, "COO_SparseArray")
class(coo3)
type(coo3)

M3C <- as(coo3, "lgCMatrix")
stopifnot(identical(M3C, as(m3, "lgCMatrix")))
coo3C <- as(M3C, "COO_SparseArray")
identical(as.matrix(coo3), as.matrix(coo3C))

M3R <- as(coo3, "lgRMatrix")
#stopifnot(identical(M3R, as(m3, "lgRMatrix")))
coo3R <- as(M3R, "COO_SparseArray")
identical(as.matrix(coo3), as.matrix(coo3R))

## ---------------------------------------------------------------------
## A BIG COO_SparseArray OBJECT
## ---------------------------------------------------------------------
nzcoo4 <- cbind(sample(25000, 600000, replace=TRUE),
                sample(195000, 600000, replace=TRUE))
nzvals4 <- runif(600000)
coo4 <- COO_SparseArray(c(25000, 195000), nzcoo4, nzvals4)
coo4
sparsity(coo4)
}
\keyword{classes}
\keyword{methods}
