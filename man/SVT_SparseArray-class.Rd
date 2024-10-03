\name{SVT_SparseArray-class}
\docType{class}

\alias{class:NULL_OR_list}
\alias{NULL_OR_list-class}
\alias{NULL_OR_list}

\alias{class:SVT_SparseArray}
\alias{SVT_SparseArray-class}
\alias{SVT_SparseArray}

\alias{class:SVT_SparseMatrix}
\alias{SVT_SparseMatrix-class}
\alias{SVT_SparseMatrix}
\alias{coerce,SVT_SparseArray,SVT_SparseMatrix-method}
\alias{coerce,SVT_SparseMatrix,SVT_SparseArray-method}
\alias{coerce,SVT_SparseMatrix,SparseArray-method}

\alias{type,SVT_SparseArray-method}
\alias{type<-,SVT_SparseArray-method}

\alias{is_nonzero,SVT_SparseArray-method}
\alias{nzcount,SVT_SparseArray-method}
\alias{nzwhich,SVT_SparseArray-method}

\alias{as.array.SVT_SparseArray}
\alias{as.array,SVT_SparseArray-method}
\alias{coerce,array,SVT_SparseArray-method}
\alias{coerce,array,SparseArray-method}
\alias{coerce,matrix,SVT_SparseMatrix-method}
\alias{coerce,matrix,SparseMatrix-method}

\alias{coerce,SVT_SparseMatrix,dgCMatrix-method}
\alias{coerce,SVT_SparseMatrix,lgCMatrix-method}
\alias{coerce,SVT_SparseMatrix,ngCMatrix-method}
\alias{coerce,CsparseMatrix,SVT_SparseMatrix-method}
\alias{coerce,CsparseMatrix,SparseMatrix-method}

\alias{coerce,Matrix,SVT_SparseArray-method}
\alias{coerce,Matrix,SparseArray-method}

\alias{coerce,SVT_SparseArray,COO_SparseArray-method}
\alias{coerce,SVT_SparseMatrix,COO_SparseMatrix-method}
\alias{coerce,COO_SparseArray,SVT_SparseArray-method}
\alias{coerce,COO_SparseMatrix,SVT_SparseMatrix-method}

\alias{coerce,ANY,SparseArray-method}
\alias{coerce,ANY,SparseMatrix-method}
\alias{coerce,ANY,SVT_SparseArray-method}
\alias{coerce,ANY,SVT_SparseMatrix-method}

\title{SVT_SparseArray objects}

\description{
  The SVT_SparseArray class is a new container for efficient in-memory
  representation of multidimensional sparse arrays. It uses the
  \emph{SVT layout} to represent the nonzero multidimensional data
  internally.

  An SVT_SparseMatrix object is an SVT_SparseArray object of 2 dimensions.

  Note that SVT_SparseArray and SVT_SparseMatrix objects replace the older
  and less efficient \link{COO_SparseArray} and COO_SparseMatrix objects.
}

\usage{
## Constructor function:
SVT_SparseArray(x, dim=NULL, dimnames=NULL, type=NA)
}

\arguments{
  \item{x}{
    If \code{dim} is \code{NULL} (the default) then \code{x} must be an
    ordinary matrix or array, or a dgCMatrix/lgCMatrix object, or any
    matrix-like or array-like object that supports coercion to SVT_SparseArray.

    If \code{dim} is provided then \code{x} can either be missing, a
    vector (atomic or list), or an array-like object. If missing, then an
    allzero SVT_SparseArray object will be constructed. Otherwise \code{x}
    will be used to fill the returned object (\code{length(x)} must be
    \code{<= prod(dim)}). Note that if \code{x} is an array-like object
    then its dimensions are ignored i.e. it's treated as a vector.
  }
  \item{dim}{
    \code{NULL} or the dimensions (supplied as an integer vector) of
    the SVT_SparseArray or SVT_SparseMatrix object to construct.
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

    Normally, the SVT_SparseArray object returned by the constructor
    function has the same \code{type()} as \code{x} but the user can use
    the \code{type} argument to request a different type. Note that doing:
    \preformatted{    svt <- SVT_SparseArray(x, type=type)}
    is equivalent to doing:
    \preformatted{    svt <- SVT_SparseArray(x)
    type(svt) <- type}
    but the former is more convenient and will generally be more efficient.

    Supported types are all R atomic types plus \code{"list"}.
  }
}

\details{
  SVT_SparseArray is a concrete subclass of the \link{SparseArray}
  virtual class. This makes SVT_SparseArray objects SparseArray derivatives.

  The nonzero data in a SVT_SparseArray object is stored in a \emph{Sparse
  Vector Tree}. We'll refer to this internal data representation as
  the \emph{SVT layout}. See the "SVT layout" section below for more
  information.

  The SVT layout is similar to the CSC layout (compressed, sparse,
  column-oriented format) used by CsparseMatrix derivatives from
  the \pkg{Matrix} package, like dgCMatrix or lgCMatrix objects,
  but with the following improvements:
  \itemize{
      \item The SVT layout supports sparse arrays of arbitrary dimensions.
      \item With the SVT layout, the sparse data can be of any type.
            Whereas CsparseMatrix derivatives only support sparse data
            of type \code{"double"} or \code{"logical"} at the moment.
      \item The SVT layout imposes no limit on the number of nonzero elements
            that can be stored. With dgCMatrix/lgCMatrix objects, this number
            must be < 2^31.
      \item Overall, the SVT layout allows more efficient operations on
            SVT_SparseArray objects.
  }
}

\value{
  An SVT_SparseArray or SVT_SparseMatrix object.
}

\section{SVT layout}{
  An SVT (Sparse Vector Tree) is a tree of depth N - 1 where N is the number
  of dimensions of the sparse array.

  The leaves in the tree can only be of two kinds: NULL or \emph{leaf vector}.
  Leaves that are leaf vectors can only be found at the deepest level in the
  tree (i.e. at depth N - 1). All leaves found at a lower depth must be NULLs.

  A leaf vector represents a sparse vector along the first dimension (a.k.a.
  innermost or fastest moving dimension) of the sparse array. It contains a
  collection of offset/value pairs sorted by strictly ascending offset.
  More precisely, a leaf vector is represented by an ordinary list of 2
  parallel dense vectors:
  \enumerate{
    \item nzvals: a vector (atomic or list) of nonzero values (zeros are
                  not allowed);
    \item nzoffs: an integer vector of offsets (i.e. 0-based positions).
  }
  The 1st vector determines the type of the leaf vector i.e. \code{"double"},
  \code{"integer"}, \code{"logical"}, etc...
  All the leaf vectors in the SVT must have the same type as the sparse array.

  It's useful to realize that a leaf vector simply represents a 1D SVT.

  In \pkg{SparseArray} 1.5.4 a new type of leaf vector was introduced called
  \emph{lacunar leaf}. A lacunar leaf is a non-empty leaf vector where the
  nzvals component is set to \code{NULL}. In this case the nonzero values are
  implicit: they're all considered to be equal to one.

  Examples:
  \itemize{
    \item An SVT_SparseArray object with 1 dimension has its nonzero data
          stored in an SVT of depth 0. Such SVT is represented by a
          single \emph{leaf vector}.

    \item An SVT_SparseArray object with 2 dimensions has its nonzero data
          stored in an SVT of depth 1. Such SVT is represented by a list of
          length the extend of the 2nd dimension (number of columns). Each
          list element is an SVT of depth 0 (as described above), or a NULL
          if the corresponding column is empty (i.e. has no nonzero data).

          For example, the nonzero data of an 8-column sparse matrix will be
          stored in an SVT that looks like this:
          \preformatted{
    .------------------list-of-length-8-----------------.
   /       /       /      |       |      \       \       \
  |       |       |       |       |       |       |       |
 leaf    leaf    NULL    leaf    leaf    leaf    leaf    NULL
vector  vector          vector  vector  vector  vector}

          The NULL leaves represent the empty columns (i.e. the columns
          with no nonzero elements).

    \item An SVT_SparseArray object with 3 dimensions has its nonzero data
          stored in an SVT of depth 2. Such SVT is represented by a list of
          length the extend of the 3rd dimension. Each list element must be
          an SVT of depth 1 (as described above) that stores the nonzero data
          of the corresponding 2D slice, or a NULL if the 2D slice is empty
          (i.e. has no nonzero data).

    \item And so on...
  }
}

\seealso{
  \itemize{
    \item The \link{SparseArray} class for the virtual parent class of
          COO_SparseArray and SVT_SparseArray.

    \item S4 classes \linkS4class{dgCMatrix} and \linkS4class{lgCMatrix}
          defined in the \pkg{Matrix} package, for the de facto standard
          of sparse matrix representations in the R ecosystem.

    \item Virtual class \linkS4class{CsparseMatrix} defined in the
          \pkg{Matrix} package for the parent class of all classes
          that use the "CSC layout".

    \item The \code{Matrix::\link[Matrix]{rsparsematrix}} function in
          the \pkg{Matrix} package.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
## ---------------------------------------------------------------------
## BASIC CONSTRUCTION
## ---------------------------------------------------------------------
SVT_SparseArray(dim=5:3)  # allzero object

SVT_SparseArray(dim=c(35000, 2e6), type="raw")  # allzero object

## Use a dgCMatrix object to fill the SVT_SparseArray object to construct:
x <- rsparsematrix(10, 16, density=0.1)  # random dgCMatrix object
SVT_SparseArray(x, dim=c(8, 5, 4))

svt1 <- SVT_SparseArray(dim=c(12, 5, 2))  # allzero object
svt1[cbind(11, 2:5, 2)] <- 22:25
svt1

svt2 <- SVT_SparseArray(dim=c(6, 4), type="integer",
                        dimnames=list(letters[1:6], LETTERS[1:4]))
svt2[c(1:2, 8, 10, 15:17, 24)] <- (1:8)*10L
svt2

## ---------------------------------------------------------------------
## CSC (Compressed Sparse Column) LAYOUT VS SVT LAYOUT
## ---------------------------------------------------------------------

## dgCMatrix objects from the Matrix package use the CSC layout:
dgcm2 <- as(svt2, "dgCMatrix")
dgcm2@x  # nonzero values
dgcm2@i  # row indices of the nonzero values
dgcm2@p  # breakpoints (0 followed by one breakpoint per column)

str(svt2)

m3 <- matrix(rpois(54e6, lambda=0.4), ncol=1200)

## Note that 'SparseArray(m3)' can also be used for this:
svt3 <- SVT_SparseArray(m3)
svt3

dgcm3 <- as(m3, "dgCMatrix")

## Compare type and memory footprint:
type(svt3)
object.size(svt3)
type(dgcm3)
object.size(dgcm3)

## Transpose:
system.time(svt <- t(t(svt3)))
system.time(dgcm <- t(t(dgcm3)))
identical(svt, svt3)
identical(dgcm, dgcm3)

## rbind():
m4 <- matrix(rpois(45e6, lambda=0.4), ncol=1200)
svt4 <- SVT_SparseArray(m4)
dgcm4 <- as(m4, "dgCMatrix")

system.time(rbind(svt3, svt4))
system.time(rbind(dgcm3, dgcm4))
}
\keyword{classes}
\keyword{methods}
