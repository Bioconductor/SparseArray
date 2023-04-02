\name{extract_sparse_array}

\alias{extract_sparse_array}
\alias{extract_sparse_array,ANY-method}

\title{extract_sparse_array}

\description{
  \code{extract_sparse_array()} is an internal generic function that
  is the workhorse behind the default \code{\link{read_block_as_sparse}()}
  method. It is not intended to be used directly by the end user.

  It is similar to the \code{\link[S4Arrays]{extract_array}()} internal
  generic function defined in the \pkg{S4Arrays} package, with the major
  difference that, in the case of \code{extract_sparse_array()}, the
  extracted array data is returned as a \link{SparseArray} object instead
  of an ordinay array.
}

\usage{
extract_sparse_array(x, index)

\S4method{extract_sparse_array}{ANY}(x, index)
}

\arguments{
  \item{x}{
    An array-like object for which \code{\link[S4Arrays]{is_sparse}(x)}
    is \code{TRUE}.
  }
  \item{index}{
    An unnamed list of integer vectors, one per dimension in \code{x}.
    Each vector is called a \emph{subscript} and can only contain
    positive integers that are valid 1-based indices along the corresponding
    dimension in \code{x}.

    Empty or missing subscripts are allowed. They must be represented
    by list elements set to \code{integer(0)} or \code{NULL}, respectively.

    The subscripts cannot contain NAs or non-positive values.

    Individual subscripts are NOT allowed to contain duplicated indices.
    This is an important difference with \code{\link[S4Arrays]{extract_array}}.
  }
}

\details{
  \code{extract_sparse_array()} should \emph{always} be called on
  an array-like object \code{x} for which \code{is_sparse(x)} is
  \code{TRUE}. Also it should \emph{never} be called with duplicated
  indices in the individual list elements of the \code{index} argument.

  For maximum efficiency, \code{extract_sparse_array()} methods should:
  \enumerate{
    \item NOT check that \code{\link[S4Arrays]{is_sparse}(x)} is \code{TRUE}.

    \item NOT check that the individual list elements in \code{index}
          contain no duplicated indices.

    \item NOT try to do anything with the dimnames on \code{x}.

    \item always operate natively on the sparse representation of the
          data in \code{x}, that is, they should never \emph{expand}
          it into a dense representation (e.g. with \emph{as.array()}).
  }

  Like for \code{\link[S4Arrays]{extract_array}()},
  \code{extract_sparse_array()} methods need to support empty or
  missing subscripts. For example, if \code{x} is an M x N matrix-like
  object for which \code{is_sparse(x)} is \code{TRUE}, then
  \code{extract_sparse_array(x, list(NULL, integer(0)))} must
  return an M x 0 \link{SparseArray} derivative, and
  \code{extract_sparse_array(x, list(integer(0), integer(0)))}
  a 0 x 0 \link{SparseArray} derivative.
}

\value{
  A \link{SparseArray} derivative (\link{COO_SparseArray} or
  \link{SVT_SparseArray}) of the same \code{type()} as \code{x}.
  For example, if \code{x} is an object representing an M x N sparse
  matrix of complex numbers (i.e. \code{type(x) == "complex"}), then
  \code{extract_sparse_array(x, list(NULL, 2L))} must return the 2nd column
  in \code{x} as an M x 1 \link{SparseArray} derivative of \code{type()}
  \code{"complex"}.
}

\seealso{
  \itemize{
    \item \code{\link[S4Arrays]{is_sparse}} in the \pkg{S4Arrays} package
          to check whether an object uses a sparse representation of the
          data or not.

    \item \link{SparseArray} objects.

    \item \code{S4Arrays::\link[S4Arrays]{type}} in the \pkg{S4Arrays}
          package to get the type of the elements of an array-like object.

    \item \code{\link{read_block_as_sparse}} to read array blocks as
          SparseArray objects.

    \item \code{\link[S4Arrays]{extract_array}} in the \pkg{S4Arrays} package.

    \item \linkS4class{dgCMatrix} objects implemented in the \pkg{Matrix}
          package.
  }
}

\examples{
extract_sparse_array
showMethods("extract_sparse_array")

## --- On a dgCMatrix object ---

m <- matrix(0L, nrow=6, ncol=4)
m[c(1:2, 8, 10, 15:17, 24)] <- (1:8)*10L
dgcm <- as(m, "dgCMatrix")
dgcm

extract_sparse_array(dgcm, list(3:6, NULL))
extract_sparse_array(dgcm, list(3:6, 2L))
extract_sparse_array(dgcm, list(3:6, integer(0)))

## --- On a SparseArray object ---

a <- array(0L, 5:3, dimnames=list(letters[1:5], NULL, LETTERS[1:3]))
a[c(1:2, 8, 10, 15:17, 20, 24, 40, 56:60)] <- (1:15)*10L
svt <- as(a, "SparseArray")
svt

extract_sparse_array(svt, list(NULL, 4:2, 1L))
extract_sparse_array(svt, list(NULL, 4:2, 2:3))
extract_sparse_array(svt, list(NULL, 4:2, integer(0)))
}
\keyword{internal}
\keyword{array}
\keyword{methods}
