\name{read_block_as_sparse}

\alias{read_block_as_sparse}
\alias{read_block_as_sparse,ANY-method}

\title{read_block_as_sparse}

\description{
  \code{read_block_as_sparse()} is an internal generic function used by
  \code{S4Arrays::\link[S4Arrays]{read_block}()} when \code{is_sparse(x)}
  is \code{TRUE}.
}

\usage{
read_block_as_sparse(x, viewport)

\S4method{read_block_as_sparse}{ANY}(x, viewport)
}

\arguments{
  \item{x}{
    An array-like object for which \code{\link[S4Arrays]{is_sparse}(x)}
    is \code{TRUE}.
  }
  \item{viewport}{
    An \link[S4Arrays]{ArrayViewport} object compatible with \code{x},
    that is, such that \code{refdim(viewport)} is identical to \code{dim(x)}.
  }
}

\details{
  Like \code{\link[S4Arrays]{read_block_as_dense}()} in the \pkg{S4Arrays}
  package, \code{read_block_as_sparse()} is not meant to be called directly
  by the end user. The end user should always call the higher-level
  user-facing \code{\link[S4Arrays]{read_block}()} function instead.
  See \code{?\link[S4Arrays]{read_block}} in the \pkg{S4Arrays} package
  for more information.

  Also, like \code{\link{extract_sparse_array}()},
  \code{read_block_as_sparse()} should \emph{always} be called on
  an array-like object \code{x} for which \code{is_sparse(x)} is
  \code{TRUE}.

  For maximum efficiency, \code{read_block_as_sparse()} methods should:
  \enumerate{
    \item NOT check that \code{\link[S4Arrays]{is_sparse}(x)} is \code{TRUE}.

    \item NOT try to do anything with the dimnames on \code{x}
          (\code{\link[S4Arrays]{read_block}()} takes care of that).

    \item always operate natively on the sparse representation of the
          data in \code{x}, that is, they should never \emph{expand}
          it into a dense representation (e.g. with \emph{as.array()}).
  }
}

\value{
  A block of data as a \link{SparseArray} derivative (\link{COO_SparseArray}
  or \link{SVT_SparseArray}) of the same \code{type()} as \code{x}.
}

\seealso{
  \itemize{
    \item \code{\link[S4Arrays]{read_block}} in the \pkg{S4Arrays} package
          for the higher-level user-facing function for reading array blocks.

    \item \link[S4Arrays]{ArrayGrid} in the \pkg{S4Arrays} package for
          ArrayGrid and ArrayViewport objects.

    \item \code{\link[S4Arrays]{is_sparse}} in the \pkg{S4Arrays} package
          to check whether an object uses a sparse representation of the
	  data or not.

    \item \link{SparseArray} objects.

    \item \code{S4Arrays::\link[S4Arrays]{type}} in the \pkg{S4Arrays}
          package to get the type of the elements of an array-like object.

    \item \code{\link{extract_sparse_array}} for the workhorse behind
          the default \code{read_block_as_sparse()} method.

    \item \linkS4class{dgCMatrix} objects implemented in the \pkg{Matrix}
          package.
  }
}

\keyword{internal}
\keyword{array}
\keyword{methods}
