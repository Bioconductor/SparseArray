\name{NaArray-matrixStats}

\alias{NaArray-matrixStats}
\alias{NaArray_matrixStats}

\alias{colSums,NaArray-method}
\alias{rowSums,NaArray-method}

\alias{colMeans,NaArray-method}
\alias{rowMeans,NaArray-method}

\alias{colSums2,NaArray-method}
\alias{rowSums2,NaArray-method}

\alias{colMeans2,NaArray-method}
\alias{rowMeans2,NaArray-method}

\title{NaArray col/row summarization}

\description{
  EXPERIMENTAL!!!

  The \pkg{SparseArray} package provides memory-efficient col/row
  summarization methods (a.k.a. matrixStats methods) for \link{NaArray}
  objects, like \code{colSums()}, \code{rowSums()}, \code{colMeans()},
  \code{rowMeans()}, etc...

  Note that these are \emph{S4 generic functions} defined in the
  \pkg{MatrixGenerics} package, with methods for ordinary matrices
  defined in the \pkg{matrixStats} package. This man page documents
  the methods defined for \link{NaArray} objects.
}

\usage{
## N.B.: Showing ONLY the col*() methods (usage of row*() methods is
## the same):

\S4method{colSums}{NaArray}(x, na.rm=FALSE, dims=1)

\S4method{colMeans}{NaArray}(x, na.rm=FALSE, dims=1)

\S4method{colSums2}{NaArray}(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)

\S4method{colMeans2}{NaArray}(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
}

\arguments{
  \item{x}{
    An \link{NaMatrix} or \link{NaArray} object.
  }
  \item{rows, cols, ...}{
    Not supported.
  }
  \item{na.rm, useNames}{
    See man pages of the corresponding generics in the \pkg{MatrixGenerics}
    package (e.g. \code{?MatrixGenerics::\link[MatrixGenerics]{colSums2}})
    for a description of these arguments.
  }
  \item{dims}{
    See \code{?base::\link[base]{colSums}} for a description of this
    argument.
  }
}

\details{
  All these methods operate \emph{natively} on the \link{NaArray}
  internal representation, for maximum efficiency.

  Note that more col/row summarization methods might be added in the future.
}

\value{
  See man pages of the corresponding generics in the \pkg{MatrixGenerics}
  package (e.g. \code{?MatrixGenerics::\link[MatrixGenerics]{colMeans2}})
  for the value returned by these methods.
}

\note{
  Most \code{col*()} methods for \link{NaArray} objects are multithreaded.
  See \code{\link{set_SparseArray_nthread}} for how to control the number
  of threads.
}

\seealso{
  \itemize{
    \item \link{NaArray} objects.

    \item The man pages of the various generic functions defined
          in the \pkg{MatrixGenerics} package e.g.
          \code{MatrixGenerics::\link[MatrixGenerics]{colSums2}} etc...
  }
}

\examples{
# COMING SOON...
}
\keyword{array}
\keyword{methods}
\keyword{algebra}
\keyword{arith}
