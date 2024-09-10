\name{NaArray-matrixStats}

\alias{NaArray-matrixStats}
\alias{NaArray_matrixStats}

\alias{colAnyNAs,NaArray-method}

\alias{colAnys,NaArray-method}

\alias{colAlls,NaArray-method}

\alias{colMins,NaArray-method}

\alias{colMaxs,NaArray-method}

\alias{colRanges,NaArray-method}

\alias{colSums,NaArray-method}
\alias{rowSums,NaArray-method}

\alias{colProds,NaArray-method}

\alias{colMeans,NaArray-method}

\alias{colSums2,NaArray-method}
\alias{rowSums2,NaArray-method}

\alias{colMeans2,NaArray-method}

\alias{colVars,NaArray-method}

\alias{colSds,NaArray-method}

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

  WORK IN PROGRESS: All the col*() methods listed below are supported.
  However, among the row*() methods, only rowSums() and rowSums2() are
  supported at the moment.
}

\usage{
## N.B.: Showing ONLY the col*() methods (usage of row*() methods is
## the same):

\S4method{colAnyNAs}{NaArray}(x, rows=NULL, cols=NULL, dims=1, ..., useNames=NA)

\S4method{colAnys}{NaArray}(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)

\S4method{colAlls}{NaArray}(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)

\S4method{colMins}{NaArray}(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)

\S4method{colMaxs}{NaArray}(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)

\S4method{colRanges}{NaArray}(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)

\S4method{colSums}{NaArray}(x, na.rm=FALSE, dims=1)

\S4method{colProds}{NaArray}(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)

\S4method{colMeans}{NaArray}(x, na.rm=FALSE, dims=1)

\S4method{colSums2}{NaArray}(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)

\S4method{colMeans2}{NaArray}(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)

\S4method{colVars}{NaArray}(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL, dims=1,
           ..., useNames=NA)

\S4method{colSds}{NaArray}(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL, dims=1,
          ..., useNames=NA)
}

\arguments{
  \item{x}{
    An \link{NaMatrix} or \link{NaArray} object.
  }
  \item{rows, cols, ...}{
    Not supported.
  }
  \item{na.rm, useNames, center}{
    See man pages of the corresponding generics in the \pkg{MatrixGenerics}
    package (e.g. \code{?MatrixGenerics::\link[MatrixGenerics]{colVars}})
    for a description of these arguments.

    Note that, unlike the methods for ordinary matrices defined in
    the \pkg{matrixStats} package, the \code{center} argument of the
    \code{colVars()}, \code{rowVars()}, \code{colSds()}, and \code{rowSds()}
    methods for \link{SparseArray} objects can only be a
    \emph{single value} (or a \code{NULL}). In particular, if \code{x}
    has more than one column, then \code{center} cannot be a vector
    with one value per column in \code{x}.
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
  package (e.g. \code{?MatrixGenerics::\link[MatrixGenerics]{colRanges}})
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
          \code{MatrixGenerics::\link[MatrixGenerics]{colVars}} etc...
  }
}

\examples{
# COMING SOON...
}
\keyword{array}
\keyword{methods}
\keyword{algebra}
\keyword{arith}
