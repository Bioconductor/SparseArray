\name{matrixStats-methods}

\alias{matrixStats-methods}

\alias{colSums}
\alias{colSums,SVT_SparseMatrix-method}
\alias{rowSums}
\alias{rowSums,SVT_SparseMatrix-method}

\alias{colMeans}
\alias{colMeans,SVT_SparseMatrix-method}
\alias{rowMeans}
\alias{rowMeans,SVT_SparseMatrix-method}

\alias{colMins}
\alias{colMins,SVT_SparseMatrix-method}
\alias{rowMins}
\alias{rowMins,SVT_SparseMatrix-method}

\alias{colMaxs}
\alias{colMaxs,SVT_SparseMatrix-method}
\alias{rowMaxs}
\alias{rowMaxs,SVT_SparseMatrix-method}

\alias{colRanges}
\alias{colRanges,SVT_SparseMatrix-method}
\alias{rowRanges}
\alias{rowRanges,SVT_SparseMatrix-method}

\alias{colMedians}
\alias{colMedians,SVT_SparseMatrix-method}
\alias{rowMedians}
\alias{rowMedians,SVT_SparseMatrix-method}

\alias{colVars}
\alias{colVars,SVT_SparseMatrix-method}
\alias{rowVars}
\alias{rowVars,SVT_SparseMatrix-method}

\alias{colSds}
\alias{colSds,SVT_SparseMatrix-method}
\alias{rowSds}
\alias{rowSds,SVT_SparseMatrix-method}

\title{SparseMatrix row/col summarization methods}

\description{
  The \pkg{SparseArray} package provides memory-efficient row/col
  summarization methods for \link{SparseMatrix} objects, like
  \code{rowSums()}, \code{colSums()}, \code{rowMedians()},
  \code{colMedians()}, \code{rowVars()}, \code{colVars()}, etc...

  Note that these are \emph{S4 generic functions} defined in the
  \pkg{MatrixGenerics} package, with methods for ordinary matrices
  defined in the \pkg{matrixStats} package. This man page documents
  the methods for \link{SVT_SparseMatrix} objects.
}

\usage{
\S4method{rowSums}{SVT_SparseMatrix}(x, na.rm=FALSE, dims=1)
\S4method{colSums}{SVT_SparseMatrix}(x, na.rm=FALSE, dims=1)

\S4method{rowMeans}{SVT_SparseMatrix}(x, na.rm=FALSE, dims=1)
\S4method{colMeans}{SVT_SparseMatrix}(x, na.rm=FALSE, dims=1)

\S4method{rowMins}{SVT_SparseMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
\S4method{colMins}{SVT_SparseMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)

\S4method{rowMaxs}{SVT_SparseMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
\S4method{colMaxs}{SVT_SparseMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)

\S4method{rowRanges}{SVT_SparseMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
\S4method{colRanges}{SVT_SparseMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)

\S4method{rowMedians}{SVT_SparseMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
\S4method{colMedians}{SVT_SparseMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)

\S4method{rowVars}{SVT_SparseMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL,
           ..., useNames=NA)
\S4method{colVars}{SVT_SparseMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL,
           ..., useNames=NA)

\S4method{rowSds}{SVT_SparseMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL,
          ..., useNames=NA)
\S4method{colSds}{SVT_SparseMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL,
          ..., useNames=NA)
}

\arguments{
  \item{x}{
    An \link{SVT_SparseMatrix} object.
  }
  \item{rows, cols, ...}{
    Not supported.
  }
  \item{na.rm, dims, center, useNames}{
    See man pages for the corresponding generics in the \pkg{MatrixGenerics}
    package (e.g. \code{?MatrixGenerics::\link[MatrixGenerics]{rowVars}})
    for a description of these arguments.
  }
}

\details{
  All these methods operate \emph{natively} on the \link{SVT_SparseMatrix}
  representation, and so are fast and memory efficient.
}

\value{
  See man pages for the corresponding generics in the \pkg{MatrixGenerics}
  package (e.g. \code{?MatrixGenerics::\link[MatrixGenerics]{colRanges}})
  for the value returned by these methods.
}

\seealso{
  \itemize{
    \item \link{SVT_SparseMatrix} objects.

    \item The man pages for the various generic functions defined in
          the \pkg{MatrixGenerics} package e.g.
          \code{MatrixGenerics::\link[MatrixGenerics]{rowVars}} etc...
  }
}

\examples{
m0 <- matrix(0L, nrow=6, ncol=4, dimnames=list(letters[1:6], LETTERS[1:4]))
m0[c(1:2, 8, 10, 15:17, 24)] <- (1:8)*10L
m0["e", "B"] <- NA
svt0 <- as(m0, "SVT_SparseMatrix")
svt0

rowSums(svt0)
rowSums(svt0, na.rm=TRUE)

colSums(svt0)
colSums(svt0, na.rm=TRUE)

rowMeans(svt0)
rowMeans(svt0, na.rm=TRUE)

colRanges(svt0)
colRanges(svt0, useNames=TRUE)
colRanges(svt0, na.rm=TRUE)
colRanges(svt0, na.rm=TRUE, useNames=TRUE)

colVars(svt0)
colVars(svt0, useNames=TRUE)

## Sanity checks:
stopifnot(
  identical(rowSums(svt0), rowSums(m0)),
  identical(rowSums(svt0, na.rm=TRUE), rowSums(m0, na.rm=TRUE)),
  identical(colSums(svt0), colSums(m0)),
  identical(colSums(svt0, na.rm=TRUE), colSums(m0, na.rm=TRUE)),
  identical(rowMeans(svt0), rowMeans(m0)),
  identical(rowMeans(svt0, na.rm=TRUE), rowMeans(m0, na.rm=TRUE)),
  identical(colRanges(svt0), colRanges(m0)),
  identical(colRanges(svt0, useNames=TRUE), colRanges(m0, useNames=TRUE)),
  identical(colRanges(svt0, na.rm=TRUE), colRanges(m0, na.rm=TRUE)),
  identical(colRanges(svt0, na.rm=TRUE, useNames=TRUE),
            colRanges(m0, na.rm=TRUE, useNames=TRUE)),
  identical(colVars(svt0), colVars(m0)),
  identical(colVars(svt0, na.rm=TRUE), colVars(m0, na.rm=TRUE))
)
}
\keyword{methods}
\keyword{algebra}
\keyword{arith}
\keyword{array}
