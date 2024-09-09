\name{SparseArray-summarization}

\alias{SparseArray-summarization}
\alias{SparseArray_summarization}

\alias{anyNA,SparseArray-method}

\alias{range.COO_SparseArray}
\alias{range,COO_SparseArray-method}
\alias{range.SVT_SparseArray}
\alias{range,SVT_SparseArray-method}

\alias{mean.SparseArray}
\alias{mean,SparseArray-method}

\alias{var,SparseArray,ANY-method}
\alias{sd,SparseArray-method}

\title{SparseArray summarization methods}

\description{
  The \pkg{SparseArray} package provides memory-efficient summarization
  methods for \link{SparseArray} objects. The following methods are
  supported at the moment: \code{anyNA()}, \code{any()}, \code{all()},
  \code{min()}, \code{max()}, \code{range()}, \code{sum()}, \code{prod()},
  \code{mean()}, \code{var()}, \code{sd()}.

  More might be added in the future.

  Note that these are \emph{S4 generic functions} defined in base R
  and in the \pkg{BiocGenerics} package, with default methods defined
  in base R. This man page documents the methods defined for
  \link{SparseArray} objects.
}

\details{
  All these methods operate \emph{natively} on the \link{COO_SparseArray}
  or \link{SVT_SparseArray} representation, for maximum efficiency.
}

\value{
  See man pages of the corresponding default methods in the
  \pkg{base} package (e.g. \code{?base::\link[base]{range}},
  \code{?base::\link[base]{mean}}, etc...) for the value returned
  by these methods.
}

\seealso{
  \itemize{
    \item \link{SparseArray} objects.

    \item The man pages of the various default methods defined in
          the \pkg{base} package e.g. \code{base::\link[base]{range}},
	  \code{base::\link[base]{mean}}, \code{base::\link[base]{anyNA}},
	  etc...
  }
}

\examples{
svt0 <- SVT_SparseArray(dim=c(4, 5, 2))
svt0[c(1:2, 8, 10, 15:17, 24:26, 28, 40)] <- (1:12)*10L
svt0[4, 3, 1] <- NA
svt0

anyNA(svt0)

range(svt0)

range(svt0, na.rm=TRUE)

sum(svt0, na.rm=TRUE)

sd(svt0, na.rm=TRUE)

## Sanity checks:
a0 <- as.array(svt0)
stopifnot(
  identical(anyNA(svt0), anyNA(a0)),
  identical(range(svt0), range(a0)),
  identical(range(svt0, na.rm=TRUE), range(a0, na.rm=TRUE)),
  identical(sum(svt0), sum(a0)),
  identical(sum(svt0, na.rm=TRUE), sum(a0, na.rm=TRUE)),
  all.equal(sd(svt0, na.rm=TRUE), sd(a0, na.rm=TRUE))
)
}
\keyword{array}
\keyword{methods}
\keyword{algebra}
\keyword{arith}
