\name{SparseArray-summarization}

\alias{SparseArray-summarization}
\alias{SparseArray_summarization}

\alias{range.COO_SparseArray}
\alias{range,COO_SparseArray-method}
\alias{range.SVT_SparseArray}
\alias{range,SVT_SparseArray-method}

\alias{mean,SparseArray-method}

\alias{anyNA,COO_SparseArray-method}
\alias{anyNA,SVT_SparseArray-method}

\alias{var,SparseArray,ANY-method}

\title{SparseArray summarization methods}

\description{
  The \pkg{SparseArray} package provides memory-efficient summarization
  methods for \link{SparseArray} objects. The following methods are
  supported at the moment: \code{min()}, \code{max()}, \code{range()},
  \code{sum()}, \code{prod()}, \code{any()}, \code{all()},
  \code{mean()}, \code{anyNA()}.
  More will be added in the future.

  Note that these are \emph{S4 generic functions} defined in base R
  and in the \pkg{BiocGenerics} package, with default methods defined
  in base R. This man page documents the methods defined for
  \link{SparseArray} objects.
}

\details{
  All these methods operate \emph{natively} on the \link{COO_SparseMatrix}
  or \link{SVT_SparseMatrix} representation, for maximum efficiency.

  Note that more summarization methods will be added in the future.
}

\value{
  See man pages for the corresponding default methods in the
  \pkg{base} package (e.g. \code{?base::\link[base]{range}},
  \code{?base::\link[base]{mean}}, etc...) for the value returned
  by these methods.
}

\seealso{
  \itemize{
    \item \link{SparseArray} objects.

    \item The man pages for the various default methods defined in
          the \pkg{base} package e.g. \code{base::\link[base]{range}},
	  \code{base::\link[base]{mean}}, \code{base::\link[base]{anyNA}},
	  etc...
  }
}

\examples{
m0 <- matrix(0L, nrow=6, ncol=4, dimnames=list(letters[1:6], LETTERS[1:4]))
m0[c(1:2, 8, 10, 15:17, 24)] <- (1:8)*10L
m0["e", "B"] <- NA
svt0 <- as(m0, "SVT_SparseMatrix")
svt0

range(svt0)
range(svt0, na.rm=TRUE)

sum(svt0)
sum(svt0, na.rm=TRUE)

anyNA(svt0)

## Sanity checks:
stopifnot(
  identical(range(svt0), range(m0)),
  identical(range(svt0, na.rm=TRUE), range(m0, na.rm=TRUE)),
  identical(sum(svt0), sum(m0)),
  identical(sum(svt0, na.rm=TRUE), sum(m0, na.rm=TRUE)),
  identical(anyNA(svt0), anyNA(m0))
)
}
\keyword{array}
\keyword{methods}
\keyword{algebra}
\keyword{arith}
