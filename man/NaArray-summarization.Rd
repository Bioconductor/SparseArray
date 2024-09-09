\name{NaArray-summarization}

\alias{NaArray-summarization}
\alias{NaArray_summarization}

\alias{anyNA,NaArray-method}

\alias{range.NaArray}
\alias{range,NaArray-method}

\alias{mean.NaArray}
\alias{mean,NaArray-method}

\alias{var,NaArray,ANY-method}
\alias{sd,NaArray-method}

\title{NaArray summarization methods}

\description{
  EXPERIMENTAL!!!

  The \pkg{SparseArray} package provides memory-efficient summarization
  methods for \link{NaArray} objects. The following methods are
  supported at the moment: \code{anyNA()}, \code{any()}, \code{all()},
  \code{min()}, \code{max()}, \code{range()}, \code{sum()}, \code{prod()},
  \code{mean()}, \code{var()}, \code{sd()}.

  More might be added in the future.

  Note that these are \emph{S4 generic functions} defined in base R
  and in the \pkg{BiocGenerics} package, with default methods defined
  in base R. This man page documents the methods defined for
  \link{NaArray} objects.
}

\details{
  All these methods operate \emph{natively} on the \link{NaArray}
  representation, for maximum efficiency.
}

\value{
  See man pages of the corresponding default methods in the
  \pkg{base} package (e.g. \code{?base::\link[base]{range}},
  \code{?base::\link[base]{mean}}, etc...) for the value returned
  by these methods.
}

\seealso{
  \itemize{
    \item \link{NaArray} objects.

    \item The man pages of the various default methods defined in
          the \pkg{base} package e.g. \code{base::\link[base]{range}},
	  \code{base::\link[base]{mean}}, \code{base::\link[base]{anyNA}},
	  etc...
  }
}

\examples{
naa <- NaArray(dim=c(4, 5, 2))
naa[c(1:2, 8, 10, 15:17, 24:26, 28, 40)] <- (1:12)*10L
naa

anyNA(naa)

range(naa, na.rm=TRUE)

sum(naa, na.rm=TRUE)

sd(naa, na.rm=TRUE)

## Sanity checks:
a0 <- as.array(naa)
stopifnot(
  identical(anyNA(naa), anyNA(a0)),
  identical(range(naa), range(a0)),
  identical(range(naa, na.rm=TRUE), range(a0, na.rm=TRUE)),
  identical(sum(naa), sum(a0)),
  identical(sum(naa, na.rm=TRUE), sum(a0, na.rm=TRUE)),
  all.equal(sd(naa, na.rm=TRUE), sd(a0, na.rm=TRUE))
)
}
\keyword{array}
\keyword{methods}
\keyword{algebra}
\keyword{arith}
