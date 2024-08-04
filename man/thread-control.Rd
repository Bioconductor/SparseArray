\name{thread-control}

\alias{thread-control}
\alias{thread_control}

\alias{get_SparseArray_nthread}
\alias{set_SparseArray_nthread}

\title{Number of threads used by SparseArray operations}

\description{
  Use \code{get_SparseArray_nthread} or \code{set_SparseArray_nthread}
  to get or set the number of threads to use by the multithreaded
  operations implemented in the \pkg{SparseArray} package.
}

\usage{
get_SparseArray_nthread()
set_SparseArray_nthread(nthread=NULL)
}

\arguments{
  \item{nthread}{
    The number of threads to use by multithreaded operations implemented
    in the \pkg{SparseArray} package.

    On systems where OpenMP is available, this must be \code{NULL}
    or an integer value >= 1. When \code{NULL} (the default),
    a "reasonable" value will be used that never exceeds one third
    of the number of logical cpus available on the machine.

    On systems where OpenMP is not available, the supplied \code{nthread}
    is ignored and \code{set_SparseArray_nthread()} is a no-op.
  }
}

\details{
  Multithreaded operations in the \pkg{SparseArray} package are implemented
  in C with OpenMP (\url{https://www.openmp.org/}).

  Note that OpenMP is not available on all systems. On systems where it's
  available, \code{get_SparseArray_nthread()} is guaranteed to return a
  value >= 1. On systems where it's not available (e.g. macOS),
  \code{get_SparseArray_nthread()} returns 0 and
  \code{set_SparseArray_nthread()} is a no-op.

  IMPORTANT: The portable way to disable multithreading is by calling
  \code{set_SparseArray_nthread(1)}, NOT \code{set_SparseArray_nthread(0)}
  (the latter returns an error on systems where OpenMP is available).
}

\value{
  \code{get_SparseArray_nthread()} returns an integer value >= 1 on systems
  where OpenMP is available, and 0 on systems where it's not.

  \code{set_SparseArray_nthread()} returns the \emph{previous} \code{nthread}
  value, that is, the value returned by \code{get_SparseArray_nthread()}
  before the call to \code{set_SparseArray_nthread()}. Note that the value
  is returned invisibly.
}

\seealso{
  \itemize{
    \item \link{SparseArray_matrixStats} for SparseArray col/row
          summarization methods.

    \item \link{SparseMatrix_mult} for SparseMatrix multiplication and
          cross-product.

    \item \link{SparseArray} objects.
  }
}

\examples{
get_SparseArray_nthread()

if (get_SparseArray_nthread() != 0) {  # multithreading is available
    svt1 <- poissonSparseMatrix(77000L, 15000L, density=0.01)

    ## 'user' time is typically N x 'elapsed' time where N is roughly the
    ## number of threads that was effectively used:
    system.time(cv1 <- colVars(svt1))

    svt2 <- poissonSparseMatrix(77000L, 300L, density=0.3) * 0.77
    system.time(cp12 <- crossprod(svt1, svt2))

    prev_nthread <- set_SparseArray_nthread(1)  # disable multithreading
    system.time(cv1 <- colVars(svt1))
    system.time(cp12 <- crossprod(svt1, svt2))

    ## Restore previous 'nthread' value:
    set_SparseArray_nthread(prev_nthread)
}
}

\keyword{utilities}
