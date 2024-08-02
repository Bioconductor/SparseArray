\name{NaArray-abind}

\alias{NaArray-abind}
\alias{NaArray_abind}
\alias{NaArray-combine}
\alias{NaArray_combine}

\alias{abind,NaArray-method}

\alias{rbind,NaArray-method}
\alias{cbind,NaArray-method}
\alias{bindROWS,NaArray-method}

\title{Combine multidimensional NaArray objects}

\description{
  Like ordinary matrices and arrays in base R, \link{NaMatrix} objects
  can be combined by rows or columns, with \code{rbind()} or \code{cbind()},
  and multidimensional \link{NaArray} objects can be bound along
  any dimension with \code{abind()}.

  Note that \code{arbind()} can also be used to combine the objects along
  their first dimension, and \code{acbind()} can be used to combine them
  along their second dimension.
}

\seealso{
  \itemize{
    \item \code{\link[base]{cbind}} in base R.

    \item \code{\link[S4Arrays]{abind}} in the \pkg{S4Arrays} package.

    \item \link{NaArray} objects.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
# COMING SOON...
}
\keyword{array}
\keyword{methods}
\keyword{manip}
