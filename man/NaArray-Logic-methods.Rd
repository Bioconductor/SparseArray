\name{NaArray-Logic-methods}

\alias{NaArray-Logic-methods}
\alias{NaArray_Logic-methods}
\alias{NaArray-Logic}
\alias{NaArray_Logic}
\alias{NaArray-logic-methods}
\alias{NaArray_logic-methods}
\alias{NaArray-logic}
\alias{NaArray_logic}

\alias{!,NaArray-method}

\alias{Logic,NaArray,vector-method}
\alias{Logic,vector,NaArray-method}
\alias{Logic,NaArray,NaArray-method}
\alias{Logic,NaArray,array-method}
\alias{Logic,array,NaArray-method}
\alias{Logic,NaArray,SVT_SparseArray-method}
\alias{Logic,SVT_SparseArray,NaArray-method}

\title{'Logic' operations on NaArray objects}

\description{
  EXPERIMENTAL!!!

  \link{NaArray} objects support operations from the \code{Logic} group
  (i.e. \code{&} and \code{|}), as well as logical negation (\code{!}).
  See \code{?\link[methods]{S4groupGeneric}} in the \pkg{methods} package
  for more information about the \code{Logic} group generic.

  Note that in base R, \code{Logic} operations support input of \code{type()}
  \code{"logical"}, \code{"integer"}, \code{"double"}, or \code{"complex"}.
  However, the corresponding methods for \link{NaArray} objects only
  support objects of \code{type()} \code{"logical"} for now.
}

\details{
  Three forms of 'Logic' operations involving \link{NaArray} objects are
  supported:
  \enumerate{
    \item Between an \link{NaArray} object \code{naa} and a single
          logical value \code{y}: \preformatted{    naa op y
    y op naa}

    \item Between two \link{NaArray} objects \code{naa1} and \code{naa2}
          of same dimensions (a.k.a. \emph{conformable arrays}):
          \preformatted{    naa1 op naa2}

    \item Between an \link{NaArray} object \code{naa} and an
          \link{SVT_SparseArray} object \code{svt} of same dimensions
          (a.k.a. \emph{conformable arrays}):
          \preformatted{    naa op svt
    svt op naa}
          Note that, in this case, \code{|} returns an \link{NaArray}
          object but \code{&} returns an \link{SVT_SparseArray} object.
  }
}

\value{
  An \link{NaArray} object of \code{type()} \code{"logical"} and same
  dimensions as the input object(s), \emph{except} when \code{&} is used
  between an \link{NaArray} object and an \link{SVT_SparseArray} object
  in which case an \link{SVT_SparseArray} object is returned.
}

\seealso{
  \itemize{
    \item \code{\link[methods]{S4groupGeneric}} in the \pkg{methods} package.

    \item \link{NaArray} objects.

    \item \link{SVT_SparseArray} objects.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
nam1 <- NaArray(dim=c(15, 6))
nam1[cbind(1:15, 2)] <- c(TRUE, FALSE, NA)
nam1[cbind(1:15, 5)] <- c(TRUE, NA, NA, FALSE, TRUE)
nam1

!nam1
nam1 & NA  # replaces all TRUE's with NA's
nam1 | NA  # replaces all FALSE's with NA's

nam2 <- NaArray(dim=dim(nam1))
nam2[c(2, 6, 12:17, 22:33, 55, 59:62, 90)] <- c(TRUE, NA)
nam2

nam1 & nam2
nam1 | nam2

## Sanity checks:
m1 <- as.matrix(nam1)
m2 <- as.matrix(nam2)
stopifnot(
  identical(as.matrix(!nam1), !m1),
  identical(as.matrix(nam1 & NA), m1 & NA),
  identical(as.matrix(nam1 | NA), m1 | NA),
  identical(as.matrix(nam1 & nam2), m1 & m2),
  identical(as.matrix(nam1 | nam2), m1 | m2)
)
}
\keyword{array}
\keyword{methods}
\keyword{algebra}
