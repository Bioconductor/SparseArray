\name{NaArray-subsetting}

\alias{NaArray-subsetting}
\alias{NaArray_subsetting}

\alias{tune_Array_dims,NaArray-method}

\alias{extract_na_array}
\alias{extract_na_array,NaArray-method}
\alias{extract_array,NaArray-method}

\title{Subsetting an NaArray object}

\description{
  EXPERIMENTAL!!!

  Like ordinary arrays in base R, \link{NaArray} objects support
  subsetting via the single bracket operator (\code{[}).
}

\seealso{
  \itemize{
    \item \code{\link[base]{drop}} in base R to drop the \emph{ineffective
          dimensions} of an array or array-like object.

    \item \code{\link[S4Arrays]{Lindex2Mindex}} in the \pkg{S4Arrays}
          package for how to convert an \emph{L-index} to an \emph{M-index}
          and vice-versa.

    \item \link{NaArray} objects.

    \item \code{\link[base]{[}} and ordinary \link[base]{array} objects
          in base R.
  }
}

\examples{
naa <- NaArray(dim=5:3)
naa[c(1:2, 8, 10, 15:17, 20, 24, 40, 56:60)] <- (1:15)*10L

## ---------------------------------------------------------------------
## N-dimensional subsetting
## ---------------------------------------------------------------------
naa[5:3, c(4,2,4), 2:3]
naa[ , c(4,2,4), 2:3]
naa[ , c(4,2,4), -1]
naa[ , c(4,2,4), 1]

naa2 <- naa[ , c(4,2,4), 1, drop=FALSE]
naa2

## Ineffective dimensions can always be dropped as a separate step:
drop(naa2)

naa[ , c(4,2,4), integer(0)]

dimnames(naa) <- list(letters[1:5], NULL, LETTERS[1:3])

naa[c("d", "a"), c(4,2,4), "C"]

naa2 <- naa["e", c(4,2,4), , drop=FALSE]
naa2

drop(naa2)

## ---------------------------------------------------------------------
## 1D-style subsetting (a.k.a. linear subsetting)
## ---------------------------------------------------------------------

## Using a numeric vector (L-index):
naa[c(60, 24, 21, 56)]

## Using a matrix subscript (M-index):
m <- rbind(c(5, 4, 3),
           c(4, 1, 2),
           c(1, 1, 2),
           c(1, 4, 3))
naa[m]

## See '?Lindex2Mindex' in the S4Arrays package for how to convert an
## L-index to an M-index and vice-versa.

## ---------------------------------------------------------------------
## Sanity checks
## ---------------------------------------------------------------------
a <- as.array(naa)
naa2 <- naa[5:3, c(4,2,4), 2:3]
a2   <- a  [5:3, c(4,2,4), 2:3]
stopifnot(identical(as.array(naa2), a2), identical(naa2, NaArray(a2)))
naa2 <- naa[ , c(4,2,4), 2:3]
a2   <- a  [ , c(4,2,4), 2:3]
stopifnot(identical(as.array(naa2), a2), identical(naa2, NaArray(a2)))
naa2 <- naa[ , c(4,2,4), -1]
a2   <- a  [ , c(4,2,4), -1]
stopifnot(identical(as.array(naa2), a2), identical(naa2, NaArray(a2)))
naa2 <- naa[ , c(4,2,4), 1]
a2   <- a  [ , c(4,2,4), 1]
stopifnot(identical(as.array(naa2), a2), identical(naa2, NaArray(a2)))
naa2 <- naa[ , c(4,2,4), 1, drop=FALSE]
a2   <- a  [ , c(4,2,4), 1, drop=FALSE]
stopifnot(identical(as.array(naa2), a2), identical(naa2, NaArray(a2)))
naa2 <- drop(naa2)
a2 <- drop(a2)
stopifnot(identical(as.array(naa2), a2), identical(naa2, NaArray(a2)))
naa2 <- naa[ , c(4,2,4), integer(0)]
a2   <- a  [ , c(4,2,4), integer(0)]
stopifnot(identical(as.array(naa2), a2),
          identical(unname(naa2), unname(NaArray(a2))))
naa2 <- naa[c("d", "a"), c(4,2,4), "C"]
a2   <- a  [c("d", "a"), c(4,2,4), "C"]
stopifnot(identical(as.array(naa2), a2), identical(naa2, NaArray(a2)))
naa2 <- naa["e", c(4,2,4), , drop=FALSE]
a2   <- a  ["e", c(4,2,4), , drop=FALSE]
stopifnot(identical(as.array(naa2), a2), identical(naa2, NaArray(a2)))
naa2 <- drop(naa2)
a2 <- drop(a2)
stopifnot(identical(as.array(naa2), a2), identical(naa2, NaArray(a2)))
stopifnot(identical(naa[c(60, 24, 21, 56)], naa[m]))
}
\keyword{array}
\keyword{methods}
