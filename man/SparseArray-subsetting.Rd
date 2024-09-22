\name{SparseArray-subsetting}

\alias{SparseArray-subsetting}
\alias{SparseArray_subsetting}

\alias{extract_sparse_array,COO_SparseArray-method}
\alias{extract_array,COO_SparseArray-method}
\alias{extract_sparse_array,SVT_SparseArray-method}

\alias{extract_array,SVT_SparseArray-method}

\title{Subsetting a SparseArray object}

\description{
  Like ordinary arrays in base R, \link{SparseArray} derivatives support
  subsetting via the single bracket operator (\code{[}).
}

\seealso{
  \itemize{
    \item \code{\link[base]{drop}} in base R to drop the \emph{ineffective
          dimensions} of an array or array-like object.

    \item \code{\link[S4Arrays]{Lindex2Mindex}} in the \pkg{S4Arrays}
          package for how to convert an \emph{L-index} to an \emph{M-index}
          and vice-versa.

    \item \link{SparseArray} objects.

    \item \code{\link[base]{[}} and ordinary \link[base]{array} objects
          in base R.
  }
}

\examples{
a <- array(0L, dim=5:3)
a[c(1:2, 8, 10, 15:17, 20, 24, 40, 56:60)] <- (1:15)*10L
svt <- SparseArray(a)
svt

## ---------------------------------------------------------------------
## N-dimensional subsetting
## ---------------------------------------------------------------------
svt[5:3, c(4,2,4), 2:3]
svt[ , c(4,2,4), 2:3]
svt[ , c(4,2,4), -1]
svt[ , c(4,2,4), 1]

svt2 <- svt[ , c(4,2,4), 1, drop=FALSE]
svt2

## Ineffective dimensions can always be dropped as a separate step:
drop(svt2)

svt[ , c(4,2,4), integer(0)]

dimnames(a) <- list(letters[1:5], NULL, LETTERS[1:3])
svt <- SparseArray(a)

svt[c("d", "a"), c(4,2,4), "C"]

svt2 <- svt["e", c(4,2,4), , drop=FALSE]
svt2

drop(svt2)

## ---------------------------------------------------------------------
## 1D-style subsetting (a.k.a. linear subsetting)
## ---------------------------------------------------------------------

## Using a numeric vector (L-index):
svt[c(60, 24, 21, 56)]

## Using a matrix subscript (M-index):
m <- rbind(c(5, 4, 3),
           c(4, 1, 2),
           c(1, 1, 2),
           c(1, 4, 3))
svt[m]

## See '?Lindex2Mindex' in the S4Arrays package for how to convert an
## L-index to an M-index and vice-versa.

## ---------------------------------------------------------------------
## Sanity checks
## ---------------------------------------------------------------------
svt2 <- svt[5:3, c(4,2,4), 2:3]
a2   <- a  [5:3, c(4,2,4), 2:3]
stopifnot(identical(as.array(svt2), a2), identical(svt2, SparseArray(a2)))
svt2 <- svt[ , c(4,2,4), 2:3]
a2   <- a  [ , c(4,2,4), 2:3]
stopifnot(identical(as.array(svt2), a2), identical(svt2, SparseArray(a2)))
svt2 <- svt[ , c(4,2,4), -1]
a2   <- a  [ , c(4,2,4), -1]
stopifnot(identical(as.array(svt2), a2), identical(svt2, SparseArray(a2)))
svt2 <- svt[ , c(4,2,4), 1]
a2   <- a  [ , c(4,2,4), 1]
stopifnot(identical(as.array(svt2), a2), identical(svt2, SparseArray(a2)))
svt2 <- svt[ , c(4,2,4), 1, drop=FALSE]
a2   <- a  [ , c(4,2,4), 1, drop=FALSE]
stopifnot(identical(as.array(svt2), a2), identical(svt2, SparseArray(a2)))
svt2 <- drop(svt2)
a2 <- drop(a2)
stopifnot(identical(as.array(svt2), a2), identical(svt2, SparseArray(a2)))
svt2 <- svt[ , c(4,2,4), integer(0)]
a2   <- a  [ , c(4,2,4), integer(0)]
stopifnot(identical(as.array(svt2), a2),
          identical(unname(svt2), unname(SparseArray(a2))))
svt2 <- svt[c("d", "a"), c(4,2,4), "C"]
a2   <- a  [c("d", "a"), c(4,2,4), "C"]
stopifnot(identical(as.array(svt2), a2), identical(svt2, SparseArray(a2)))
svt2 <- svt["e", c(4,2,4), , drop=FALSE]
a2   <- a  ["e", c(4,2,4), , drop=FALSE]
stopifnot(identical(as.array(svt2), a2), identical(svt2, SparseArray(a2)))
svt2 <- drop(svt2)
a2 <- drop(a2)
stopifnot(identical(as.array(svt2), a2), identical(svt2, SparseArray(a2)))
stopifnot(identical(svt[c(60, 24, 21, 56)], svt[m]))
}
\keyword{array}
\keyword{methods}
