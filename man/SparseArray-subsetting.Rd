\name{SparseArray-subsetting}

\alias{SparseArray-subsetting}
\alias{SparseArray_subsetting}

\alias{extract_sparse_array,COO_SparseArray-method}
\alias{extract_array,COO_SparseArray-method}
\alias{extract_sparse_array,SVT_SparseArray-method}

\alias{drop,SVT_SparseArray-method}
\alias{[,SVT_SparseArray,ANY,ANY,ANY-method}
\alias{extract_array,SVT_SparseArray-method}

\title{Subsetting a SparseArray object}

\description{
  Like ordinary arrays in R, \link{SparseArray} derivatives support
  subsetting via the single bracket operator (\code{[}).
}

\seealso{
  \itemize{
    \item \code{\link[base]{[}} in base R.

    \item \link{SparseArray} objects.

    \item \link[base]{array} objects in base R.
  }
}

\examples{
a <- array(0L, 5:3)
a[c(1:2, 8, 10, 15:17, 20, 24, 40, 56:60)] <- (1:15)*10L
svt <- SparseArray(a)
svt

svt[5:3, c(4L,2L,4L), 2:3]

svt[ , c(4L,2L,4L), 2:3]

svt[ , c(4L,2L,4L), -1]

svt[ , c(4L,2L,4L), 1L]

svt2 <- svt[ , c(4L,2L,4L), 1L, drop=FALSE]
svt2

## Ineffective dimensions can always be dropped as a separate step:
drop(svt2)

svt[ , c(4L,2L,4L), integer(0)]

dimnames(a) <- list(letters[1:5], NULL, LETTERS[1:3])
svt <- SparseArray(a)

svt[c("d", "a"), c(4L,2L,4L), "C"]

svt2 <- svt["e", c(4L,2L,4L), , drop=FALSE]
svt2

drop(svt2)

## Sanity checks:
svt2 <- svt[5:3, c(4L,2L,4L), 2:3]
a2   <- a  [5:3, c(4L,2L,4L), 2:3]
stopifnot(identical(as.array(svt2), a2), identical(svt2, SparseArray(a2)))
svt2 <- svt[ , c(4L,2L,4L), 2:3]
a2   <- a  [ , c(4L,2L,4L), 2:3]
stopifnot(identical(as.array(svt2), a2), identical(svt2, SparseArray(a2)))
svt2 <- svt[ , c(4L,2L,4L), -1]
a2   <- a  [ , c(4L,2L,4L), -1]
stopifnot(identical(as.array(svt2), a2), identical(svt2, SparseArray(a2)))
svt2 <- svt[ , c(4L,2L,4L), 1L]
a2   <- a  [ , c(4L,2L,4L), 1L]
stopifnot(identical(as.array(svt2), a2), identical(svt2, SparseArray(a2)))
svt2 <- svt[ , c(4L,2L,4L), 1L, drop=FALSE]
a2   <- a  [ , c(4L,2L,4L), 1L, drop=FALSE]
stopifnot(identical(as.array(svt2), a2), identical(svt2, SparseArray(a2)))
svt2 <- drop(svt2)
a2 <- drop(a2)
stopifnot(identical(as.array(svt2), a2), identical(svt2, SparseArray(a2)))
svt2 <- svt[ , c(4L,2L,4L), integer(0)]
a2   <- a  [ , c(4L,2L,4L), integer(0)]
stopifnot(identical(as.array(svt2), a2),
          identical(unname(svt2), unname(SparseArray(a2))))
svt2 <- svt[c("d", "a"), c(4L,2L,4L), "C"]
a2   <- a  [c("d", "a"), c(4L,2L,4L), "C"]
stopifnot(identical(as.array(svt2), a2), identical(svt2, SparseArray(a2)))
svt2 <- svt["e", c(4L,2L,4L), , drop=FALSE]
a2   <- a  ["e", c(4L,2L,4L), , drop=FALSE]
stopifnot(identical(as.array(svt2), a2), identical(svt2, SparseArray(a2)))
svt2 <- drop(svt2)
a2 <- drop(a2)
stopifnot(identical(as.array(svt2), a2), identical(svt2, SparseArray(a2)))
}
\keyword{methods}
\keyword{array}
