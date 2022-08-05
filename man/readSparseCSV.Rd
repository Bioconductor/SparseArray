\name{readSparseCSV}

\alias{writeSparseCSV}
\alias{readSparseCSV}
\alias{readSparseTable}

\title{Write/read a sparse matrix from/to a CSV file}

\description{
  Write or read a sparse matrix from or to a CSV (comma-separated values) file.
}

\usage{
writeSparseCSV(x, filepath, sep=",", transpose=FALSE, write.zeros=FALSE,
                  chunknrow=250)

readSparseCSV(filepath, sep=",", transpose=FALSE)
}

\arguments{
  \item{x}{
    A matrix-like object, typically sparse. IMPORTANT: The object must have
    rownames and colnames! These will be written to the file.

    Another requirement is that the object must be subsettable. More
    precisely: it must support 2D-style subsetting of the kind \code{x[i, ]}
    and \code{x[ , j]} where \code{i} and \code{j} are integer vectors
    of valid row and column indices.
  }
  \item{filepath}{
    The path (as a single string) to the file where to write the
    matrix-like object or to read it from. Compressed files are
    supported.

    If \code{""}, \code{writeSparseCSV()} will write the data to the
    standard output connection.

    Note that \code{filepath} can also be a connection.
  }
  \item{sep}{
    The field separator character. Values on each line of the
    file are separated by this character.
  }
  \item{transpose}{
    \code{TRUE} or \code{FALSE}. By default, rows in the matrix-like object
    correspond to lines in the CSV file. Set \code{transpose} to \code{TRUE}
    to transpose the matrix-like object on-the-fly, that is, to have its
    columns written to or read from the lines in the CSV file.

    Note that using \code{transpose=TRUE} is semantically equivalent to
    calling \code{t()} on the object before writting it or after reading
    it, but it will tend to be more efficient. Also it will work even if
    \code{x} does not support \code{t()} (not all matrix-like objects are
    guaranteed to be transposables).
  }
  \item{write.zeros}{
    \code{TRUE} or \code{FALSE}. By default, the zero values in \code{x}
    are not written to the file. Set \code{write.zeros} to \code{TRUE} to
    write them.
  }
  \item{chunknrow}{
    \code{writeSparseCSV()} uses a block-processing strategy to try to
    speed up things. By default blocks of 250 rows (or columns if
    \code{transpose=TRUE}) are used. In our experience trying to increase
    this (e.g. to 500 or more) will generally not produce siginificant
    benefits while it will increase memory usage, so use carefully.
  }
}

\value{
  \code{writeSparseCSV} returns an invisible \code{NULL}.

  \code{readSparseCSV} returns a \link{SparseMatrix} object of class
  \link{SVT_SparseMatrix}.
}

\seealso{
  \itemize{
    \item \link{SparseArray} objects.

    \item \linkS4class{dgCMatrix} objects in the \pkg{Matrix} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## writeSparseCSV()
## ---------------------------------------------------------------------

## Prepare toy matrix 'm0':
rownames0 <- LETTERS[1:6]
colnames0 <- letters[1:4]
m0 <- matrix(0L, nrow=length(rownames0), ncol=length(colnames0),
                 dimnames=list(rownames0, colnames0))
m0[c(1:2, 8, 10, 15:17, 24)] <- (1:8)*10L
m0

## writeSparseCSV():
writeSparseCSV(m0, filepath="", sep="\t")
writeSparseCSV(m0, filepath="", sep="\t", write.zeros=TRUE)
writeSparseCSV(m0, filepath="", sep="\t", transpose=TRUE)

## Note that writeSparseCSV() will automatically (and silently) coerce
## non-integer values to integer by passing them thru as.integer().

## Example where type(x) is "double":
m1 <- m0 * runif(length(m0))
m1
type(m1)
writeSparseCSV(m1, filepath="", sep="\t")

## Example where type(x) is "logical":
writeSparseCSV(m0 != 0, filepath="", sep="\t")

## Example where type(x) is "raw":
m2 <- m0
type(m2) <- "raw"
m2
writeSparseCSV(m2, filepath="", sep="\t")

## ---------------------------------------------------------------------
## readSparseCSV()
## ---------------------------------------------------------------------

outfile <- tempfile()
writeSparseCSV(m0, outfile)

svt1 <- readSparseCSV(outfile)
svt1

svt2 <- readSparseCSV(outfile, transpose=TRUE)
svt2

## If you need the sparse data as a dgCMatrix object, just coerce the
## returned object:
as(svt1, "dgCMatrix")
as(svt2, "dgCMatrix")

## Sanity checks:
stopifnot(identical(m0, as.matrix(svt1)))
stopifnot(identical(t(m0), as.matrix(svt2)))
}
\keyword{utilities}
