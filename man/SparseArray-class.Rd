\name{SparseArray}
\docType{class}

\alias{class:SparseArray}
\alias{SparseArray-class}
\alias{SparseArray}

\alias{class:SparseMatrix}
\alias{SparseMatrix-class}
\alias{SparseMatrix}

\alias{dim,SparseArray-method}
\alias{dimnames,SparseArray-method}
\alias{dimnames<-,SparseArray,ANY-method}
\alias{is_sparse,SparseArray-method}
\alias{nzcount}
\alias{nzcount,CsparseMatrix-method}
\alias{nzcount,RsparseMatrix-method}
\alias{sparsity}

\alias{show,SparseArray-method}

\title{SparseArray objects}

\description{
  The \pkg{SparseArray} package defines the SparseArray virtual class
  whose purpose is to be extended by other S4 classes that aim
  at representing in-memory multidimensional sparse arrays.

  It has currently two concrete subclasses, \link{COO_SparseArray}
  and \link{SVT_SparseArray}, both also defined in this package.
  Each subclass uses its own internal representation for the nonzero
  multidimensional data, the \emph{COO layout} for \link{COO_SparseArray},
  and the \emph{SVT layout} for \link{SVT_SparseArray}. The two layouts
  are described in the \link{COO_SparseArray} and \link{SVT_SparseArray}
  man pages, respectively.

  Finally, the package also defines the SparseMatrix virtual class, as
  a subclass of the SparseArray class, for the specific 2D case.
}

\details{
  The SparseArray class extends the \link[S4Arrays]{Array} virtual class
  defined in the \pkg{S4Arrays} package. Here is the full SparseArray
  sub-hierarchy as defined in the \pkg{SparseArray} package (virtual
  classes are marked with an asterisk):
  \preformatted{
: Array class :                 Array*
: hierarchy   :                   ^
                                  |
- - - - - - - - - - - - - - - - - | - - - - - - - - - - - - - - -
: SparseArray   :            SparseArray*
: sub-hierarchy :            ^    ^    ^
                             |    |    |
                 COO_SparseArray  |  SVT_SparseArray
                        ^         |         ^
- - - - - - - - - - - - | - - - - | - - - - | - - - - - - - - - -
: SparseMatrix      :   |    SparseMatrix*  |
: sub-sub-hierarchy :   |    ^         ^    |
                        |    |         |    |
                 COO_SparseMatrix    SVT_SparseMatrix}

  Any object that belongs to a class that extends SparseArray e.g. (a
  \link{SVT_SparseArray} or \link{SVT_SparseMatrix} object) is called
  a \emph{SparseArray derivative}.

  Most of the \emph{standard array API} defined in base R should work
  on SparseArray derivatives, including \code{dim()}, \code{length()},
  \code{dimnames()}, \code{`dimnames<-`()}, \code{[}, \code{drop()},
  \code{`[<-`} (subassignment), \code{t()}, \code{rbind()}, \code{cbind()},
  etc...

  SparseArray derivatives also support \code{type()}, \code{`type<-`()},
  \code{is_sparse()}, \code{nzcount()}, \code{sparsity()}, \code{arbind()}
  and \code{acbind()}.

  \code{sparsity(x)} returns the ratio between the number of zero-valued
  elements in array-like object \code{x} and its total number of elements
  (\code{length(x)} or \code{prod(dim(x))}).
  More precisely, \code{sparsity(x)} is \code{1 - nzcount(x)/length(x)}.
}

\seealso{
  \itemize{
    \item The \link{COO_SparseArray} and \link{SVT_SparseArray} classes.

    \item \link{SparseArray_subsetting} for subsetting a SparseArray object.

    \item \link{SparseArray_subassignment} for SparseArray subassignment.

    \item \link{SparseArray_combine} for combining SparseArray objects
          by rows or columns.

    \item \link{SparseArray_summarization} for SparseArray summarization
          methods.

    \item \link{SparseArray_misc_ops} for miscellaneous operations on a
          SparseArray object.

    \item \link{SparseMatrix_mult} for SparseMatrix multiplication and
          cross-product.

    \item \link{matrixStats_methods} for SparseMatrix row/col summarization
          methods.

    \item \code{\link{randomSparseArray}} to generate a random SparseArray
          object.

    \item \code{\link{readSparseCSV}} to read/write a sparse matrix
          from/to a CSV (comma-separated values) file.

    \item S4 classes \linkS4class{dgCMatrix} and \linkS4class{lgCMatrix}
          defined in the \pkg{Matrix} package, for the de facto standard
          for sparse matrix representations in the R ecosystem.

    \item \code{\link[S4Arrays]{is_sparse}} in the \pkg{S4Arrays} package.

    \item The \link[S4Arrays]{Array} class defined in the \pkg{S4Arrays}
          package.

    \item Ordinary \link[base]{array} objects in base R.
  }
}

\examples{
showClass("SparseArray")
}

\keyword{methods}
\keyword{classes}
