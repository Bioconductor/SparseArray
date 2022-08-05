### =========================================================================
### matrixStats methods for SparseMatrix objects
### -------------------------------------------------------------------------
###
### matrixStats usage in Bioconductor: Based on some quick grep-based
### inspection, the matrixStats operations used by Bioconductor software
### packages are (looking at the col* functions only):
###   (1) Heavily used: colSums, colMeans, colMedians, colVars, colSds,
###       colMaxs, colMins, colMeans2, colSums2
###   (2) Not so heavily used: colRanges, colRanks, colQuantiles, colMads,
###       colIQRs
###   (3) Marginally used: colAlls, colCumsums, colWeightedMeans, colAnyNAs
###
### Notes:
### - colSums() and colMeans() are base functions (i.e. from the base
###   package), and the colSums() and colMeans() generics are defined
###   in the BiocGenerics package.
### - All other matrix row/column summarization operations are from the
###   matrixStats package with corresponding generics defined in the
###   MatrixGenerics package.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colSums/rowSums
###
### Like base::colSums() and base::rowSums(), the colSums() and rowSums()
### generics in BiocGenerics have a 'dims' argument. We do NOT support it.
### base::colSums() and base::rowSums() propagate the colnames and rownames.

.check_dims <- function(dims, method)
{
    if (!identical(dims, 1))
        stop(wmsg("\"", method, "\" method for SVT_SparseMatrix ",
                  "objects does not support the 'dims' argument"))
}

.colSums_SVT_SparseMatrix <- function(x, na.rm=FALSE, dims=1)
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    .check_dims(dims, "colSums")
    if (is.null(x@SVT)) {
        ans <- numeric(ncol(x))
    } else {
        ans <- vapply(x@SVT,
            function(lv) {
                if (is.null(lv))
                    return(0)
                sum(lv[[2L]], na.rm=na.rm)
            }, numeric(1), USE.NAMES=FALSE)
    }
    setNames(ans, colnames(x))
}
setMethod("colSums", "SVT_SparseMatrix", .colSums_SVT_SparseMatrix)

.rowSums_SVT_SparseMatrix <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "rowSums")
    .colSums_SVT_SparseMatrix(t(x), na.rm=na.rm)
}
setMethod("rowSums", "SVT_SparseMatrix", .rowSums_SVT_SparseMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colCountNAs/rowCountNAs
###
### Not part of the matrixStats API!

.colCountNAs_SVT_SparseMatrix <- function(x)
{
    if (is.null(x@SVT)) {
        ans <- integer(ncol(x))
    } else {
        ans <- vapply(x@SVT,
            function(lv) {
                if (is.null(lv))
                    return(0L)
                sum(is.na(lv[[2L]]))
            }, integer(1), USE.NAMES=FALSE)
    }
    setNames(ans, colnames(x))
}
#setMethod("colCountNAs", "SVT_SparseMatrix", .colCountNAs_SVT_SparseMatrix)

.rowCountNAs_SVT_SparseMatrix <- function(x)
{
    .colCountNAs_SVT_SparseMatrix(t(x))
}
#setMethod("rowCountNAs", "SVT_SparseMatrix", .rowCountNAs_SVT_SparseMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colMeans/rowMeans
###
### Like base::colMeans() and base::rowMeans(), the colMeans() and rowMeans()
### generics in BiocGenerics have a 'dims' argument. We do NOT support it.
### base::colMeans() and base::rowMeans() propagate the colnames and rownames.

.colMeans_SVT_SparseMatrix <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "colMeans")
    sums <- .colSums_SVT_SparseMatrix(x, na.rm=na.rm)
    nvals <- rep.int(nrow(x), ncol(x))
    if (na.rm)
        nvals <- nvals - .colCountNAs_SVT_SparseMatrix(x)
    sums / nvals
}
setMethod("colMeans", "SVT_SparseMatrix", .colMeans_SVT_SparseMatrix)

.rowMeans_SVT_SparseMatrix <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "rowMeans")
    .colMeans_SVT_SparseMatrix(t(x), na.rm=na.rm)
}
setMethod("rowMeans", "SVT_SparseMatrix", .rowMeans_SVT_SparseMatrix)

