### =========================================================================
### rowsum()/colsum() methods for SparseMatrix and dsparseMatrix derivatives
### -------------------------------------------------------------------------
###


.rowsum_method <- function(x, group, reorder=TRUE, na.rm=FALSE)
{
    stopifnot(is(x, "SparseMatrix") || is(x, "dsparseMatrix"))
    ugroup <- S4Arrays:::compute_ugroup(group, nrow(x), reorder)
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    group <- match(group, ugroup)
    if (is(x, "SparseMatrix")) {
        if (is(x, "SVT_SparseMatrix")) {
            check_svt_version(x)
        } else {
            x <- as(x, "SVT_SparseMatrix")
        }
        ans <- SparseArray.Call("C_rowsum_SVT", x@dim, x@type, x@SVT,
                                group, length(ugroup), na.rm)
    } else {
        if (!is(x, "dgCMatrix"))
            x <- as(x, "CsparseMatrix")
        ans <- SparseArray.Call("C_rowsum_dgCMatrix", x,
                                group, length(ugroup), na.rm)
    }
    S4Arrays:::set_dimnames(ans, list(as.character(ugroup), colnames(x)))
}

.colsum_method <- function(x, group, reorder=TRUE, na.rm=FALSE)
{
    stopifnot(is(x, "SparseMatrix") || is(x, "dsparseMatrix"))
    ugroup <- S4Arrays:::compute_ugroup(group, ncol(x), reorder)
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    group <- match(group, ugroup)
    if (is(x, "SparseMatrix")) {
        if (is(x, "SVT_SparseMatrix")) {
            check_svt_version(x)
        } else {
            x <- as(x, "SVT_SparseMatrix")
        }
        ans <- SparseArray.Call("C_colsum_SVT", x@dim, x@type, x@SVT,
                                group, length(ugroup), na.rm)
    } else {
        if (!is(x, "dgCMatrix"))
            x <- as(x, "CsparseMatrix")
        ans <- SparseArray.Call("C_colsum_dgCMatrix", x,
                                group, length(ugroup), na.rm)
    }
    S4Arrays:::set_dimnames(ans, list(rownames(x), as.character(ugroup)))
}

### S3/S4 combo for rowsum.SparseMatrix
rowsum.SparseMatrix <-
    function(x, group, reorder=TRUE, na.rm=FALSE, ...)
        .rowsum_method(x, group, reorder=reorder, na.rm=na.rm, ...)
setMethod("rowsum", "SparseMatrix",
    function(x, group, reorder=TRUE, ...)
        .rowsum_method(x, group, reorder=reorder, ...)
)

### S3/S4 combo for rowsum.dsparseMatrix
rowsum.dsparseMatrix <- rowsum.SparseMatrix
setMethod("rowsum", "dsparseMatrix",
    function(x, group, reorder=TRUE, ...)
        .rowsum_method(x, group, reorder=reorder, ...)
)

setMethod("colsum", "SparseMatrix",
    function(x, group, reorder=TRUE, ...)
        .colsum_method(x, group, reorder=reorder, ...)
)

setMethod("colsum", "dsparseMatrix",
    function(x, group, reorder=TRUE, ...)
        .colsum_method(x, group, reorder=reorder, ...)
)

