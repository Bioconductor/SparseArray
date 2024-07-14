### =========================================================================
### rowsum() methods for SparseMatrix and dgCMatrix objects
### -------------------------------------------------------------------------
###


.rowsum_method <- function(x, group, reorder=TRUE, na.rm=FALSE)
{
    stopifnot(is(x, "SparseMatrix") || is(x, "dgCMatrix"))
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
        ans <- SparseArray.Call("C_rowsum_dgCMatrix", x,
                                group, length(ugroup), na.rm)
    }
    dimnames(ans) <- list(as.character(ugroup), colnames(x))
    ans
}

### S3/S4 combo for rowsum.SparseMatrix
rowsum.SparseMatrix <-
    function(x, group, reorder=TRUE, na.rm=FALSE, ...)
        .rowsum_method(x, group, reorder=reorder, na.rm=na.rm, ...)
setMethod("rowsum", "SparseMatrix",
    function(x, group, reorder=TRUE, ...)
        .rowsum_method(x, group, reorder=reorder, ...)
)

### S3/S4 combo for rowsum.dgCMatrix
rowsum.dgCMatrix <- rowsum.SparseMatrix
setMethod("rowsum", "dgCMatrix",
    function(x, group, reorder=TRUE, ...)
        .rowsum_method(x, group, reorder=reorder, ...)
)

