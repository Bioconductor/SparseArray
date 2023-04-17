### =========================================================================
### rowsum() methods for SparseArray and dgCMatrix objects
### -------------------------------------------------------------------------
###


.rowsum_method <- function(x, group, reorder=TRUE, na.rm=FALSE)
{
    stopifnot(is(x, "SVT_SparseMatrix") || is(x, "dgCMatrix"))
    ugroup <- S4Arrays:::compute_ugroup(group, nrow(x), reorder)
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    group <- match(group, ugroup)
    if (is(x, "SVT_SparseMatrix")) {
        ans <- .Call2("C_rowsum_SVT", x@dim, x@type, x@SVT,
                      group, length(ugroup), na.rm,
                      PACKAGE="SparseArray")
    } else {
        ans <- .Call2("C_rowsum_dgCMatrix", x,
                      group, length(ugroup), na.rm,
                      PACKAGE="SparseArray")
    }
    dimnames(ans) <- list(as.character(ugroup), colnames(x))
    ans
}

### S3/S4 combo for rowsum.SVT_SparseMatrix
rowsum.SVT_SparseMatrix <-
    function(x, group, reorder=TRUE, na.rm=FALSE, ...)
        .rowsum_method(x, group, reorder=reorder, na.rm=na.rm, ...)
setMethod("rowsum", "SVT_SparseMatrix",
    function(x, group, reorder=TRUE, ...)
        .rowsum_method(x, group, reorder=reorder, ...)
)

### S3/S4 combo for rowsum.dgCMatrix
rowsum.dgCMatrix <- rowsum.SVT_SparseMatrix
setMethod("rowsum", "dgCMatrix",
    function(x, group, reorder=TRUE, ...)
        .rowsum_method(x, group, reorder=reorder, ...)
)

