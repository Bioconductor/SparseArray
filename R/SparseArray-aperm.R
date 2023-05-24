### =========================================================================
### Transposition of a SparseArray object
### -------------------------------------------------------------------------
###

.transpose_SVT <- function(x)
{
    stopifnot(is(x, "SVT_SparseArray"))
    new_SVT <- .Call2("C_transpose_SVT",
                      x@dim, x@type, x@SVT, PACKAGE="SparseArray")
    BiocGenerics:::replaceSlots(x, dim=rev(x@dim),
                                   dimnames=rev(x@dimnames),
                                   SVT=new_SVT,
                                   check=FALSE)
}

