### =========================================================================
### NaArray subassignment
### -------------------------------------------------------------------------
###


.subassign_NaSVT_by_Lindex <- function(x, Lindex, value)
{
    x <- adjust_left_type(x, value)
    stopifnot(is.vector(Lindex), is.numeric(Lindex))

    ## No-op (except for type adjustment above) if selection is empty.
    if (length(Lindex) == 0L)
        return(x)

    value <- .normalize_right_value(value, type(x), length(Lindex))

    new_NaSVT <- SparseArray.Call("C_subassign_SVT_by_Lindex",
                                  x@dim, x@type, x@NaSVT, Lindex, value, TRUE)
    BiocGenerics:::replaceSlots(x, NaSVT=new_NaSVT, check=FALSE)
}

setMethod("subassign_Array_by_Lindex", "NaArray",
    function(x, Lindex, value) .subassign_NaSVT_by_Lindex(x, Lindex, value)
)

