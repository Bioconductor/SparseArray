### =========================================================================
### 'Complex' methods for SparseArray objects
### -------------------------------------------------------------------------
###
### The 'Complex' group consists of the following methods:
###   Re, Im, Mod, Arg, Conj
###
### See '?S4groupGeneric' for more information.
###
### The corresponding base functions all accept an array of type "complex",
### "double", or "integer", and they all return an array of same dimensions
### as the input array.
### Re(), Im(), Mod(), Arg() return an array of type "double" (whatever
### the type of the input array).
### Conj() returns an array of type "complex" if the input array is of type
### "complex', and an array of type "double" otherwise.
### 


.SVT_SparseArray_Complex <- function(op, z)
{
    stopifnot(isSingleString(op), is(z, "SVT_SparseArray"))
    if (type(z) != "complex")
        stop(wmsg("the ", op, "() method for SVT_SparseArray objects ",
                  "only supports input of type \"complex\" at the moment"))

    ## Returns 'ans_type' and 'ans_SVT' in a list of length 2.
    C_ans <- .Call2("C_Complex_SVT", z@dim, z@type, z@SVT, op,
                    PACKAGE="SparseArray")
    ans_type <- C_ans[[1L]]
    ans_SVT <- C_ans[[2L]]

    new_SVT_SparseArray(z@dim, z@dimnames, ans_type, ans_SVT, check=FALSE)
}

setMethod("Complex", "SVT_SparseArray",
    function(z) .SVT_SparseArray_Complex(.Generic, z)
)

