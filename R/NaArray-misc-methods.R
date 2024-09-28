### =========================================================================
### Miscellaneous operations on NaArray objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Various "unary isometric" array transformations
###
### A "unary isometric" array transformation is a transformation that returns
### an array-like object with the same dimensions as the input and where each
### element is the result of applying a function to the corresponding element
### in the input.
###
### Note that:
### - Some "unary isometric" transformations preserve sparsity when applied
###   to an NaArray object (e.g. is.nan(), is.finite(), nchar(), etc...) and
###   others don't (e.g. is.na()). NaArray objects only need to support the
###   former.
### - All operations from the 'Math' and 'Math2' groups are "unary isometric"
###   transformations (see '?S4groupGeneric'). The corresponding methods for
###   NaArray objects are implemented in R/NaArray-Math-methods.R
### - All the "unary isometric" methods implemented below return an array-like
###   object of the same class as the input (endomorphism).

### Returns a "logical" **SVT_SparseArray** object!
.isFUN_NaSVT <- function(isFUN, x)
{
    stopifnot(is(x, "NaArray"))
    check_svt_version(x)
    ans_SVT <- SparseArray.Call("C_SVT_apply_isFUN",
                                x@dim, x@type, x@NaSVT, isFUN)
    new_SVT_SparseArray(x@dim, x@dimnames, "logical", ans_SVT, check=FALSE)
}

setMethod("is.na", "NaArray",
    function(x) stop(wmsg("is.na() is not supported on NaArray objects ",
                          "(result wouldn't be sparse in general). Maybe ",
                          "try is_nonna() instead (see '?is_nonna')."))
)
setMethod("is.nan", "NaArray",
    function(x) .isFUN_NaSVT("is.nan", x)
)
setMethod("is.infinite", "NaArray",
    function(x) .isFUN_NaSVT("is.infinite", x)
)

### TODO: Support more methods!


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Various "N-ary isometric" array transformations
###
### An "N-ary isometric" array transformation is a transformation that takes
### one or more array-like objects of the same dimensions (a.k.a. conformable
### arrays) and returns an array-like object of the same dimensions.
###

### COMING SOON...

