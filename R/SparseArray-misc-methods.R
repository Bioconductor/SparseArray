### =========================================================================
### Miscellaneous operations on SparseArray objects
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
### - Some "unary isometric" transformations preserve sparsity (e.g. is.na(),
###   nchar(), round(), sqrt(), log1p(), etc...) and others don't (e.g.
###   is.finite(), !, log(), etc..). SparseArray objects only need to support
###   the former.
### - All operations from the 'Math' and 'Math2' groups are "unary isometric"
###   transformations (see '?S4groupGeneric') and the corresponding methods
###   for SparseArray objects are implemented in R/SparseArray-Math-methods.R
### - All the "unary isometric" methods implemented below return an array-like
###   object of the same class as the input (endomorphism).

### --- Methods for COO_SparseArray objects ---

.isoFUN_COO <- function(isoFUN, x, ...)
{
    GENERIC <- match.fun(isoFUN)
    new_nzdata <- GENERIC(x@nzdata, ...)
    BiocGenerics:::replaceSlots(x, nzdata=new_nzdata, check=FALSE)
}

setMethod("is.na", "COO_SparseArray",
    function(x) .isoFUN_COO("is.na", x)
)
setMethod("is.nan", "COO_SparseArray",
    function(x) .isoFUN_COO("is.nan", x)
)
setMethod("is.infinite", "COO_SparseArray",
    function(x) .isoFUN_COO("is.infinite", x)
)
setMethod("tolower", "COO_SparseArray",
    function(x) .isoFUN_COO("tolower", x)
)
setMethod("toupper", "COO_SparseArray",
    function(x) .isoFUN_COO("toupper", x)
)
setMethod("nchar", "COO_SparseArray",
    function(x, type="chars", allowNA=FALSE, keepNA=NA)
        .isoFUN_COO("nchar", x, type=type, allowNA=allowNA, keepNA=keepNA)
)

### --- Methods for SVT_SparseArray objects ---

.isFUN_SVT <- function(isFUN, x)
{
    ans_SVT <- SparseArray.Call("C_SVT_apply_isFUN",
                                x@dim, x@type, x@SVT, isFUN)
    new_SVT_SparseArray(x@dim, x@dimnames, "logical", ans_SVT, check=FALSE)
}

setMethod("is.na", "SVT_SparseArray",
    function(x) .isFUN_SVT("is.na", x)
)
setMethod("is.nan", "SVT_SparseArray",
    function(x) .isFUN_SVT("is.nan", x)
)
setMethod("is.infinite", "SVT_SparseArray",
    function(x) .isFUN_SVT("is.infinite", x)
)

### TODO: Support more methods!

