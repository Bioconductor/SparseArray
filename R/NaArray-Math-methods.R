### =========================================================================
### 'Math' and 'Math2' methods for NaArray objects
### -------------------------------------------------------------------------
###
### See '?S4groupGeneric' for which functions belong to the 'Math'
### and 'Math2' groups.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Math' group
###

### NaArray objects only support functions from the 'Math' group that
### propagate NAs.
.SUPPORTED_NAARRAY_MATH_OPS <- c(SUPPORTED_MATH_OPS,
    "log", "log10", "log2", "exp",
    "cos", "acos", "cosh", "acosh", "cospi",
    "gamma", "lgamma", "digamma", "trigamma"
)

.check_NaArray_Math_op <- function(op)
{
    if (!(op %in% .SUPPORTED_NAARRAY_MATH_OPS))
        stop(wmsg(op, "() is not supported on NaArray objects ",
                  "(result wouldn't be \"non-NA sparse\" in general)"))
}

.Math_NaSVT <- function(op, x)
{
    stopifnot(isSingleString(op), is(x, "NaArray"))
    check_svt_version(x)
    .check_NaArray_Math_op(op)
    if (type(x) != "double")
        stop(wmsg("the ", op, "() method for NaArray objects ",
                  "only supports input of type \"double\" at the moment"))
    new_NaSVT <- SparseArray.Call("C_Math_SVT",
                                  x@dim, x@type, x@NaSVT, TRUE, op, 0.0)
    BiocGenerics:::replaceSlots(x, type="double", NaSVT=new_NaSVT, check=FALSE)
}

setMethod("Math", "NaArray", function(x) .Math_NaSVT(.Generic, x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Math2' group
###

.Math2_NaSVT <- function(op, x, digits)
{
    stopifnot(isSingleString(op), is(x, "NaArray"))
    check_svt_version(x)
    if (type(x) != "double")
        stop(wmsg("the ", op, "() method for NaArray objects ",
                  "only supports input of type \"double\" at the moment"))
    if (!isSingleNumber(digits))
        stop(wmsg("'digits' must be a single number"))
    if (!is.double(digits))
        digits <- as.double(digits)
    new_NaSVT <- SparseArray.Call("C_Math_SVT",
                                  x@dim, x@type, x@NaSVT, TRUE, op, digits)
    BiocGenerics:::replaceSlots(x, type="double", NaSVT=new_NaSVT, check=FALSE)
}

setMethod("round", "NaArray",
    function(x, digits=0) .Math2_NaSVT("round", x, digits)
)

setMethod("signif", "NaArray",
    function(x, digits=6) .Math2_NaSVT("signif", x, digits)
)

