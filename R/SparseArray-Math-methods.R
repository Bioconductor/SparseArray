### =========================================================================
### 'Math' and 'Math2' methods for SparseArray objects
### -------------------------------------------------------------------------
###
### The 'Math' group consists of the following methods:
### - abs, sign, sqrt, floor, ceiling, trunc
### - cummax, cummin, cumprod, cumsum
### - log, log10, log2, log1p, exp, expm1
### - sin, asin, sinh, asinh, sinpi,
### - cos, acos, cosh, acosh, cospi,
### - tan, atan, tanh, atanh, tanpi,
### - gamma, lgamma, digamma, trigamma
###
### The 'Math2' group consists of the following methods: round, signif
###
### See '?S4groupGeneric' for more information.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Math' group
###

### SparseArray objects only support functions from the 'Math' group that
### propagate zeros.
SUPPORTED_MATH_OPS <- c(
    "abs", "sign", "sqrt", "floor", "ceiling", "trunc",
    "log1p", "expm1",
    "sin", "asin", "sinh", "asinh", "sinpi",
    "tan", "atan", "tanh", "atanh", "tanpi"
)

.check_Math_op <- function(op)
{
    if (!(op %in% SUPPORTED_MATH_OPS))
        stop(wmsg(op, "() is not supported on SparseArray ",
                  "objects (result wouldn't be sparse in general)"))
}

.Math_SVT <- function(op, x)
{
    stopifnot(isSingleString(op), is(x, "SVT_SparseArray"))
    check_svt_version(x)
    .check_Math_op(op)
    if (type(x) != "double")
        stop(wmsg("the ", op, "() method for SVT_SparseArray objects ",
                  "only supports input of type \"double\" at the moment"))
    ans_SVT <- SparseArray.Call("C_Math_SVT",
                                x@dim, x@type, x@SVT, FALSE, op, 0.0)
    new_SVT_SparseArray(x@dim, x@dimnames, "double", ans_SVT, check=FALSE)
}

setMethod("Math", "SVT_SparseArray", function(x) .Math_SVT(.Generic, x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Math2' group
###

.Math2_SVT <- function(op, x, digits)
{
    stopifnot(isSingleString(op), is(x, "SVT_SparseArray"))
    check_svt_version(x)
    if (type(x) != "double")
        stop(wmsg("the ", op, "() method for SVT_SparseArray objects ",
                  "only supports input of type \"double\" at the moment"))
    if (!isSingleNumber(digits))
        stop(wmsg("'digits' must be a single number"))
    if (!is.double(digits))
        digits <- as.double(digits)
    ans_SVT <- SparseArray.Call("C_Math_SVT",
                                x@dim, x@type, x@SVT, FALSE, op, digits)
    new_SVT_SparseArray(x@dim, x@dimnames, "double", ans_SVT, check=FALSE)
}

setMethod("round", "SVT_SparseArray",
    function(x, digits=0) .Math2_SVT("round", x, digits)
)

setMethod("signif", "SVT_SparseArray",
    function(x, digits=6) .Math2_SVT("signif", x, digits)
)

