### =========================================================================
### 'Math' and 'Math2' methods for SparseArray objects
### -------------------------------------------------------------------------
###
### The 'Math' group consists of the following methods:
### - abs, sign, sqrt, ceiling, floor, trunc
### - cummax, cummin, cumprod, cumsum
### - log, log10, log2, log1p, exp, expm1
### - cos, cospi, acos, sin, sinpi, asin, tan, tanpi, atan
### - cosh, acosh, sinh, asinh, tanh, atanh
### - gamma, lgamma, digamma, trigamma
###
### The 'Math2' group consists of the following methods: round, signif
###
### See '?S4groupGeneric' for more information.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Math' group
###

.Math_SVT <- function(op, x)
{
    stopifnot(isSingleString(op), is(x, "SVT_SparseArray"))
    if (type(x) != "double")
        stop(wmsg("the ", op, "() method for SVT_SparseArray objects ",
                  "only supports input of type \"double\" at the moment"))

    ans_SVT <- .Call2("C_Math_SVT", x@dim, x@type, x@SVT, op,
                      PACKAGE="SparseArray")

    new_SVT_SparseArray(x@dim, x@dimnames, "double", ans_SVT, check=FALSE)
}

setMethod("Math", "SVT_SparseArray", function(x) .Math_SVT(.Generic, x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Math2' group
###

.Math2_SVT <- function(op, x, digits)
{
    stopifnot(isSingleString(op), is(x, "SVT_SparseArray"))
    if (type(x) != "double")
        stop(wmsg("the ", op, "() method for SVT_SparseArray objects ",
                  "only supports input of type \"double\""))

    if (!isSingleNumber(digits))
        stop(wmsg("'digits' must be a single number"))
    if (!is.double(digits))
        digits <- as.double(digits)

    ans_SVT <- .Call2("C_Math2_SVT", x@dim, x@type, x@SVT, op, digits,
                      PACKAGE="SparseArray")

    new_SVT_SparseArray(x@dim, x@dimnames, "double", ans_SVT, check=FALSE)
}

setMethod("round", "SVT_SparseArray",
    function(x, digits=0) .Math2_SVT("round", x, digits)
)

setMethod("signif", "SVT_SparseArray",
    function(x, digits=6) .Math2_SVT("signif", x, digits)
)

