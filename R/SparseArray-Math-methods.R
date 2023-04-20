### =========================================================================
### Math and Math2 methods for SparseArray objects
### -------------------------------------------------------------------------
###
### The 'Math' group consists of the following methods:
### - abs, sign, sqrt, ceiling, floor, trunc
### - cummax, cummin, cumprod, cumsum
### - log, log10, log2, log1p, acos, acosh
### - asin, asinh, atan, atanh, exp, expm1
### - cos, cosh, cospi, sin, sinh, sinpi, tan, tanh, tanpi
### - gamma, lgamma, digamma, trigamma
###
### The 'Math2' group consists of the following methods: round, signif
###
### See '?S4groupGeneric' for more information.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Math' group
###

.SVT_SparseArray_Math <- function(op, x)
{
    stopifnot(isSingleString(op), is(x, "SVT_SparseArray"))
    if (type(x) != "double")
        stop(wmsg("the ", op, "() method for SVT_SparseArray objects ",
                  "only supports input of type \"double\" at the moment"))

    ans_SVT <- .Call2("C_SVT_Math", x@dim, x@type, x@SVT, op,
                      PACKAGE="SparseArray")

    new_SVT_SparseArray(x@dim, x@dimnames, "double", ans_SVT, check=FALSE)
}

setMethod("Math", "SVT_SparseArray",
    function(x) .SVT_SparseArray_Math(.Generic, x)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Math2' group
###

.SVT_SparseArray_Math2 <- function(op, x, digits)
{
    stopifnot(isSingleString(op), is(x, "SVT_SparseArray"))
    if (type(x) != "double")
        stop(wmsg("the ", op, "() method for SVT_SparseArray objects ",
                  "only supports input of type \"double\""))

    if (!isSingleNumber(digits))
        stop(wmsg("'digits' must be a single number"))
    if (!is.integer(digits))
        digits <- as.integer(digits)

    ans_SVT <- .Call2("C_SVT_Math2", x@dim, x@type, x@SVT, op, digits,
                      PACKAGE="SparseArray")

    new_SVT_SparseArray(x@dim, x@dimnames, "double", ans_SVT, check=FALSE)
}

setMethod("round", c("SVT_SparseArray", "ANY"),
    function(x, digits=0) .SVT_SparseArray_Math2("round", x, digits)
)

setMethod("signif", c("SVT_SparseArray", "ANY"),
    function(x, digits=6) .SVT_SparseArray_Math2("signif", x, digits)
)

