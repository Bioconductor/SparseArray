### =========================================================================
### Ops methods for SparseArray objects
### -------------------------------------------------------------------------
###
### The 'Ops' group of methods consists of three sub groups:
### - 'Arith' group:   "+", "-", "*", "/", "^", "%%", "%/%"
### - 'Compare' group: "==", "!=", "<=", ">=", "<", ">"
### - 'Logic' group:   "&", "|"
###
### See '?S4groupGeneric' for more information.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Arith' group
###

.SVT_SparseArray_Arith <- function(op, e1, e2)
{
    stopifnot(isSingleString(op),
              is(e1, "SVT_SparseArray"),
              is(e2, "SVT_SparseArray"))
    e1_dim <- dim(e1)
    e2_dim <- dim(e2)
    if (!identical(e1_dim, e2_dim))
        stop(wmsg("non-conformable arrays"))

    if (!all(c(type(e1), type(e2)) %in% c("integer", "double", "complex")))
        stop(wmsg("non-numeric argument to arithmetic operation"))

    ## Compute 'ans_dimnames'.
    ans_dimnames <- S4Arrays:::get_first_non_NULL_dimnames(list(e1, e2))

    ## Compute 'ans_type'.
    ans_type <- type(c(vector(type(e1)), vector(type(e2))))
    if (op %in% c("%%", "%/%") && ans_type == "complex")
        stop(wmsg("unimplemented complex operation"))
    if (op %in% c("^", "/") && ans_type == "integer")
        ans_type <- "double"

    ## Pretty inefficient to do this coercion upfront.
    ## TODO: Operate directly on the original types at the C level.
    ## This should be significantly more efficient (but at the cost of some
    ## complication of the C code).
    type(e1) <- type(e2) <- ans_type

    ans_SVT <- .Call2("C_SVT_Arith",
                      e1_dim, e1@type, e1@SVT, e2_dim, e2@type, e2@SVT, op,
                      PACKAGE="SparseArray")

    new_SVT_SparseArray(e1_dim, ans_dimnames, ans_type, ans_SVT, check=FALSE)
}

setMethod("Arith", c("SVT_SparseArray", "SVT_SparseArray"),
    function(e1, e2) .SVT_SparseArray_Arith(.Generic, e1, e2)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Compare' group
###

.SVT_SparseArray_Compare <- function(op, e1, e2)
{
    stopifnot(isSingleString(op),
              is(e1, "SVT_SparseArray"),
              is(e2, "SVT_SparseArray"))
    e1_dim <- dim(e1)
    e2_dim <- dim(e2)
    if (!identical(e1_dim, e2_dim))
        stop(wmsg("non-conformable arrays"))

    ## Compute 'ans_dimnames'.
    ans_dimnames <- S4Arrays:::get_first_non_NULL_dimnames(list(e1, e2))

    ## Pretty inefficient to do this coercion upfront.
    ## TODO: Operate directly on the original types at the C level.
    ## This should be significantly more efficient (but at the cost of some
    ## complication of the C code).
    biggest_type <- type(c(vector(type(e1)), vector(type(e2))))
    type(e1) <- type(e2) <- biggest_type

    ans_SVT <- .Call2("C_SVT_Compare",
                      e1_dim, e1@type, e1@SVT, e2_dim, e2@type, e2@SVT, op,
                      PACKAGE="SparseArray")

    new_SVT_SparseArray(e1_dim, ans_dimnames, "logical", ans_SVT, check=FALSE)
}

setMethod("Compare", c("SVT_SparseArray", "SVT_SparseArray"),
    function(e1, e2) .SVT_SparseArray_Compare(.Generic, e1, e2)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Logic' group
###

.SVT_SparseArray_Logic <- function(op, e1, e2)
{
    stopifnot(isSingleString(op),
              is(e1, "SVT_SparseArray"),
              is(e2, "SVT_SparseArray"))
    e1_dim <- dim(e1)
    e2_dim <- dim(e2)
    if (!identical(e1_dim, e2_dim))
        stop(wmsg("non-conformable arrays"))

    if (type(e1) != "logical" || type(e2) != "logical")
        stop(wmsg("the \"", op, "\" method for SVT_SparseArray objects ",
                  "only supports input objects of type \"logical\" at ",
                  "the moment"))

    ## Compute 'ans_dimnames'.
    ans_dimnames <- S4Arrays:::get_first_non_NULL_dimnames(list(e1, e2))

    ans_SVT <- .Call2("C_SVT_Logic",
                      e1_dim, e1@type, e1@SVT, e2_dim, e2@type, e2@SVT, op,
                      PACKAGE="SparseArray")

    new_SVT_SparseArray(e1_dim, ans_dimnames, "logical", ans_SVT, check=FALSE)
}

setMethod("Logic", c("SVT_SparseArray", "SVT_SparseArray"),
    function(e1, e2) .SVT_SparseArray_Logic(.Generic, e1, e2)
)

