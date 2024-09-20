### =========================================================================
### 'Logic' operations on NaArray objects
### -------------------------------------------------------------------------
###
### 'Logic' operations:   "&", "|"
###
### See '?S4groupGeneric' for more information.
###
### We also implement a logical negation ("!") method for NaArray objects.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Logical negation
###

.logical_neg_NaArray <- function(x)
{
    stopifnot(is(x, "NaArray"))
    check_svt_version(x)

    ## Check types.
    x_type <- type(x)
    if (!(x_type %in% LOGIC_INPUT_TYPES))
        stop(wmsg("logical negation (\"!\") of an NaArray object ",
                  "of type() \"", x_type, "\" is not supported"))

    new_NaSVT <- SparseArray.Call("C_logical_neg_NaSVT", x@dim, x@type, x@NaSVT)
    BiocGenerics:::replaceSlots(x, type="logical", NaSVT=new_NaSVT, check=FALSE)
}

setMethod("!", "NaArray", .logical_neg_NaArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Logic' group
###

.Logic_NaSVT1_v2 <- function(op, x, y)
{
    stopifnot(isSingleString(op), is(x, "NaArray"))
    check_svt_version(x)

    ## Check types.
    x_type <- type(x)
    check_Logic_input_type(x_type, "NaArray object")
    if (!(type(y) %in% LOGIC_INPUT_TYPES))
        stop(wmsg("\"", op, "\" between an NaArray object ",
                  "and a ", class(y), " vector is not supported"))

    ## Check 'y'.
    if (length(y) != 1L)
        stop(wmsg("\"", op, "\" between an NaArray object ",
                  "and a vector of length != 1 is not supported"))

    if (is.na(y)) {
        new_NaSVT <- SparseArray.Call("C_Logic_NaSVT1_na",
                                      x@dim, x@type, x@NaSVT, op)
        BiocGenerics:::replaceSlots(x, type="logical", NaSVT=new_NaSVT,
                                    check=FALSE)
    } else {
        if (op == "&" && isFALSE(y))
            error_on_left_NAsparsity_not_preserved(op, "y is FALSE")
        if (op == "|" && isTRUE(y))
            error_on_left_NAsparsity_not_preserved(op, "y is TRUE")
        x
    }
}

setMethod("Logic", c("NaArray", "vector"),
    function(e1, e2) .Logic_NaSVT1_v2(.Generic, e1, e2)
)

setMethod("Logic", c("vector", "NaArray"),
    function(e1, e2) .Logic_NaSVT1_v2(.Generic, e2, e1)
)

### Returns an NaArray object.
.Logic_NaSVT1_NaSVT2 <- function(op, x, y)
{
    stopifnot(isSingleString(op), is(x, "NaArray"), is(y, "NaArray"))
    check_svt_version(x)
    check_svt_version(y)

    ## Check types.
    check_Logic_input_type(type(x), "NaArray object")
    check_Logic_input_type(type(y), "NaArray object")

    ## Check array conformability.
    x_dim <- dim(x)
    y_dim <- dim(y)
    if (!identical(x_dim, y_dim))
        stop(wmsg("non-conformable arrays"))

    ## Compute 'ans_dimnames'.
    ans_dimnames <- S4Arrays:::get_first_non_NULL_dimnames(list(x, y))

    ans_NaSVT <- SparseArray.Call("C_Logic_SVT1_SVT2",
                                  x_dim, x@type, x@NaSVT, TRUE,
                                  y_dim, y@type, y@NaSVT, TRUE, op)
    new_NaArray(x_dim, ans_dimnames, NaSVT=ans_NaSVT, check=FALSE)
}

setMethod("Logic", c("NaArray", "NaArray"),
    function(e1, e2) .Logic_NaSVT1_NaSVT2(.Generic, e1, e2)
)

setMethod("Logic", c("NaArray", "array"),
    function(e1, e2) .Logic_NaSVT1_NaSVT2(.Generic, e1, as(e2, "NaArray"))
)

setMethod("Logic", c("array", "NaArray"),
    function(e1, e2) .Logic_NaSVT1_NaSVT2(.Generic, as(e1, "NaArray"), e2)
)

### Returns an NaArray object if 'op' is "|" (because NA | FALSE is NA),
### and a SparseArray if 'op' is "&" (because NA & FALSE is FALSE).
.Logic_NaSVT1_SVT2 <- function(op, x, y)
{
    stopifnot(isSingleString(op), is(x, "NaArray"), is(y, "SVT_SparseArray"))
    check_svt_version(x)
    check_svt_version(y)

    ## Check types.
    check_Logic_input_type(type(x), "NaArray object")
    check_Logic_input_type(type(y), "SparseArray object")

    ## Check array conformability.
    x_dim <- dim(x)
    y_dim <- dim(y)
    if (!identical(x_dim, y_dim))
        stop(wmsg("non-conformable arrays"))

    ## Compute 'ans_dimnames'.
    ans_dimnames <- S4Arrays:::get_first_non_NULL_dimnames(list(x, y))

    ans_SVT <- SparseArray.Call("C_Logic_SVT1_SVT2",
                                x_dim, x@type, x@NaSVT, TRUE,
                                y_dim, y@type, y@SVT, FALSE, op)
    if (op == "|")
        return(new_NaArray(x_dim, ans_dimnames, NaSVT=ans_SVT, check=FALSE))
    new_SVT_SparseArray(x_dim, ans_dimnames, SVT=ans_SVT, check=FALSE)
}

### Returns an NaArray object if 'op' is "|" (because FALSE | NA is NA),
### and a SparseArray if 'op' is "&" (because FALSE & NA is FALSE).
.Logic_SVT1_NaSVT2 <- function(op, x, y)
{
    stopifnot(isSingleString(op), is(x, "SVT_SparseArray"), is(y, "NaArray"))
    check_svt_version(x)
    check_svt_version(y)

    ## Check types.
    check_Logic_input_type(type(x), "SparseArray object")
    check_Logic_input_type(type(y), "NaArray object")

    ## Check array conformability.
    x_dim <- dim(x)
    y_dim <- dim(y)
    if (!identical(x_dim, y_dim))
        stop(wmsg("non-conformable arrays"))

    ## Compute 'ans_dimnames'.
    ans_dimnames <- S4Arrays:::get_first_non_NULL_dimnames(list(x, y))

    ans_SVT <- SparseArray.Call("C_Logic_SVT1_SVT2",
                                x_dim, x@type, x@SVT, FALSE,
                                y_dim, y@type, y@NaSVT, TRUE, op)
    if (op == "|")
        return(new_NaArray(x_dim, ans_dimnames, NaSVT=ans_SVT, check=FALSE))
    new_SVT_SparseArray(x_dim, ans_dimnames, SVT=ans_SVT, check=FALSE)
}

setMethod("Logic", c("NaArray", "SVT_SparseArray"),
    function(e1, e2) .Logic_NaSVT1_SVT2(.Generic, e1, e2)
)

setMethod("Logic", c("SVT_SparseArray", "NaArray"),
    function(e1, e2) .Logic_SVT1_NaSVT2(.Generic, e1, e2)
)

