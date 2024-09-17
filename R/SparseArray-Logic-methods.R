### =========================================================================
### 'Logic' operations on SparseArray objects
### -------------------------------------------------------------------------
###
### 'Logic' operations:   "&", "|"
###
### See '?S4groupGeneric' for more information.
###
### We also implement a dummy logical negation ("!") method for SparseArray
### objects.
###


### In base R, "&" and "|" support input of type() "logical", "integer",
### "double", and "complex". We only support "logical" for now.
LOGIC_INPUT_TYPES <- "logical"

check_Logic_input_type <- function(type, what)
{
    if (!(type %in% LOGIC_INPUT_TYPES))
        stop(wmsg("'Logic' operation \"&\" and \"|\" not supported ",
                  "on ", what, " of type() \"", type , "\""))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Dummy logical negation
###

setMethod("!", "SparseArray",
    function(x)
    {
        stop(wmsg("logical negation (\"!\") is not supported on SparseArray ",
                  "objects (result wouldn't be sparse in general)"))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Logic' group
###

.Logic_SVT1_v2 <- function(op, x, y)
{
    stopifnot(isSingleString(op), is(x, "SVT_SparseArray"))
    check_svt_version(x)

    ## Check types.
    x_type <- type(x)
    check_Logic_input_type(x_type, "SparseArray object")
    if (!(type(y) %in% LOGIC_INPUT_TYPES))
        stop(wmsg("\"", op, "\" between a SparseArray object ",
                  "and a ", class(y), " vector is not supported"))

    ## Check 'y'.
    if (length(y) != 1L)
        stop(wmsg("\"", op, "\" between a SparseArray object ",
                  "and a vector of length != 1 is not supported"))
    if (is.na(y))
        error_on_left_sparsity_not_preserved(op, "y is NA")

    if (op == "&" && isFALSE(y))
        return(BiocGenerics:::replaceSlots(x, SVT=NULL, check=FALSE))
    if (op == "|" && isTRUE(y))
        error_on_left_sparsity_not_preserved(op, "y is TRUE")
    x
}

setMethod("Logic", c("SVT_SparseArray", "vector"),
    function(e1, e2) .Logic_SVT1_v2(.Generic, e1, e2)
)

setMethod("Logic", c("vector", "SVT_SparseArray"),
    function(e1, e2) .Logic_SVT1_v2(.Generic, e2, e1)
)

.Logic_SVT1_SVT2 <- function(op, x, y)
{
    stopifnot(isSingleString(op),
              is(x, "SVT_SparseArray"),
              is(y, "SVT_SparseArray"))
    check_svt_version(x)
    check_svt_version(y)

    ## Check types.
    check_Logic_input_type(type(x), "SparseArray object")
    check_Logic_input_type(type(y), "SparseArray object")

    ## Check array conformability.
    x_dim <- dim(x)
    y_dim <- dim(y)
    if (!identical(x_dim, y_dim))
        stop(wmsg("non-conformable arrays"))

    ## Compute 'ans_dimnames'.
    ans_dimnames <- S4Arrays:::get_first_non_NULL_dimnames(list(x, y))

    ans_SVT <- SparseArray.Call("C_Logic_SVT1_SVT2",
                                x_dim, x@type, x@SVT,
                                y_dim, y@type, y@SVT, op)

    new_SVT_SparseArray(x_dim, ans_dimnames, "logical", ans_SVT, check=FALSE)
}

setMethod("Logic", c("SVT_SparseArray", "SVT_SparseArray"),
    function(e1, e2) .Logic_SVT1_SVT2(.Generic, e1, e2)
)

setMethod("Logic", c("SVT_SparseArray", "array"),
    function(e1, e2) .Logic_SVT1_SVT2(.Generic, e1, as(e2, "SVT_SparseArray"))
)

setMethod("Logic", c("array", "SVT_SparseArray"),
    function(e1, e2) .Logic_SVT1_SVT2(.Generic, as(e1, "SVT_SparseArray"), e2)
)

