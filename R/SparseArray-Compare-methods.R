### =========================================================================
### 'Compare' operations on SparseArray objects
### -------------------------------------------------------------------------
###
### 'Compare' operations: "==", "!=", "<=", ">=", "<", ">"
###
### See '?S4groupGeneric' for more information.
###


### All the atomic vector types (i.e. all vector types except "list").
COMPARE_INPUT_TYPES <- c("logical", "integer", "double", "complex",
                         "character", "raw")

check_Compare_input_type <- function(type, what)
{
    if (!(type %in% COMPARE_INPUT_TYPES))
        stop(wmsg("comparison operation not supported ",
                  "on ", what, " of type() \"", type , "\""))
}

flip_Compare_op <- function(op)
    switch(op, `<=`=">=", `>=`="<=", `<`=">", `>`="<", op)

check_Compare_op_on_complex_vals <- function(op, x_type, y_type)
{
    if ((x_type == "complex" || y_type == "complex")
     && op %in% c("<=", ">=", "<", ">"))
        stop(wmsg("invalid comparison with complex values"))
}

must_homogenize_for_Compare <- function(x_type, y_type)
{
    if (x_type == "raw" && y_type == "logical" ||
        x_type == "logical" && y_type == "raw")
    {
        ## This is a case where C-level Compare_Rbyte_int() function
        ## (defined in src/leaf_vector_Compare.c) won't compare the
        ## Rbyte values in one object with the int values in the other
        ## object in a meaningful way. That's because the nonzero Rbyte
        ## values can be anything between 1 and 255 while the nonzero
        ## int values are always 1.
        ## An easy workaround is to set the type() of both objects
        ## to "logical".
        return(TRUE)
    }
    if (x_type == "character" || y_type == "character") {
        ## Temporary.
        stop(wmsg("comparison operations are not implemented yet between ",
                  "SVT_SparseArray objects, or between an SVT_SparseArray ",
                  "object and a single value, when one or the other is of ",
                  "type() \"character\""))
        return(TRUE)
    }
    FALSE
}

### Supports all 'Compare' ops: "==", "!=", "<=", ">=", "<", ">"
.Compare_SVT1_v2 <- function(op, x, y)
{
    stopifnot(isSingleString(op), is(x, "SVT_SparseArray"))
    check_svt_version(x)

    ## Check types.
    x_type <- type(x)
    check_Compare_input_type(x_type, "SparseArray object")
    if (!(type(y) %in% COMPARE_INPUT_TYPES))
        stop(wmsg("comparison operations between SparseArray objects ",
                  "and ", class(y), " vectors are not supported"))
    check_Compare_op_on_complex_vals(op, x_type, type(y))

    ## Check 'y'.
    if (length(y) != 1L)
        stop(wmsg("comparison operations are not supported between a ",
                  "SparseArray object and a vector of length != 1"))
    if (is.na(y))
        error_on_sparsity_not_preserved(op,
                    "y is NA or NaN")
    if (type(y) %in% c("logical", "raw") && op %in% c("<=", "<"))
        error_on_sparsity_not_preserved(op,
                    "y is a logical or raw value")

    biggest_type <- type(c(vector(x_type), y))
    if (biggest_type == "character" && op %in% c("<=", "<"))
        error_on_sparsity_not_preserved(op,
                    "type(x) is \"character\" or y is a string")

    type(y) <- biggest_type
    zero <- vector(type(y), length=1L)
    if (op == "==" && y == zero)
        error_on_sparsity_not_preserved(op,
                    "y is 0 or FALSE or the empty string")
    if (op == "!=" && y != zero)
        error_on_sparsity_not_preserved(op,
                    "y is not 0, FALSE, or the empty string")
    if (op == "<=" && y >= zero)
        error_on_sparsity_not_preserved(op,
                    "y is >= 0")
    if (op == ">=" && y <= zero)
        error_on_sparsity_not_preserved(op,
                    "y is <= 0, or FALSE, or the empty string")
    if (op == "<" && y > zero)
        error_on_sparsity_not_preserved(op,
                    "y is > 0")
    if (op == ">" && y < zero)
        error_on_sparsity_not_preserved(op,
                    "y is < 0")

    ## Handle situations where we need to change the type() of 'x' to
    ## the type() of 'y'. This is possibly expensive so we do it only
    ## after all the above checks have passed.
    if (must_homogenize_for_Compare(type(x), type(y)))
        type(x) <- type(y)

    ## 'type(y)' is guaranteed to be the same as 'type(x)' or a "bigger" type,
    ## considering raw < logical < integer < double < complex < character.
    new_SVT <- SparseArray.Call("C_Compare_SVT1_v2",
                                x@dim, x@type, x@SVT, FALSE, y, op)
    BiocGenerics:::replaceSlots(x, type="logical", SVT=new_SVT, check=FALSE)
}

setMethod("Compare", c("SVT_SparseArray", "vector"),
    function(e1, e2) .Compare_SVT1_v2(.Generic, e1, e2)
)

setMethod("Compare", c("vector", "SVT_SparseArray"),
    function(e1, e2) .Compare_SVT1_v2(flip_Compare_op(.Generic), e2, e1)
)

### Supports: "!=", "<", ">"
.Compare_SVT1_SVT2 <- function(op, x, y)
{
    stopifnot(isSingleString(op),
              is(x, "SVT_SparseArray"),
              is(y, "SVT_SparseArray"))
    check_svt_version(x)
    check_svt_version(y)

    ## Check types.
    check_Compare_input_type(type(x), "SparseArray object")
    check_Compare_input_type(type(y), "SparseArray object")
    check_Compare_op_on_complex_vals(op, type(x), type(y))

    ## Check 'op'.
    if (!(op %in% c("!=", "<", ">"))) {
        suggest <- switch(op, `==`="!=", `<=`="<", `>=`=">")
        suggest <- if (is.null(suggest)) "" else
                       paste0(", but \"", suggest, "\" is")
        stop(wmsg("\"", op, "\" is not supported between SparseArray ",
                  "objects (result wouldn't be sparse in general)", suggest))
    }

    ## Check array conformability.
    x_dim <- dim(x)
    y_dim <- dim(y)
    if (!identical(x_dim, y_dim))
        stop(wmsg("non-conformable arrays"))

    ## Compute 'ans_dimnames'.
    ans_dimnames <- S4Arrays:::get_first_non_NULL_dimnames(list(x, y))

    ## Homogenization is possibly expensive so we do it only after all
    ## the above checks have passed.
    if (must_homogenize_for_Compare(type(x), type(y)))
        type(x) <- type(y) <- type(c(vector(type(x)), vector(type(y))))

    ans_SVT <- SparseArray.Call("C_Compare_SVT1_SVT2",
                                x_dim, x@type, x@SVT, y_dim, y@type, y@SVT, op)

    new_SVT_SparseArray(x_dim, ans_dimnames, "logical", ans_SVT, check=FALSE)
}

setMethod("Compare", c("SVT_SparseArray", "SVT_SparseArray"),
    function(e1, e2) .Compare_SVT1_SVT2(.Generic, e1, e2)
)

setMethod("Compare", c("SVT_SparseArray", "array"),
    function(e1, e2) .Compare_SVT1_SVT2(.Generic, e1, as(e2, "SVT_SparseArray"))
)

setMethod("Compare", c("array", "SVT_SparseArray"),
    function(e1, e2) .Compare_SVT1_SVT2(.Generic, as(e1, "SVT_SparseArray"), e2)
)

