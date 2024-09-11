### =========================================================================
### 'Ops' methods for NaArray objects
### -------------------------------------------------------------------------
###
### The 'Ops' group of methods consists of three sub groups:
### - 'Arith' group:   "+", "-", "*", "/", "^", "%%", "%/%"
### - 'Compare' group: "==", "!=", "<=", ">=", "<", ">"
### - 'Logic' group:   "&", "|"
###
### See '?S4groupGeneric' for more information.
###
### NaArray objects only support the 'Compare' group for now.

### NOT used at the moment!
.error_on_NaArray_sparsity_not_preserved <- function(op, when)
{
    flipped_op <- flip_Compare_op(op)
    show_flipped_op <- flipped_op != op || op_is_commutative(op)
    if (show_flipped_op) {
        msg <- c("'x ", op, " y' and 'y ", flipped_op, " x': operations")
    } else {
        msg <- c("x ", op, " y: operation")
    }
    stop(wmsg(msg, " not supported on NaArray object x ",
              "when ", when, " (result wouldn't be \"NA-sparse\")"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Compare' group
###

### Supports all 'Compare' ops: "==", "!=", "<=", ">=", "<", ">"
.Compare_NaSVT1_v2 <- function(op, x, y)
{
    stopifnot(isSingleString(op), is(x, "NaArray"))
    check_svt_version(x)

    ## Check types.
    x_type <- type(x)
    check_Compare_input_type(x_type, "NaArray object")
    if (!(type(y) %in% COMPARE_INPUT_TYPES))
        stop(wmsg("comparison operations between NaArray objects ",
                  "and ", class(y), " vectors are not supported"))
    check_Compare_op_on_complex_vals(op, x_type, type(y))

    ## Check 'y'.
    if (length(y) != 1L)
        stop(wmsg("comparison operations are not supported between an ",
                  "NaArray object and a vector of length != 1"))
    #if (is.na(y))
    #    .error_on_NaArray_sparsity_not_preserved(op,
    #                "y is NA or NaN")
    #if (type(y) %in% c("logical", "raw") && op %in% c("<=", "<"))
    #    .error_on_NaArray_sparsity_not_preserved(op,
    #                "y is a logical or raw value")

    biggest_type <- type(c(vector(x_type), y))
    #if (biggest_type == "character" && op %in% c("<=", "<"))
    #    .error_on_NaArray_sparsity_not_preserved(op,
    #                "type(x) is \"character\" or y is a string")

    type(y) <- biggest_type
    #zero <- vector(type(y), length=1L)
    #if (op == "==" && y == zero)
    #    .error_on_NaArray_sparsity_not_preserved(op,
    #                "y is 0 or FALSE or the empty string")
    #if (op == "!=" && y != zero)
    #    .error_on_NaArray_sparsity_not_preserved(op,
    #                "y is not 0, FALSE, or the empty string")
    #if (op == "<=" && y >= zero)
    #    .error_on_NaArray_sparsity_not_preserved(op,
    #                "y is >= 0")
    #if (op == ">=" && y <= zero)
    #    .error_on_NaArray_sparsity_not_preserved(op,
    #                "y is <= 0, or FALSE, or the empty string")
    #if (op == "<" && y > zero)
    #    .error_on_NaArray_sparsity_not_preserved(op,
    #                "y is > 0")
    #if (op == ">" && y < zero)
    #    .error_on_NaArray_sparsity_not_preserved(op,
    #                "y is < 0")

    ## Handle situations where we need to change the type() of 'x' to
    ## the type() of 'y'. This is possibly expensive so we do it only
    ## after all the above checks have passed.
    if (must_homogenize_for_Compare(type(x), type(y)))
        type(x) <- type(y)

    ## 'type(y)' is guaranteed to be the same as 'type(x)' or a "bigger" type,
    ## considering raw < logical < integer < double < complex < character.
    new_NaSVT <- SparseArray.Call("C_Compare_SVT1_v2",
                                  x@dim, x@type, x@NaSVT, TRUE, y, op)
    BiocGenerics:::replaceSlots(x, type="logical", NaSVT=new_NaSVT, check=FALSE)
}

setMethod("Compare", c("NaArray", "vector"),
    function(e1, e2) .Compare_NaSVT1_v2(.Generic, e1, e2)
)

setMethod("Compare", c("vector", "NaArray"),
    function(e1, e2) .Compare_NaSVT1_v2(flip_Compare_op(.Generic), e2, e1)
)

### Supports: "!=", "<", ">"
.Compare_NaSVT1_NaSVT2 <- function(op, x, y)
{
    stopifnot(isSingleString(op), is(x, "NaArray"), is(y, "NaArray"))
    check_svt_version(x)
    check_svt_version(y)

    ## Check types.
    check_Compare_input_type(type(x), "NaArray object")
    check_Compare_input_type(type(y), "NaArray object")
    check_Compare_op_on_complex_vals(op, type(x), type(y))

    ## Check 'op'.
    if (!(op %in% c("!=", "<", ">"))) {
        suggest <- switch(op, `==`="!=", `<=`="<", `>=`=">")
        suggest <- if (is.null(suggest)) "" else
                       paste0(", but \"", suggest, "\" is")
        stop(wmsg("\"", op, "\" is not supported between NaArray ",
                  "objects (result wouldn't be \"NA-sparse\" in general)",
                  suggest))
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

    ans_NaSVT <- SparseArray.Call("C_Compare_SVT1_SVT2",
                                  x_dim, x@type, x@NaSVT,
                                  y_dim, y@type, y@NaSVT, op)

    new_NaArray(x_dim, ans_dimnames, "logical", ans_NaSVT, check=FALSE)
}

#setMethod("Compare", c("NaArray", "NaArray"),
#    function(e1, e2) .Compare_NaSVT1_NaSVT2(.Generic, e1, e2)
#)

#setMethod("Compare", c("NaArray", "array"),
#    function(e1, e2)
#        .Compare_NaSVT1_NaSVT2(.Generic, e1, as(e2, "NaArray"))
#)

#setMethod("Compare", c("array", "NaArray"),
#    function(e1, e2)
#        .Compare_NaSVT1_NaSVT2(.Generic, as(e1, "NaArray"), e2)
#)

