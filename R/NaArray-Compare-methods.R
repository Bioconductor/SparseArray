### =========================================================================
### 'Compare' operations on NaArray objects
### -------------------------------------------------------------------------
###
### 'Compare' operations: "==", "!=", "<=", ">=", "<", ">"
###
### See '?S4groupGeneric' for more information.
###


### Supports all 'Compare' ops: "==", "!=", "<=", ">=", "<", ">"
### Returns a "logical" NaArray object.
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

    biggest_type <- type(c(vector(x_type), y))
    type(y) <- biggest_type

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

### Supports all 'Compare' ops: "==", "!=", "<=", ">=", "<", ">"
### Returns a "logical" NaArray object.
.Compare_NaSVT1_NaSVT2 <- function(op, x, y)
{
    stopifnot(isSingleString(op), is(x, "NaArray"), is(y, "NaArray"))
    check_svt_version(x)
    check_svt_version(y)

    ## Check types.
    check_Compare_input_type(type(x), "NaArray object")
    check_Compare_input_type(type(y), "NaArray object")
    check_Compare_op_on_complex_vals(op, type(x), type(y))

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
                                  x_dim, x@type, x@NaSVT, TRUE,
                                  y_dim, y@type, y@NaSVT, TRUE, op)
    new_NaArray(x_dim, ans_dimnames, "logical", ans_NaSVT, check=FALSE)
}

setMethod("Compare", c("NaArray", "NaArray"),
    function(e1, e2) .Compare_NaSVT1_NaSVT2(.Generic, e1, e2)
)

setMethod("Compare", c("NaArray", "array"),
    function(e1, e2) .Compare_NaSVT1_NaSVT2(.Generic, e1, as(e2, "NaArray"))
)

setMethod("Compare", c("array", "NaArray"),
    function(e1, e2) .Compare_NaSVT1_NaSVT2(.Generic, as(e1, "NaArray"), e2)
)

### Supports all 'Compare' ops: "==", "!=", "<=", ">=", "<", ">"
### Returns a "logical" NaArray object.
.Compare_NaSVT1_SVT2 <- function(op, x, y)
{
    stopifnot(isSingleString(op), is(x, "NaArray"), is(y, "SVT_SparseArray"))
    check_svt_version(x)
    check_svt_version(y)

    ## Check types.
    check_Compare_input_type(type(x), "NaArray object")
    check_Compare_input_type(type(y), "SparseArray object")
    check_Compare_op_on_complex_vals(op, type(x), type(y))

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
                                  x_dim, x@type, x@NaSVT, TRUE,
                                  y_dim, y@type, y@SVT, FALSE, op)
    new_NaArray(x_dim, ans_dimnames, "logical", ans_NaSVT, check=FALSE)
}

### Supports all 'Compare' ops: "==", "!=", "<=", ">=", "<", ">"
### Returns a "logical" NaArray object.
.Compare_SVT1_NaSVT2 <- function(op, x, y)
{
    stopifnot(isSingleString(op), is(x, "SVT_SparseArray"), is(y, "NaArray"))
    check_svt_version(x)
    check_svt_version(y)

    ## Check types.
    check_Compare_input_type(type(x), "SparseArray object")
    check_Compare_input_type(type(y), "NaArray object")
    check_Compare_op_on_complex_vals(op, type(x), type(y))

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
                                  x_dim, x@type, x@SVT, FALSE,
                                  y_dim, y@type, y@NaSVT, TRUE, op)
    new_NaArray(x_dim, ans_dimnames, "logical", ans_NaSVT, check=FALSE)
}

setMethod("Compare", c("NaArray", "SVT_SparseArray"),
    function(e1, e2) .Compare_NaSVT1_SVT2(.Generic, e1, e2)
)

setMethod("Compare", c("SVT_SparseArray", "NaArray"),
    function(e1, e2) .Compare_SVT1_NaSVT2(.Generic, e1, e2)
)

