### =========================================================================
### 'Arith' operations on NaArray objects
### -------------------------------------------------------------------------
###
### 'Arith' operations: "+", "-", "*", "/", "^", "%%", "%/%"
###
### See '?S4groupGeneric' for more information.
###
### We also implement unary "+" and "-" for NaArray objects.
###


error_on_left_NAsparsity_not_preserved <- function(op, when)
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

error_on_right_NAsparsity_not_preserved <- function(op, when)
{
    flipped_op <- flip_Compare_op(op)
    show_flipped_op <- flipped_op != op || op_is_commutative(op)
    if (show_flipped_op) {
        msg <- c("'x ", op, " y' and 'y ", flipped_op, " x': operations")
    } else {
        msg <- c("x ", op, " y: operation")
    }
    stop(wmsg(msg, " not supported on NaArray object y ",
              "when ", when, " (result wouldn't be \"NA-sparse\")"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Unary "+" and "-"
###

.unary_plus_NaArray <- function(x)
{
    stopifnot(is(x, "NaArray"))
    check_Arith_input_type(type(x), "NaArray object")
    x  # no-op
}

.unary_minus_NaArray <- function(x)
{
    stopifnot(is(x, "NaArray"))
    check_Arith_input_type(type(x), "NaArray object")
    check_svt_version(x)
    new_NaSVT <- SparseArray.Call("C_unary_minus_SVT", x@dim, x@type, x@NaSVT)
    BiocGenerics:::replaceSlots(x, NaSVT=new_NaSVT, check=FALSE)
}

setMethod("+", c("NaArray", "missing"),
    function(e1, e2) .unary_plus_NaArray(e1)
)

setMethod("-", c("NaArray", "missing"),
    function(e1, e2) .unary_minus_NaArray(e1)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Arith' group
###

### Supports all 'Arith' ops: "+", "-", "*", "/", "^", "%%", "%/%"
### Returns an NaArray object.
.Arith_NaSVT1_v2 <- function(op, x, y)
{
    stopifnot(isSingleString(op), is(x, "NaArray"))
    check_svt_version(x)

    ## Check types.
    check_Arith_input_type(type(x), "NaArray object")
    if (!(type(y) %in% ARITH_INPUT_TYPES))
        stop(wmsg("arithmetic operations between NaArray objects ",
                  "and ", class(y), " vectors are not supported"))

    ## Check 'y'.
    if (length(y) != 1L)
        stop(wmsg("arithmetic operations are not supported between an ",
                  "NaArray object and a vector of length != 1"))
    if ((op == "^") && (y %in% c(0, NaN)))
        error_on_left_NAsparsity_not_preserved(op, "y is 0 or NaN")
    if ((op == "%%") && !is.na(y) && y == 0)
        error_on_left_NAsparsity_not_preserved(op, "y == 0")

    ## Compute 'ans_type'.
    if (type(x) == "double" && type(y) == "integer" || op %in% c("/", "^"))
        type(y) <- "double"  # cheap
    ans_type <- get_Arith_output_type(op, type(x), type(y))

    new_NaSVT <- SparseArray.Call("C_Arith_SVT1_v2",
                                  x@dim, x@type, x@NaSVT, TRUE, y, op, ans_type)
    BiocGenerics:::replaceSlots(x, type=ans_type, NaSVT=new_NaSVT, check=FALSE)
}

### Supports all 'Arith' ops: "+", "-", "*", "/", "^", "%%", "%/%"
### Returns an NaArray object.
.Arith_v1_NaSVT2 <- function(op, x, y)
{
    stopifnot(isSingleString(op), is(y, "NaArray"))
    check_svt_version(y)

    ## Check types.
    check_Arith_input_type(type(y), "NaArray object")
    if (!(type(x) %in% ARITH_INPUT_TYPES))
        stop(wmsg("arithmetic operations between NaArray objects ",
                  "and ", class(x), " vectors are not supported"))

    ## Check 'x'.
    if (length(x) != 1L)
        stop(wmsg("arithmetic operations are not supported between an ",
                  "NaArray object and a vector of length != 1"))
    if (op == "^" && !is.na(x) && x == 1)
        error_on_right_NAsparsity_not_preserved(op, "x == 1")

    ## Compute 'ans_type'.
    if (type(x) == "integer" && type(y) == "double" || op %in% c("/", "^"))
        type(x) <- "double"  # cheap
    ans_type <- get_Arith_output_type(op, type(x), type(y))

    new_NaSVT <- SparseArray.Call("C_Arith_v1_SVT2",
                                  x, y@dim, y@type, y@NaSVT, TRUE, op, ans_type)
    BiocGenerics:::replaceSlots(y, type=ans_type, NaSVT=new_NaSVT, check=FALSE)
}

setMethod("Arith", c("NaArray", "vector"),
    function(e1, e2) .Arith_NaSVT1_v2(.Generic, e1, e2)
)

setMethod("Arith", c("vector", "NaArray"),
    function(e1, e2) .Arith_v1_NaSVT2(.Generic, e1, e2)
)

### Supports all 'Arith' ops: "+", "-", "*", "/", "^", "%%", "%/%"
### Returns an NaArray object.
.Arith_NaSVT1_NaSVT2 <- function(op, x, y)
{
    stopifnot(isSingleString(op), is(x, "NaArray"), is(y, "NaArray"))
    check_svt_version(x)
    check_svt_version(y)

    ## Check types.
    check_Arith_input_type(type(x), "NaArray object")
    check_Arith_input_type(type(y), "NaArray object")

    ## Check array conformability.
    x_dim <- dim(x)
    y_dim <- dim(y)
    if (!identical(x_dim, y_dim))
        stop(wmsg("non-conformable arrays"))

    ## Compute 'ans_dimnames'.
    ans_dimnames <- S4Arrays:::get_first_non_NULL_dimnames(list(x, y))

    ## Compute 'ans_type'.
    ans_type <- get_Arith_output_type(op, type(x), type(y))

    ans_NaSVT <- SparseArray.Call("C_Arith_SVT1_SVT2",
                                  x_dim, x@type, x@NaSVT, TRUE,
                                  y_dim, y@type, y@NaSVT, TRUE, op, ans_type)
    new_NaArray(x_dim, ans_dimnames, ans_type, ans_NaSVT, check=FALSE)
}

setMethod("Arith", c("NaArray", "NaArray"),
    function(e1, e2) .Arith_NaSVT1_NaSVT2(.Generic, e1, e2)
)

setMethod("Arith", c("NaArray", "array"),
    function(e1, e2) .Arith_NaSVT1_NaSVT2(.Generic, e1, as(e2, "NaArray"))
)

setMethod("Arith", c("array", "NaArray"),
    function(e1, e2) .Arith_NaSVT1_NaSVT2(.Generic, as(e1, "NaArray"), e2)
)

### Supports all 'Arith' ops: "+", "-", "*", "/", "^", "%%", "%/%"
### Returns an NaArray object.
.Arith_NaSVT1_SVT2 <- function(op, x, y)
{
    stopifnot(isSingleString(op), is(x, "NaArray"), is(y, "SVT_SparseArray"))
    check_svt_version(x)
    check_svt_version(y)

    ## Check types.
    check_Arith_input_type(type(x), "NaArray object")
    check_Arith_input_type(type(y), "SparseArray object")

    ## Make sure that result will be NA-sparse.
    if (op == "^")
        stop(wmsg("'x ^ y' is not supported when 'x' is an NaArray object ",
                  "and 'y' a SparseArray object (result wouldn't be ",
                  "\"NA-sparse\" in general)"))
    if (op == "%%" && (type(x) == "double" || type(y) == "double"))
        stop(wmsg("'x %% y' is not supported when 'x' is an NaArray object ",
                  "and 'y' a SparseArray object, and when 'x' or 'y' is of ",
                  "type \"double\" (result wouldn't be \"NA-sparse\" in ",
                  "general)"))

    ## Check array conformability.
    x_dim <- dim(x)
    y_dim <- dim(y)
    if (!identical(x_dim, y_dim))
        stop(wmsg("non-conformable arrays"))

    ## Compute 'ans_dimnames'.
    ans_dimnames <- S4Arrays:::get_first_non_NULL_dimnames(list(x, y))

    ## Compute 'ans_type'.
    ans_type <- get_Arith_output_type(op, type(x), type(y))

    ans_NaSVT <- SparseArray.Call("C_Arith_SVT1_SVT2",
                                  x_dim, x@type, x@NaSVT, TRUE,
                                  y_dim, y@type, y@SVT, FALSE, op, ans_type)
    new_NaArray(x_dim, ans_dimnames, ans_type, ans_NaSVT, check=FALSE)
}

### Supports all 'Arith' ops: "+", "-", "*", "/", "^", "%%", "%/%"
### Returns an NaArray object.
.Arith_SVT1_NaSVT2 <- function(op, x, y)
{
    stopifnot(isSingleString(op), is(x, "SVT_SparseArray"), is(y, "NaArray"))
    check_svt_version(x)
    check_svt_version(y)

    ## Check types.
    check_Arith_input_type(type(x), "SparseArray object")
    check_Arith_input_type(type(y), "NaArray object")

    ## Check array conformability.
    x_dim <- dim(x)
    y_dim <- dim(y)
    if (!identical(x_dim, y_dim))
        stop(wmsg("non-conformable arrays"))

    ## Compute 'ans_dimnames'.
    ans_dimnames <- S4Arrays:::get_first_non_NULL_dimnames(list(x, y))

    ## Compute 'ans_type'.
    ans_type <- get_Arith_output_type(op, type(x), type(y))

    ans_NaSVT <- SparseArray.Call("C_Arith_SVT1_SVT2",
                                  x_dim, x@type, x@SVT, FALSE,
                                  y_dim, y@type, y@NaSVT, TRUE, op, ans_type)
    new_NaArray(x_dim, ans_dimnames, ans_type, ans_NaSVT, check=FALSE)
}

setMethod("Arith", c("NaArray", "SVT_SparseArray"),
    function(e1, e2) .Arith_NaSVT1_SVT2(.Generic, e1, e2)
)

setMethod("Arith", c("SVT_SparseArray", "NaArray"),
    function(e1, e2) .Arith_SVT1_NaSVT2(.Generic, e1, e2)
)

