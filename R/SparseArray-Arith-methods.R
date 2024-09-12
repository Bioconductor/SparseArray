### =========================================================================
### 'Arith' operations on SparseArray objects
### -------------------------------------------------------------------------
###
### 'Arith' operations: "+", "-", "*", "/", "^", "%%", "%/%"
###
### See '?S4groupGeneric' for more information.
###
### We also implement unary "+" and "-" for SparseArray objects.
###


ARITH_INPUT_TYPES <- c("integer", "double", "complex")

check_Arith_input_type <- function(type, what)
{
    if (!(type %in% ARITH_INPUT_TYPES))
        stop(wmsg("arithmetic operation not supported ",
                  "on ", what, " of type() \"", type , "\""))
}

op_is_commutative <- function(op)
    (op %in% c("+", "*", "==", "!=", "&", "|"))

error_on_sparsity_not_preserved <- function(op, when)
{
    flipped_op <- flip_Compare_op(op)
    show_flipped_op <- flipped_op != op || op_is_commutative(op)
    if (show_flipped_op) {
        msg <- c("'x ", op, " y' and 'y ", flipped_op, " x': operations")
    } else {
        msg <- c("x ", op, " y: operation")
    }
    stop(wmsg(msg, " not supported on SparseArray object x ",
              "when ", when, " (result wouldn't be sparse)"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Unary "+" and "-"
###

.unary_plus_SparseArray <- function(x)
{
    check_Arith_input_type(type(x), "SparseArray object")
    x  # no-op
}

.unary_minus_SparseArray <- function(x)
{
    check_Arith_input_type(type(x), "SparseArray object")
    if (is(x, "COO_SparseArray")) {
        ans <- BiocGenerics:::replaceSlots(x, nzdata=-x@nzdata, check=FALSE)
    } else if (is(x, "SVT_SparseArray")) {
        check_svt_version(x)
        new_SVT <- SparseArray.Call("C_unary_minus_SVT", x@dim, x@type, x@SVT)
        ans <- BiocGenerics:::replaceSlots(x, SVT=new_SVT, check=FALSE)
    } else {
        stop(wmsg("unary \"-\" is not supported on ", class(x), " objects"))
    }
    ans
}

setMethod("+", c("SparseArray", "missing"),
    function(e1, e2) .unary_plus_SparseArray(e1)
)

setMethod("-", c("SparseArray", "missing"),
    function(e1, e2) .unary_minus_SparseArray(e1)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Arith' group
###

### Supports: "*", "/", "^", "%%", "%/%"
.Arith_SVT1_v2 <- function(op, x, y)
{
    stopifnot(isSingleString(op), is(x, "SVT_SparseArray"))
    check_svt_version(x)

    ## Check types.
    check_Arith_input_type(type(x), "SparseArray object")
    if (!(type(y) %in% ARITH_INPUT_TYPES))
        stop(wmsg("arithmetic operations between SparseArray objects ",
                  "and ", class(y), " vectors are not supported"))

    ## Check 'op'.
    if (!(op %in% c("*", "/", "^", "%%", "%/%")))
        stop(wmsg("\"", op, "\" is not supported between a SparseArray ",
                  "object and a ", class(y), " vector (result wouldn't ",
                  "be sparse in general)"))

    ## Check 'y'.
    if (length(y) != 1L)
        stop(wmsg("arithmetic operations are not supported between a ",
                  "SparseArray object and a vector of length != 1"))
    if (is.na(y))
        error_on_sparsity_not_preserved(op, "y is NA or NaN")
    if (op == "*" && is.infinite(y))
        error_on_sparsity_not_preserved(op, "y is Inf or -Inf")
    if (op == "^" && y <= 0)
        error_on_sparsity_not_preserved(op, "y is non-positive")
    if (op != "*" && y == 0)
        error_on_sparsity_not_preserved(op, "y == 0")

    ## Compute 'ans_type'.
    if (type(y) == "integer" && (type(x) == "double") || op %in% c("/", "^")) {
        ans_type <- type(y) <- "double"
    } else {
        ans_type <- type(c(vector(type(x)), y))
        if (ans_type == "complex")  # temporary
            stop(wmsg("\"", op, "\" is not implemented yet between an ",
                      "SVT_SparseArray object and a single value when ",
                      "one or the other is of type() \"", ans_type, "\""))
        if (op %in% c("%%", "%/%") && ans_type == "complex")
            stop(wmsg("unimplemented complex operation"))
    }

    new_SVT <- SparseArray.Call("C_Arith_SVT1_v2",
                                x@dim, x@type, x@SVT, y, op, ans_type)
    BiocGenerics:::replaceSlots(x, type=ans_type, SVT=new_SVT, check=FALSE)
}

setMethod("Arith", c("SVT_SparseArray", "vector"),
    function(e1, e2) .Arith_SVT1_v2(.Generic, e1, e2)
)

setMethod("Arith", c("vector", "SVT_SparseArray"),
    function(e1, e2) {
        if (.Generic != "*")
            stop(wmsg("\"", .Generic, "\" is not supported between ",
                      "a ", class(e1), " vector on the left and an ",
                      "SVT_SparseArray object on the right (result ",
                      "wouldn't be sparse in general)"))
        .Arith_SVT1_v2(.Generic, e2, e1)
    }
)

### Supports: "+", "-", "*"
.Arith_SVT1_SVT2 <- function(op, x, y)
{
    stopifnot(isSingleString(op),
              is(x, "SVT_SparseArray"),
              is(y, "SVT_SparseArray"))
    check_svt_version(x)
    check_svt_version(y)

    ## Check types.
    check_Arith_input_type(type(x), "SparseArray object")
    check_Arith_input_type(type(y), "SparseArray object")

    ## Check 'op'.
    if (!(op %in% c("+", "-", "*")))
        stop(wmsg("\"", op, "\" is not supported between SparseArray ",
                  "objects (result wouldn't be sparse in general)"))

    ## Check array conformability.
    x_dim <- dim(x)
    y_dim <- dim(y)
    if (!identical(x_dim, y_dim))
        stop(wmsg("non-conformable arrays"))

    ## Compute 'ans_dimnames'.
    ans_dimnames <- S4Arrays:::get_first_non_NULL_dimnames(list(x, y))

    ## Compute 'ans_type'.
    ans_type <- type(c(vector(type(x)), vector(type(y))))
    if (ans_type == "complex")  # temporary
        stop(wmsg("\"", op, "\" is not implemented yet between ",
                  "SVT_SparseArray objects of type() \"", ans_type, "\""))

    ans_SVT <- SparseArray.Call("C_Arith_SVT1_SVT2",
                                x_dim, x@type, x@SVT, y_dim, y@type, y@SVT,
                                op, ans_type)

    new_SVT_SparseArray(x_dim, ans_dimnames, ans_type, ans_SVT, check=FALSE)
}

setMethod("Arith", c("SVT_SparseArray", "SVT_SparseArray"),
    function(e1, e2) .Arith_SVT1_SVT2(.Generic, e1, e2)
)

### We could either:
###  (a) Turn 'e2' into an SVT_SparseArray object and return an
###      SVT_SparseArray object. This is the dgCMatrix approach.
###  (b) Turn 'e1' into an ordinary array and return an ordinary array.
### We choose to do (a). Note that there's no best choice in general as it
### really depends on the 'Arith' operation (i.e. "+", "-", or "*") and
### whether 'e2' has a lot of zeros or not. If it has no or very little
### zeros, then (a) will tend to be less memory efficient than (b) when
### doing "+" or "-". When doing "*", (a) should be always more memory
### efficient than (b), no matter what.
### The cautious user would typically make that choice upfront anyway, by
### coercing one or the other object before calling the 'Arith' op on them.
setMethod("Arith", c("SVT_SparseArray", "array"),
    function(e1, e2) .Arith_SVT1_SVT2(.Generic, e1, as(e2, "SVT_SparseArray"))
)

setMethod("Arith", c("array", "SVT_SparseArray"),
    function(e1, e2) .Arith_SVT1_SVT2(.Generic, as(e1, "SVT_SparseArray"), e2)
)

