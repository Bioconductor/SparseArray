### =========================================================================
### 'Ops' methods for SparseArray objects
### -------------------------------------------------------------------------
###
### The 'Ops' group of methods consists of three sub groups:
### - 'Arith' group:   "+", "-", "*", "/", "^", "%%", "%/%"
### - 'Compare' group: "==", "!=", "<=", ">=", "<", ">"
### - 'Logic' group:   "&", "|"
###
### See '?S4groupGeneric' for more information.
###
### We also implement:
### - unary "+" and "-" methods for SparseArray objects
### - a logical negation ("!") method for SparseArray objects

.ARITH_INPUT_TYPES <- c("integer", "double", "complex")

### All the atomic vector types (i.e. all vector types except "list").
.COMPARE_INPUT_TYPES <- c("logical", "integer", "double", "complex",
                          "character", "raw")

### In base R, "&" and "|" support input of type() "logical", "integer",
### "double", and "complex". We only support "logical" for now.
.LOGIC_INPUT_TYPES <- "logical"

.check_Arith_input_type <- function(type)
{
    if (!(type %in% .ARITH_INPUT_TYPES))
        stop(wmsg("arithmetic operation not supported on SparseArray ",
                  "object of type() \"", type , "\""))
}

.check_Compare_input_type <- function(type)
{
    if (!(type %in% .COMPARE_INPUT_TYPES))
        stop(wmsg("comparison operation not supported on SparseArray ",
                  "object of type() \"", type , "\""))
}

.check_Logic_input_type <- function(type)
{
    if (!(type %in% .LOGIC_INPUT_TYPES))
        stop(wmsg("'Logic' operation \"&\" and \"|\" not supported ",
                  "on SparseArray object of type() \"", type , "\""))
}

.flip_Compare_op <- function(op)
    switch(op, `<=`=">=", `>=`="<=", `<`=">", `>`="<", op)

.op_is_commutative <- function(op)
    (op %in% c("+", "*", "==", "!=", "&", "|"))

.error_on_sparsity_not_preserved <- function(op, when)
{
    flipped_op <- .flip_Compare_op(op)
    show_flipped_op <- flipped_op != op || .op_is_commutative(op)
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
    .check_Arith_input_type(type(x))
    x
}

.unary_minus_SparseArray <- function(x)
{
    .check_Arith_input_type(type(x))
    if (is(x, "COO_SparseArray")) {
        ans <- BiocGenerics:::replaceSlots(x, nzvals=-x@nzvals, check=FALSE)
    } else if (is(x, "SVT_SparseArray")) {
        if (type(x) == "complex")
            stop(wmsg("unary \"-\" is not implemented yet on an ",
                      "SVT_SparseArray object of type \"", type(x), "\""))
        new_SVT <- .Call2("C_unary_minus_SVT", x@dim, x@type, x@SVT,
                          PACKAGE="SparseArray")
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

    ## Check types.
    .check_Arith_input_type(type(x))
    if (!(type(y) %in% .ARITH_INPUT_TYPES))
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
        .error_on_sparsity_not_preserved(op, "y is NA or NaN")
    if (op == "*" && is.infinite(y))
        .error_on_sparsity_not_preserved(op, "y is Inf or -Inf")
    if (op == "^" && y <= 0)
        .error_on_sparsity_not_preserved(op, "y is non-positive")
    if (op != "*" && y == 0)
        .error_on_sparsity_not_preserved(op, "y == 0")

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

    new_SVT <- .Call2("C_Arith_SVT1_v2",
                      x@dim, x@type, x@SVT, y, op, ans_type,
                      PACKAGE="SparseArray")
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

    ## Check types.
    .check_Arith_input_type(type(x))
    .check_Arith_input_type(type(y))

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

    ans_SVT <- .Call2("C_Arith_SVT1_SVT2",
                      x_dim, x@type, x@SVT, y_dim, y@type, y@SVT,
                      op, ans_type,
                      PACKAGE="SparseArray")

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Compare' group
###

.check_Compare_op_on_complex_vals <- function(op, x_type, y_type)
{
    if ((x_type == "complex" || y_type == "complex")
     && op %in% c("<=", ">=", "<", ">"))
        stop(wmsg("invalid comparison with complex values"))
}

.must_homogenize_for_Compare <- function(x_type, y_type)
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

    ## Check types.
    x_type <- type(x)
    .check_Compare_input_type(x_type)
    if (!(type(y) %in% .COMPARE_INPUT_TYPES))
        stop(wmsg("comparison operations between SparseArray objects ",
                  "and ", class(y), " vectors are not supported"))
    .check_Compare_op_on_complex_vals(op, x_type, type(y))

    ## Check 'y'.
    if (length(y) != 1L)
        stop(wmsg("comparison operations are not supported between a ",
                  "SparseArray object and a vector of length != 1"))
    if (is.na(y))
        .error_on_sparsity_not_preserved(op,
                    "y is NA or NaN")
    if (type(y) %in% c("logical", "raw") && op %in% c("<=", "<"))
        .error_on_sparsity_not_preserved(op,
                    "y is a logical or raw value")

    biggest_type <- type(c(vector(x_type), y))
    if (biggest_type == "character" && op %in% c("<=", "<"))
        .error_on_sparsity_not_preserved(op,
                    "type(x) is \"character\" or y is a string")

    type(y) <- biggest_type
    zero <- vector(type(y), length=1L)
    if (op == "==" && y == zero)
        .error_on_sparsity_not_preserved(op,
                    "y is 0 or FALSE or the empty string")
    if (op == "!=" && y != zero)
        .error_on_sparsity_not_preserved(op,
                    "y is not 0, FALSE, or the empty string")
    if (op == "<=" && y >= zero)
        .error_on_sparsity_not_preserved(op,
                    "y is >= 0")
    if (op == ">=" && y <= zero)
        .error_on_sparsity_not_preserved(op,
                    "y is <= 0, or FALSE, or the empty string")
    if (op == "<" && y > zero)
        .error_on_sparsity_not_preserved(op,
                    "y is > 0")
    if (op == ">" && y < zero)
        .error_on_sparsity_not_preserved(op,
                    "y is < 0")

    ## Handle situations where we need to change the type() of 'x' to
    ## the type() of 'y'. This is possibly expensive so we do it only
    ## after all the above checks have passed.
    if (.must_homogenize_for_Compare(type(x), type(y)))
        type(x) <- type(y)

    ## 'type(y)' is guaranteed to be the same as 'type(x)' or a "bigger" type,
    ## considering raw < logical < integer < double < complex < character.
    new_SVT <- .Call2("C_Compare_SVT1_v2",
                      x@dim, x@type, x@SVT, y, op,
                      PACKAGE="SparseArray")
    BiocGenerics:::replaceSlots(x, type="logical", SVT=new_SVT, check=FALSE)
}

setMethod("Compare", c("SVT_SparseArray", "vector"),
    function(e1, e2) .Compare_SVT1_v2(.Generic, e1, e2)
)

setMethod("Compare", c("vector", "SVT_SparseArray"),
    function(e1, e2) .Compare_SVT1_v2(.flip_Compare_op(.Generic), e2, e1)
)

### Supports: "!=", "<", ">"
.Compare_SVT1_SVT2 <- function(op, x, y)
{
    stopifnot(isSingleString(op),
              is(x, "SVT_SparseArray"),
              is(y, "SVT_SparseArray"))

    ## Check types.
    .check_Compare_input_type(type(x))
    .check_Compare_input_type(type(y))
    .check_Compare_op_on_complex_vals(op, type(x), type(y))

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
    if (.must_homogenize_for_Compare(type(x), type(y)))
        type(x) <- type(y) <- type(c(vector(type(x)), vector(type(y))))

    ans_SVT <- .Call2("C_Compare_SVT1_SVT2",
                      x_dim, x@type, x@SVT, y_dim, y@type, y@SVT, op,
                      PACKAGE="SparseArray")

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Logical negation
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

    ## Check types.
    x_type <- type(x)
    .check_Logic_input_type(x_type)
    if (!(type(y) %in% .LOGIC_INPUT_TYPES))
        stop(wmsg("\"", op, "\" between a SparseArray object ",
                  "and a ", class(y), " vector is not supported"))

    ## Check 'y'.
    if (length(y) != 1L)
        stop(wmsg("\"", op, "\" between a SparseArray object ",
                  "and a vector of length != 1 is not supported"))
    if (is.na(y))
        .error_on_sparsity_not_preserved(op, "y is NA")

    if (op == "&" && isFALSE(y))
        return(BiocGenerics:::replaceSlots(x, SVT=NULL, check=FALSE))
    if (op == "|" && isTRUE(y))
        .error_on_sparsity_not_preserved(op, "y is TRUE")
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

    ## Check types.
    .check_Logic_input_type(type(x))
    .check_Logic_input_type(type(y))

    ## Check array conformability.
    x_dim <- dim(x)
    y_dim <- dim(y)
    if (!identical(x_dim, y_dim))
        stop(wmsg("non-conformable arrays"))

    ## Compute 'ans_dimnames'.
    ans_dimnames <- S4Arrays:::get_first_non_NULL_dimnames(list(x, y))

    ans_SVT <- .Call2("C_Logic_SVT1_SVT2",
                      x_dim, x@type, x@SVT, y_dim, y@type, y@SVT, op,
                      PACKAGE="SparseArray")

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

