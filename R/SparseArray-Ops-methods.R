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

### Supports: "*", "/", "^", "%%", "%/%"
.Arith_SVT1_v2 <- function(op, x, y)
{
    stopifnot(isSingleString(op), is(x, "SVT_SparseArray"))

    ## Check types.
    if (!(type(x) %in% c("integer", "double", "complex")))
        stop(wmsg("arithmetic operations are not suported on ",
                  "SVT_SparseArray objects of type \"", type(x), "\""))
    if (!(is.numeric(y) || is.complex(y)))
        stop(wmsg("arithmetic operations between SVT_SparseArray objects ",
                  "and ", class(y), " vectors are not supported"))

    ## Check 'op'.
    if (!(op %in% c("*", "/", "^", "%%", "%/%")))
        stop(wmsg("\"", op, "\" is not supported between an SVT_SparseArray ",
                  "object and a ", class(y), " vector"))

    ## Check 'y'.
    if (length(y) != 1L)
        stop(wmsg("arithmetic operations are not supported between an ",
                  "SVT_SparseArray object and a ", class(y), " vector ",
                  "of length != 1"))
    if (!is.finite(y))
        stop(wmsg("\"", op, "\" is not supported between an SVT_SparseArray ",
                  "object and a non-finite value (NA, NaN, Inf, or -Inf)"))
    if (op == "^" && y <= 0)
        stop(wmsg("x ", op, " y: operation not supported on ",
                  "SVT_SparseArray object 'x' when 'y' is non-positive"))
    if (op != "*" && y == 0)
        stop(wmsg("x ", op, " 0: operation not supported on ",
                  "SVT_SparseArray object 'x'"))
                
    ## Compute 'ans_type'.
    if (type(y) == "integer" && (type(x) == "double") || op %in% c("/", "^")) {
        ans_type <- type(y) <- "double"
    } else {
        ans_type <- type(c(vector(type(x)), y))
        if (ans_type == "complex")
            stop(wmsg("\"", op, "\" is not implemented yet between an ",
                      "SVT_SparseArray object and a single value when ",
                      "one or the other is of type \"", ans_type, "\""))
        if (op %in% c("%%", "%/%") && ans_type == "complex")
            stop(wmsg("unimplemented complex operation"))
    }

    ans_SVT <- .Call2("C_Arith_SVT1_v2",
                      x@dim, x@type, x@SVT, y, op, ans_type,
                      PACKAGE="SparseArray")

    new_SVT_SparseArray(x@dim, x@dimnames, ans_type, ans_SVT, check=FALSE)
}

setMethod("Arith", c("SVT_SparseArray", "vector"),
    function(e1, e2) .Arith_SVT1_v2(.Generic, e1, e2)
)

setMethod("Arith", c("vector", "SVT_SparseArray"),
    function(e1, e2) {
        if (.Generic != "*")
            stop(wmsg("\"", .Generic, "\" is not supported between ",
                      "a ", class(e1), " vector (on the left) and ",
                      "an SVT_SparseArray object (on the right)"))
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
    if (!(type(x) %in% c("integer", "double", "complex")))
        stop(wmsg("arithmetic operations are not suported on ",
                  "SVT_SparseArray objects of type \"", type(x), "\""))
    if (!(type(y) %in% c("integer", "double", "complex")))
        stop(wmsg("arithmetic operations are not suported on ",
                  "SVT_SparseArray objects of type \"", type(y), "\""))

    ## Check 'op'.
    if (!(op %in% c("+", "-", "*")))
        stop(wmsg("\"", op, "\" is not supported between SVT_SparseArray ",
                  "objects"))

    ## Check array conformability.
    x_dim <- dim(x)
    y_dim <- dim(y)
    if (!identical(x_dim, y_dim))
        stop(wmsg("non-conformable arrays"))

    ## Compute 'ans_dimnames'.
    ans_dimnames <- S4Arrays:::get_first_non_NULL_dimnames(list(x, y))

    ## Compute 'ans_type'.
    ans_type <- type(c(vector(type(x)), vector(type(y))))
    if (ans_type == "complex")
        stop(wmsg("\"", op, "\" is not implemented yet between ",
                  "SVT_SparseArray objects of type \"", ans_type, "\""))

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

.Compare_SVT1_SVT2 <- function(op, x, y)
{
    stopifnot(isSingleString(op),
              is(x, "SVT_SparseArray"),
              is(y, "SVT_SparseArray"))
    x_dim <- dim(x)
    y_dim <- dim(y)
    if (!identical(x_dim, y_dim))
        stop(wmsg("non-conformable arrays"))

    ## Compute 'ans_dimnames'.
    ans_dimnames <- S4Arrays:::get_first_non_NULL_dimnames(list(x, y))

    ## Pretty inefficient to do this coercion upfront.
    ## TODO: Operate directly on the original types at the C level.
    ## This should be significantly more efficient (but at the cost of some
    ## complication of the C code).
    biggest_type <- type(c(vector(type(x)), vector(type(y))))
    type(x) <- type(y) <- biggest_type

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
### 'Logic' group
###

.Logic_SVT1_SVT2 <- function(op, x, y)
{
    stopifnot(isSingleString(op),
              is(x, "SVT_SparseArray"),
              is(y, "SVT_SparseArray"))
    x_dim <- dim(x)
    y_dim <- dim(y)
    if (!identical(x_dim, y_dim))
        stop(wmsg("non-conformable arrays"))

    if (type(x) != "logical" || type(y) != "logical")
        stop(wmsg("the \"", op, "\" method for SVT_SparseArray objects ",
                  "only supports input objects of type \"logical\" at ",
                  "the moment"))

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

