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
.Arith_SVT_vector <- function(op, e1, e2)
{
    stopifnot(isSingleString(op), is(e1, "SVT_SparseArray"))

    ## Check types.
    if (!(type(e1) %in% c("integer", "double", "complex")))
        stop(wmsg("arithmetic operations are not suported on ",
                  "SVT_SparseArray objects of type \"", type(e1), "\""))
    if (!(is.numeric(e2) || is.complex(e2)))
        stop(wmsg("arithmetic operations between SVT_SparseArray objects ",
                  "and ", class(e2), " vectors are not supported"))

    ## Check 'op'.
    if (!(op %in% c("*", "/", "^", "%%", "%/%")))
        stop(wmsg("\"", op, "\" is not supported between an SVT_SparseArray ",
                  "object and a ", class(e2), " vector"))

    ## Check 'e2'.
    if (length(e2) != 1L)
        stop(wmsg("arithmetic operations are not supported between an ",
                  "SVT_SparseArray object and a ", class(e2), " vector ",
                  "of length != 1"))
    if (!is.finite(e2))
        stop(wmsg("\"", op, "\" is not supported between an SVT_SparseArray ",
                  "object and a non-finite value (NA, NaN, Inf, or -Inf)"))
    if (e2 == 0 && op != "*")
         stop(wmsg("x ", op, " 0: operation not supported on ",
                   "SVT_SparseArray object 'x'"))
                
    ## Compute 'ans_type'.
    ans_type <- type(c(vector(type(e1)), e2))
    if (ans_type == "complex")
        stop(wmsg("\"", op, "\" is not implemented yet between an ",
                  "SVT_SparseArray object and a single value when ",
                  "one or the other is of type \"", ans_type, "\""))
    if (op %in% c("%%", "%/%") && ans_type == "complex")
        stop(wmsg("unimplemented complex operation"))
    if (op %in% c("/", "^") && ans_type == "integer")
        ans_type <- "double"

    ans_SVT <- .Call2("C_Arith_SVT_num",
                      e1@dim, e1@type, e1@SVT, e2, op, ans_type,
                      PACKAGE="SparseArray")

    new_SVT_SparseArray(e1@dim, e1@dimnames, ans_type, ans_SVT, check=FALSE)
}

setMethod("Arith", c("SVT_SparseArray", "vector"),
    function(e1, e2) .Arith_SVT_vector(.Generic, e1, e2)
)

setMethod("Arith", c("vector", "SVT_SparseArray"),
    function(e1, e2) {
        if (.Generic != "*")
            stop(wmsg("\"", .Generic, "\" is not supported between ",
                      "a ", class(e1), " vector (on the left) and ",
                      "an SVT_SparseArray object (on the right)"))
        .Arith_SVT_vector(.Generic, e2, e1)
    }
)

### Supports: "+", "-", "*"
.Arith_SVT_SVT <- function(op, e1, e2)
{
    stopifnot(isSingleString(op),
              is(e1, "SVT_SparseArray"),
              is(e2, "SVT_SparseArray"))

    ## Check types.
    if (!(type(e1) %in% c("integer", "double", "complex")))
        stop(wmsg("arithmetic operations are not suported on ",
                  "SVT_SparseArray objects of type \"", type(e1), "\""))
    if (!(type(e2) %in% c("integer", "double", "complex")))
        stop(wmsg("arithmetic operations are not suported on ",
                  "SVT_SparseArray objects of type \"", type(e1), "\""))

    ## Check 'op'.
    if (!(op %in% c("+", "-", "*")))
        stop(wmsg("\"", op, "\" is not supported between SVT_SparseArray ",
                  "objects"))

    ## Check array conformability.
    e1_dim <- dim(e1)
    e2_dim <- dim(e2)
    if (!identical(e1_dim, e2_dim))
        stop(wmsg("non-conformable arrays"))

    ## Compute 'ans_dimnames'.
    ans_dimnames <- S4Arrays:::get_first_non_NULL_dimnames(list(e1, e2))

    ## Compute 'ans_type'.
    ans_type <- type(c(vector(type(e1)), vector(type(e2))))
    if (ans_type == "complex")
        stop(wmsg("\"", op, "\" is not implemented yet between ",
                  "SVT_SparseArray objects of type \"", ans_type, "\""))

    ans_SVT <- .Call2("C_Arith_SVT_SVT",
                      e1_dim, e1@type, e1@SVT, e2_dim, e2@type, e2@SVT,
                      op, ans_type,
                      PACKAGE="SparseArray")

    new_SVT_SparseArray(e1_dim, ans_dimnames, ans_type, ans_SVT, check=FALSE)
}

setMethod("Arith", c("SVT_SparseArray", "SVT_SparseArray"),
    function(e1, e2) .Arith_SVT_SVT(.Generic, e1, e2)
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
    function(e1, e2) .Arith_SVT_SVT(.Generic, e1, as(e2, "SVT_SparseArray"))
)

setMethod("Arith", c("array", "SVT_SparseArray"),
    function(e1, e2) .Arith_SVT_SVT(.Generic, as(e1, "SVT_SparseArray"), e2)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Compare' group
###

.Compare_SVT_SVT <- function(op, e1, e2)
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

    ans_SVT <- .Call2("C_Compare_SVT_SVT",
                      e1_dim, e1@type, e1@SVT, e2_dim, e2@type, e2@SVT, op,
                      PACKAGE="SparseArray")

    new_SVT_SparseArray(e1_dim, ans_dimnames, "logical", ans_SVT, check=FALSE)
}

setMethod("Compare", c("SVT_SparseArray", "SVT_SparseArray"),
    function(e1, e2) .Compare_SVT_SVT(.Generic, e1, e2)
)

setMethod("Compare", c("SVT_SparseArray", "array"),
    function(e1, e2) .Compare_SVT_SVT(.Generic, e1, as(e2, "SVT_SparseArray"))
)

setMethod("Compare", c("array", "SVT_SparseArray"),
    function(e1, e2) .Compare_SVT_SVT(.Generic, as(e1, "SVT_SparseArray"), e2)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Logic' group
###

.Logic_SVT_SVT <- function(op, e1, e2)
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

    ans_SVT <- .Call2("C_Logic_SVT_SVT",
                      e1_dim, e1@type, e1@SVT, e2_dim, e2@type, e2@SVT, op,
                      PACKAGE="SparseArray")

    new_SVT_SparseArray(e1_dim, ans_dimnames, "logical", ans_SVT, check=FALSE)
}

setMethod("Logic", c("SVT_SparseArray", "SVT_SparseArray"),
    function(e1, e2) .Logic_SVT_SVT(.Generic, e1, e2)
)

setMethod("Logic", c("SVT_SparseArray", "array"),
    function(e1, e2) .Logic_SVT_SVT(.Generic, e1, as(e2, "SVT_SparseArray"))
)

setMethod("Logic", c("array", "SVT_SparseArray"),
    function(e1, e2) .Logic_SVT_SVT(.Generic, as(e1, "SVT_SparseArray"), e2)
)

