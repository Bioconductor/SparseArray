### =========================================================================
### SVT_SparseArray objects
### -------------------------------------------------------------------------
###
### Use SVT layout to store the sparse data.
###
### An SVT_SparseArray object stores its nonzero data in a "Sparse Vector
### Tree" (SVT).
###
### If the sparse array is empty (i.e. has no nonzero data), the 'SVT' slot
### must be set to NULL.
###
### IMPORTANT NOTES:
### - All the "leaf vectors" in the SVT are guaranteed to have a
###   length <= the first dimension of the SVT_SparseArray object, which
###   itself is guaranteed to be <= INT_MAX (2^31 - 1).
### - The cumulated length of the "leaf vectors" in the SVT is the number
###   of nonzero values (i.e. nzdata length) in the SVT_SparseArray object.
###   There is no upper limit to this number.
###   In other words, unlike dgCMatrix objects where this number is
###   limited to INT_MAX, an SVT_SparseArray can store an arbitrary number
###   of nonzero values.
###

setClassUnion("NULL_OR_list", c("NULL", "list"))

setClass("SVT_SparseArray",
    contains="SparseArray",
    representation(
        type="character",
        SVT="NULL_OR_list"  # NULL or Sparse Vector Tree (SVT)
    ),
    prototype(
        type="logical"
    )
)

setClass("SVT_SparseMatrix",
    contains=c("SVT_SparseArray", "SparseMatrix"),
    prototype=prototype(
        dim=c(0L, 0L),
        dimnames=list(NULL, NULL)
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Block going back and forth between SVT_SparseArray and SVT_SparseMatrix
###

### --- From SVT_SparseArray to SVT_SparseMatrix ---

### The user should NOT be able to promote an SVT_SparseArray object to
### SVT_SparseMatrix. Problem is that the automatic coercion method from
### SVT_SparseArray to SVT_SparseMatrix silently returns a broken object
### (unfortunately these dummy automatic coercion methods don't bother to
### validate the object they return). So we overwrite it with a method that
### will fail (as expected) thanks to the validity method for SparseMatrix
### objects.
setAs("SVT_SparseArray", "SVT_SparseMatrix",
    function(from) new("SVT_SparseMatrix", from)
)

### --- From SVT_SparseMatrix to SVT_SparseArray ---

### The user should NOT be able to demote an SVT_SparseMatrix object to
### SVT_SparseArray, so 'as(x, "SVT_SparseArray")' and 'as(x, "SparseArray")'
### should fail or do nothing when 'x' is an SVT_SparseMatrix object, even
### when called with 'strict=TRUE'. Making these coercions behave like no-ops
### seems to be the easiest (and safest) way to go.

setAs("SVT_SparseMatrix", "SVT_SparseArray", function(from) from)  # no-op

### Do NOT use setAs() here! setAs() does really bad things if used to define
### this coercion method e.g. for some reason it calls setIs() internally to
### make SVT_SparseMatrix a **direct** extension of SparseArray, thus
### altering (and breaking) our class hierarchy. This is not only conceptually
### wrong but it also seems to break dispatch e.g. calling 'show(x)' on
### SVT_SparseMatrix object 'x' does not find the method for SparseArray
### objects despite 'is(x, "SparseArray")' being TRUE.
### Worst part is that this seems to be a "feature" (apparently setAs() tries
### to be really smart here!) but it's just a big mess.
setMethod("coerce", c("SVT_SparseMatrix", "SparseArray"),
    function(from, to, strict=TRUE) from  # no-op
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level constructor
###

new_SVT_SparseArray <- function(dim, dimnames=NULL,
                                type="logical", SVT=NULL, check=TRUE)
{
    stopifnot(is.integer(dim))
    if (length(dim) == 2L) {
        ans_class <- "SVT_SparseMatrix"
    } else {
        ans_class <- "SVT_SparseArray"
    }
    dimnames <- S4Arrays:::normarg_dimnames(dimnames, dim)
    new2(ans_class, dim=dim, dimnames=dimnames,
                    type=type, SVT=SVT, check=check)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setMethod("type", "SVT_SparseArray", function(x) x@type)

### Note that like for the length of atomic vectors in base R, the returned
### length will be a double if it's > .Machine$integer.max
.get_SVT_SparseArray_nzcount <- function(x)
{
    stopifnot(is(x, "SVT_SparseArray"))
    .Call2("C_get_SVT_SparseArray_nzcount",
           x@dim, x@SVT, PACKAGE="SparseArray")
}

setMethod("nzcount", "SVT_SparseArray", .get_SVT_SparseArray_nzcount)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### type() setter
###

.set_SVT_SparseArray_type <- function(x, value)
{
    stopifnot(is(x, "SVT_SparseArray"))

    value <- S4Arrays:::normarg_array_type(value, "the supplied type")
    x_type <- type(x)
    if (value == x_type)
        return(x)

    new_SVT <- .Call2("C_set_SVT_SparseArray_type",
                      x@dim, x@type, x@SVT, value,
                      PACKAGE="SparseArray")
    BiocGenerics:::replaceSlots(x, type=value, SVT=new_SVT, check=FALSE)
}

setReplaceMethod("type", "SVT_SparseArray", .set_SVT_SparseArray_type)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going back and forth between SVT_SparseArray objects and ordinary arrays
###

.from_SVT_SparseArray_to_array <- function(from)
{
    stopifnot(is(from, "SVT_SparseArray"))
    .Call2("C_from_SVT_SparseArray_to_Rarray",
           from@dim, dimnames(from), from@type, from@SVT,
           PACKAGE="SparseArray")
}

### S3/S4 combo for as.array.SVT_SparseArray
as.array.SVT_SparseArray <- function(x, ...) .from_SVT_SparseArray_to_array(x)
setMethod("as.array", "SVT_SparseArray", as.array.SVT_SparseArray)

.build_SVT_SparseArray_from_array <- function(x, type=NA)
{
    stopifnot(is.array(x))
    if (identical(type, NA))
        type <- type(x)
    ans_SVT <- .Call2("C_build_SVT_from_Rarray",
                      x, type, PACKAGE="SparseArray")
    new_SVT_SparseArray(dim(x), dimnames(x), type, ans_SVT, check=FALSE)
}

setAs("array", "SVT_SparseArray",
    function(from) .build_SVT_SparseArray_from_array(from)
)
setAs("matrix", "SVT_SparseMatrix",
    function(from) .build_SVT_SparseArray_from_array(from)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going back and forth between SVT_SparseMatrix and [d|l]gCMatrix objects
###

.make_CsparseMatrix_from_SVT_SparseMatrix <- function(from, to_type)
{
    stopifnot(is(from, "SVT_SparseMatrix"))

    ## Late type switching tends to be slightly more memory efficient.
    ## However, switching to a smaller type (e.g. from "complex" to "double"
    ## or from "integer" to "logical") can introduce zeros. In this case,
    ## we must switch the type early. Otherwise we will end up with zeros
    ## in the "x" slot of the resulting dgCMatrix or lgCMatrix object.
    switch_type_early <- coercion_can_introduce_zeros(from@type, to_type)
    if (switch_type_early)
        type(from) <- to_type  # early type switching

    ## Returns 'ans_p', 'ans_i', and 'ans_x', in a list of length 3.
    C_ans <- .Call2("C_from_SVT_SparseMatrix_to_CsparseMatrix",
                    from@dim, from@type, from@SVT, PACKAGE="SparseArray")
    ans_p <- C_ans[[1L]]
    ans_i <- C_ans[[2L]]
    ans_x <- C_ans[[3L]]  # same type as 'from'

    ## This type switching is safe only if it does not introduce zeros.
    if (!switch_type_early)
        storage.mode(ans_x) <- to_type  # late type switching

    new_CsparseMatrix(from@dim, ans_p, ans_i, ans_x, dimnames=from@dimnames)
}

.from_SVT_SparseMatrix_to_dgCMatrix <- function(from)
    .make_CsparseMatrix_from_SVT_SparseMatrix(from, "double")

.from_SVT_SparseMatrix_to_lgCMatrix <- function(from)
    .make_CsparseMatrix_from_SVT_SparseMatrix(from, "logical")

setAs("SVT_SparseMatrix", "dgCMatrix", .from_SVT_SparseMatrix_to_dgCMatrix)
setAs("SVT_SparseMatrix", "lgCMatrix", .from_SVT_SparseMatrix_to_lgCMatrix)

.build_SVT_SparseMatrix_from_CsparseMatrix <- function(x, type=NA)
{
    stopifnot(is(x, "CsparseMatrix"))
    if (identical(type, NA))
        type <- type(x)
    ans_SVT <- .Call2("C_build_SVT_from_CsparseMatrix",
                      x, type, PACKAGE="SparseArray")
    new_SVT_SparseArray(dim(x), dimnames(x), type, ans_SVT, check=FALSE)
}

setAs("CsparseMatrix", "SVT_SparseMatrix",
    function(from) .build_SVT_SparseMatrix_from_CsparseMatrix(from)
)

setAs("Matrix", "SVT_SparseArray", function(from) as(from, "SVT_SparseMatrix"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going back and forth between SVT_SparseArray and COO_SparseArray objects
###

.from_SVT_SparseArray_to_COO_SparseArray <- function(from)
{
    stopifnot(is(from, "SVT_SparseArray"))
    ## Returns 'ans_nzcoo' and 'ans_nzvals' in a list of length 2.
    C_ans <- .Call2("C_from_SVT_SparseArray_to_COO_SparseArray",
                    from@dim, from@type, from@SVT, PACKAGE="SparseArray")
    ans_nzcoo <- C_ans[[1L]]
    ans_nzvals <- C_ans[[2L]]
    new_COO_SparseArray(from@dim, from@dimnames,
                        ans_nzcoo, ans_nzvals, check=FALSE)
}

setAs("SVT_SparseArray", "COO_SparseArray",
    .from_SVT_SparseArray_to_COO_SparseArray
)
setAs("SVT_SparseMatrix", "COO_SparseMatrix",
    .from_SVT_SparseArray_to_COO_SparseArray
)

.build_SVT_SparseArray_from_COO_SparseArray <- function(x, type=NA)
{
    stopifnot(is(x, "COO_SparseArray"))
    if (identical(type, NA)) {
        type <- type(x)
    } else {
        ## Some quick testing/benchmarking seemed to suggest that it's
        ## slightly more efficient to change the type of the input
        ## COO_SparseArray object than that of the output SVT_SparseArray
        ## object.
        type(x) <- type
    }
    ## We start with an allzero SVT_SparseArray object and subassign
    ## the nonzero data to it.
    ans <- new_SVT_SparseArray(x@dim, x@dimnames, type, check=FALSE)
    ans[x@nzcoo] <- x@nzvals
    ans
}

setAs("COO_SparseArray", "SVT_SparseArray",
    function(from) .build_SVT_SparseArray_from_COO_SparseArray(from)
)
setAs("COO_SparseMatrix", "SVT_SparseMatrix",
    function(from) .build_SVT_SparseArray_from_COO_SparseArray(from)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SVT_SparseArray() constructor
###

.new_empty_SVT_SparseArray <- function(type=NA)
{
    if (identical(type, NA))
        type <- "logical"
    new2("SVT_SparseArray", type=type, check=FALSE)
}

.SVT_SparseArray <- function(x, type=NA)
{
    if (is.array(x))
        return(.build_SVT_SparseArray_from_array(x, type=type))

    if (is(x, "CsparseMatrix"))
        return(.build_SVT_SparseMatrix_from_CsparseMatrix(x, type=type))

    if (is(x, "COO_SparseArray"))
        return(.build_SVT_SparseArray_from_COO_SparseArray(x, type=type))

    ans <- as(x, "SVT_SparseArray")
    if (!identical(type, NA))
        type(ans) <- type
    ans
}

SVT_SparseArray <- function(x, type=NA)
{
    if (!identical(type, NA))
        type <- S4Arrays:::normarg_array_type(type, "the requested type")

    if (missing(x))
        return(.new_empty_SVT_SparseArray(type))

    .SVT_SparseArray(x, type=type)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to SparseArray, and the SparseArray() constructor
###
### The SVT_SparseArray representation is preferred over COO_SparseArray when
### coercing to SparseArray or when using the SparseArray() constructor. With
### one exception: when the object to coerce is an RsparseMatrix derivative!
### This is because we don't have an efficient way to coerce these objects to
### SVT_SparseArray at the moment.

setAs("ANY", "SparseArray", function(from) as(from, "SVT_SparseArray"))
setAs("ANY", "SparseMatrix", function(from) as(from, "SVT_SparseMatrix"))
setAs("RsparseMatrix", "SparseArray",
    function(from) as(from, "COO_SparseArray")
)
setAs("RsparseMatrix", "SparseMatrix",
    function(from) as(from, "COO_SparseMatrix")
)

SparseArray <- function(x, type=NA)
{
    if (!identical(type, NA))
        type <- S4Arrays:::normarg_array_type(type, "the requested type")

    if (missing(x))
        return(.new_empty_SVT_SparseArray(type))

    if (is(x, "SparseArray")) {
        if (!identical(type, NA))
            type(x) <- type
        return(x)
    }

    if (is(x, "RsparseMatrix")) {
        ans <- as(x, "COO_SparseArray")
        if (!identical(type, NA))
            type(ans) <- type
        return(ans)
    }

    .SVT_SparseArray(x, type=type)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Transposition
###

### S3/S4 combo for t.SVT_SparseMatrix
t.SVT_SparseMatrix <- function(x)
{
    new_SVT <- .Call2("C_transpose_SVT_SparseMatrix",
                      x@dim, x@type, x@SVT, PACKAGE="SparseArray")
    BiocGenerics:::replaceSlots(x, dim=rev(x@dim),
                                   dimnames=rev(x@dimnames),
                                   SVT=new_SVT,
                                   check=FALSE)
}
setMethod("t", "SVT_SparseMatrix", t.SVT_SparseMatrix)

