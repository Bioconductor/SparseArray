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
###   of nonzero elements (i.e. nzcount) in the SVT_SparseArray object.
###   There is no upper limit to this number.
###   In other words, unlike *gCMatrix objects where this number is
###   limited to INT_MAX, an SVT_SparseArray can store an arbitrary number
###   of nonzero elements.
###

setClassUnion("NULL_OR_list", c("NULL", "list"))

SVT_VERSION <- 1L

setClass("SVT_SparseArray",
    contains="SparseArray",
    representation(
        type="character",
        SVT="NULL_OR_list",  # NULL or Sparse Vector Tree (SVT)
        .svt_version="integer"
    ),
    prototype(
        type="logical",
        .svt_version=SVT_VERSION
    )
)

setClass("SVT_SparseMatrix",
    contains=c("SVT_SparseArray", "SparseMatrix"),
    prototype=prototype(
        dim=c(0L, 0L),
        dimnames=list(NULL, NULL)
    )
)

### Not exported (for internal use only).
svt_version <- function(x)
{
    stopifnot(is(x, "SVT_SparseArray") || is(x, "NaArray"))
    if (.hasSlot(x, ".svt_version")) x@.svt_version else 0L
}

check_svt_version <- function(x)
{
    if (svt_version(x) != 0L ||
        is(x, "SVT_SparseArray") && is.null(x@SVT) ||
        is(x, "NaArray") && is.null(x@NaSVT))
    {
        return(invisible(NULL))
    }
    pkg_version <- as.character(packageVersion("SparseArray"))
    stop(wmsg("Old ", class(x)[[1L]], " object detected: object uses ",
              "version 0 of the SVT internal layout which is not ",
              "compatible with versions >= 1.5.0 of the SparseArray ",
              "package (your version is ", pkg_version, ")."))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.validate_SVT_SparseArray <- function(x)
{
    if (!isSingleString(x@type))
        return("'type' slot must be a single string")
    TRUE
}
setValidity2("SVT_SparseArray", .validate_SVT_SparseArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going back and forth between SVT_SparseArray and SVT_SparseMatrix
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
### type() getter and setter
###

setMethod("type", "SVT_SparseArray", function(x) x@type)

.set_SVT_SparseArray_type <- function(x, value)
{
    stopifnot(is(x, "SVT_SparseArray"))
    check_svt_version(x)

    value <- S4Arrays:::normarg_array_type(value, "the supplied type")
    x_type <- type(x)
    if (value == x_type)
        return(x)

    new_SVT <- SparseArray.Call("C_set_SVT_type",
                                x@dim, x@type, x@SVT, FALSE, value)
    BiocGenerics:::replaceSlots(x, type=value, SVT=new_SVT, check=FALSE)
}

setReplaceMethod("type", "SVT_SparseArray", .set_SVT_SparseArray_type)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### is_nonzero(), nzcount(), nzwhich(), nzvals(), `nzvals<-`()
###

### Returns a "logical" SVT_SparseArray object.
.is_nonzero_SVT <- function(x)
{
    stopifnot(is(x, "SVT_SparseArray"))
    check_svt_version(x)
    new_SVT <- SparseArray.Call("C_is_nonzero_SVT", x@dim, x@SVT)
    BiocGenerics:::replaceSlots(x, type="logical", SVT=new_SVT, check=FALSE)
}

setMethod("is_nonzero", "SVT_SparseArray", .is_nonzero_SVT)

### Note that like for the length of atomic vectors in base R, the "nonzero
### count" will be returned as a double if it's > .Machine$integer.max
.nzcount_SVT <- function(x)
{
    stopifnot(is(x, "SVT_SparseArray"))
    check_svt_version(x)
    SparseArray.Call("C_nzcount_SVT", x@dim, x@SVT)
}

setMethod("nzcount", "SVT_SparseArray", .nzcount_SVT)

### Returns an integer vector of length nzcount(x) if 'arr.ind=FALSE', or
### a matrix with nzcount(x) rows if 'arr.ind=TRUE'.
.nzwhich_SVT <- function(x, arr.ind=FALSE)
{
    stopifnot(is(x, "SVT_SparseArray"))
    check_svt_version(x)
    if (!isTRUEorFALSE(arr.ind))
        stop(wmsg("'arr.ind' must be TRUE or FALSE"))
    SparseArray.Call("C_nzwhich_SVT", x@dim, x@SVT, arr.ind)
}

setMethod("nzwhich", "SVT_SparseArray", .nzwhich_SVT)

### TODO: Implement optimized nzvals() and `nzvals<-`() methods for
### SVT_SparseArray objects.


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
### Going back and forth between SVT_SparseArray objects and ordinary arrays
###

.from_SVT_SparseArray_to_array <- function(from)
{
    stopifnot(is(from, "SVT_SparseArray"))
    check_svt_version(from)
    SparseArray.Call("C_from_SVT_SparseArray_to_Rarray",
                     from@dim, dimnames(from), from@type, from@SVT, FALSE)
}

### S3/S4 combo for as.array.SVT_SparseArray
as.array.SVT_SparseArray <- function(x, ...) .from_SVT_SparseArray_to_array(x)
setMethod("as.array", "SVT_SparseArray", as.array.SVT_SparseArray)

.build_SVT_SparseArray_from_array <- function(x, dimnames=NULL, type=NA)
{
    stopifnot(is.array(x))
    if (is.null(dimnames)) {
        ans_dimnames <- dimnames(x)
    } else {
        ans_dimnames <- S4Arrays:::normarg_dimnames(dimnames, dim(x))
    }
    if (identical(type, NA))
        type <- type(x)
    ans_SVT <- SparseArray.Call("C_build_SVT_from_Rarray", x, type, FALSE)
    new_SVT_SparseArray(dim(x), ans_dimnames, type, ans_SVT, check=FALSE)
}

setAs("array", "SVT_SparseArray",
    function(from) .build_SVT_SparseArray_from_array(from)
)
setAs("array", "SparseArray",
    function(from) .build_SVT_SparseArray_from_array(from)
)
setAs("matrix", "SVT_SparseMatrix",
    function(from) .build_SVT_SparseArray_from_array(from)
)
setAs("matrix", "SparseMatrix",
    function(from) .build_SVT_SparseArray_from_array(from)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Make an SVT_SparseMatrix object from CSC components
###

### NOT exported but used in the HDF5Array package!
### 'row_indices' must be an integer vector containing 0-based or 1-based
### row indices. Note that the indices in 'row_indices' are not required
### to be increasing within columns. However, 'row_indices' should NOT
### contain duplicates within columns. This is NOT checked!
make_SVT_SparseMatrix_from_CSC <- function(dim, indptr, data, row_indices,
                                           indices.are.1based=FALSE,
                                           dimnames=NULL)
{
    stopifnot(is.integer(row_indices), isTRUEorFALSE(indices.are.1based))
    ans_SVT <- SparseArray.Call("C_build_SVT_from_CSC",
                                dim, indptr, data, row_indices,
                                indices.are.1based)
    new_SVT_SparseArray(dim, dimnames, type(data), ans_SVT, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going back and forth between SVT_SparseMatrix and CsparseMatrix
###

.make_CsparseMatrix_from_SVT_SparseMatrix <- function(from, to_type=NULL)
{
    stopifnot(is(from, "SVT_SparseMatrix"))
    check_svt_version(from)

    if (!is.null(to_type)) {
        ## Late type switching tends to be slightly more memory efficient.
        ## However, switching to a smaller type (e.g. from "complex" to "double"
        ## or from "integer" to "logical") can introduce zeros. In this case,
        ## we must switch the type early. Otherwise we will end up with zeros
        ## in the "x" slot of the resulting dgCMatrix or lgCMatrix object.
        switch_type_early <- coercion_can_introduce_zeros(from@type, to_type)
        if (switch_type_early)
            type(from) <- to_type  # early type switching
    }

    ## Returns 'ans_p', 'ans_i', and 'ans_x', in a list of length 3.
    C_ans <- SparseArray.Call("C_from_SVT_SparseMatrix_to_CsparseMatrix",
                              from@dim, from@type, from@SVT, is.null(to_type))
    ans_p <- C_ans[[1L]]
    ans_i <- C_ans[[2L]]
    ans_x <- C_ans[[3L]]  # NULL (if 'is.null(to_type)') or same type as 'from'

    if (!is.null(to_type)) {
        ## This type switching is safe only if it does not introduce zeros.
        if (!switch_type_early)
            storage.mode(ans_x) <- to_type  # late type switching
    }

    new_CsparseMatrix(from@dim, ans_p, ans_i, ans_x, dimnames=from@dimnames)
}

.from_SVT_SparseMatrix_to_dgCMatrix <- function(from)
    .make_CsparseMatrix_from_SVT_SparseMatrix(from, "double")
.from_SVT_SparseMatrix_to_lgCMatrix <- function(from)
    .make_CsparseMatrix_from_SVT_SparseMatrix(from, "logical")
.from_SVT_SparseMatrix_to_ngCMatrix <- function(from)
    .make_CsparseMatrix_from_SVT_SparseMatrix(from)

setAs("SVT_SparseMatrix", "dgCMatrix", .from_SVT_SparseMatrix_to_dgCMatrix)
setAs("SVT_SparseMatrix", "lgCMatrix", .from_SVT_SparseMatrix_to_lgCMatrix)
setAs("SVT_SparseMatrix", "ngCMatrix", .from_SVT_SparseMatrix_to_ngCMatrix)

.build_SVT_SparseMatrix_from_CsparseMatrix <- function(x, dimnames=NULL,
                                                          type=NA)
{
    stopifnot(is(x, "CsparseMatrix"))

    ## Turn any [d|l|n]gCMatrix derivative (e.g. TestColMatrix object 'y'
    ## defined in alabaster.matrix/tests/testthat/test-SparseMatrix.R) into
    ## a [d|l|n]gCMatrix **instance**. We should not need to do this. Only
    ## reason we do it is because we don't know how to test for inheritance
    ## at the C level (Rf_inherits() doesn't seem to work properly on S4
    ## objects). More precisely, without this coercion, C function
    ## get_gCMatrix_subtype() defined in src/SVT_SparseArray_class.c won't
    ## be able to recognize a [d|l|n]gCMatrix derivative so will reject it.
    x <- as(x, "CsparseMatrix")

    if (is.null(dimnames)) {
        ans_dimnames <- dimnames(x)
    } else {
        ans_dimnames <- S4Arrays:::normarg_dimnames(dimnames, dim(x))
    }
    if (identical(type, NA))
        type <- type(x)
    ans_SVT <- SparseArray.Call("C_build_SVT_from_CsparseMatrix", x, type)
    new_SVT_SparseArray(dim(x), ans_dimnames, type, ans_SVT, check=FALSE)
}

setAs("CsparseMatrix", "SVT_SparseMatrix",
    function(from) .build_SVT_SparseMatrix_from_CsparseMatrix(from)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going back and forth between SVT_SparseMatrix and TsparseMatrix
###

setAs("SVT_SparseMatrix", "dgTMatrix",
    function(from) as(as(from, "dgCMatrix"), "TsparseMatrix")
)
setAs("SVT_SparseMatrix", "lgTMatrix",
    function(from) as(as(from, "lgCMatrix"), "TsparseMatrix")
)
setAs("SVT_SparseMatrix", "ngTMatrix",
    function(from) as(as(from, "ngCMatrix"), "TsparseMatrix")
)

setAs("TsparseMatrix", "SVT_SparseMatrix",
    function(from) as(as(from, "CsparseMatrix"), "SVT_SparseMatrix")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion from a Matrix derivative to SparseMatrix or SparseArray
###

### Coercing a sparseMatrix derivative (e.g. a CsparseMatrix or TsparseMatrix
### derivative) to SparseMatrix produces an SVT_SparseMatrix object.
### RsparseMatrix is the exception (see COO_SparseArray-class.R).
setAs("sparseMatrix", "SparseMatrix",
    function(from) as(from, "SVT_SparseMatrix")
)

setAs("Matrix", "SparseArray",
    function(from) as(from, "SparseMatrix")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going back and forth between SVT_SparseArray and COO_SparseArray objects
###

.from_SVT_SparseArray_to_COO_SparseArray <- function(from)
{
    stopifnot(is(from, "SVT_SparseArray"))
    check_svt_version(from)
    ## Returns 'ans_nzcoo' and 'ans_nzdata' in a list of length 2.
    C_ans <- SparseArray.Call("C_from_SVT_SparseArray_to_COO_SparseArray",
                              from@dim, from@type, from@SVT)
    ans_nzcoo <- C_ans[[1L]]
    ans_nzdata <- C_ans[[2L]]
    new_COO_SparseArray(from@dim, from@dimnames,
                        ans_nzcoo, ans_nzdata, check=FALSE)
}

setAs("SVT_SparseArray", "COO_SparseArray",
    .from_SVT_SparseArray_to_COO_SparseArray
)
setAs("SVT_SparseMatrix", "COO_SparseMatrix",
    .from_SVT_SparseArray_to_COO_SparseArray
)

.build_SVT_SparseArray_from_COO_SparseArray <- function(x, dimnames=NULL,
                                                           type=NA)
{
    stopifnot(is(x, "COO_SparseArray"))
    if (is.null(dimnames)) {
        ans_dimnames <- dimnames(x)
    } else {
        ans_dimnames <- S4Arrays:::normarg_dimnames(dimnames, dim(x))
    }
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
    ans <- new_SVT_SparseArray(x@dim, ans_dimnames, type, check=FALSE)
    ans[x@nzcoo] <- x@nzdata
    ans
}

setAs("COO_SparseArray", "SVT_SparseArray",
    function(from) .build_SVT_SparseArray_from_COO_SparseArray(from)
)
setAs("COO_SparseMatrix", "SVT_SparseMatrix",
    function(from) .build_SVT_SparseArray_from_COO_SparseArray(from)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Default coercions to SparseArray or SVT_SparseArray
###
### Given a DelayedArray object or any out-of-memory array-like object 'x':
### - as.array(x) is the standard way to realize it as an ordinary array, that
###   is, as a **dense** array);
### - as(x, "SparseArray") is the standard way to realize it in memory as a
###   SparseArray derivative (SVT_SparseArray or COO_SparseArray object), that
###   is as a **sparse** array;
### - as(x, "SVT_SparseArray") is the standard way to realize it in memory
###   as an SVT_SparseArray object.
###

### Similar to as.array() method for Array object (defined in the S4Arrays
### package) but based on extract_sparse_array() instead of extract_array(),
### and without the 'drop' argument.
### Returns an SVT_SparseArray or COO_SparseArray object.
.as_SparseArray <- function(x)
{
    if (!is_sparse(x)) {
        ## Go thru .dense2sparse().
        return(as(x, "COO_SparseArray"))
    }
    index <- vector("list", length=length(dim(x)))
    ans <- extract_sparse_array(x, index)
    S4Arrays:::set_dimnames(ans, dimnames(x))
}
setAs("ANY", "SparseArray", function(from) .as_SparseArray(from))

.as_SparseMatrix <- function(x)
{
    x_ndim <- length(dim(x))
    if (x_ndim != 2L)
        stop(wmsg("cannot coerce ", class(x)[[1L]], " object ",
                  "with ", x_ndim, " dimensions to SparseMatrix ",
                  "(object to coerce must have 2 dimensions)"))
    .as_SparseArray(x)
}
setAs("ANY", "SparseMatrix", function(from) .as_SparseMatrix(from))

setAs("ANY", "SVT_SparseArray",
    function(from) as(as(from, "SparseArray"), "SVT_SparseArray")
)
setAs("ANY", "SVT_SparseMatrix",
    function(from) as(as(from, "SparseMatrix"), "SVT_SparseMatrix")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The SVT_SparseArray() constructor
###

.new_empty_SVT_SparseArray <- function(type=NA)
{
    if (identical(type, NA))
        type <- "logical"
    new2("SVT_SparseArray", type=type, check=FALSE)
}

.SVT_SparseArray <- function(x, dimnames=NULL, type=NA)
{
    if (is.array(x))
        return(.build_SVT_SparseArray_from_array(x,
                                      dimnames=dimnames, type=type))

    if (is(x, "CsparseMatrix"))
        return(.build_SVT_SparseMatrix_from_CsparseMatrix(x,
                                      dimnames=dimnames, type=type))

    if (is(x, "COO_SparseArray"))
        return(.build_SVT_SparseArray_from_COO_SparseArray(x,
                                      dimnames=dimnames, type=type))

    ans <- as(x, "SVT_SparseArray")
    ans <- S4Arrays:::set_dimnames(ans, dimnames)
    if (!identical(type, NA))
        type(ans) <- type
    ans
}

SVT_SparseArray <- function(x, dim=NULL, dimnames=NULL, type=NA)
{
    if (!identical(type, NA))
        type <- S4Arrays:::normarg_array_type(type, "the requested type")

    if (is.null(dim)) {
        if (missing(x))
            return(.new_empty_SVT_SparseArray(type))
        return(.SVT_SparseArray(x, dimnames=dimnames, type=type))
    }

    dim <- S4Arrays:::normarg_dim(dim)
    ans <- new_SVT_SparseArray(dim, dimnames=dimnames, check=FALSE)
    if (!missing(x)) {
        nzidx <- nzwhich(x)
        ans[nzidx] <- as.vector(x[nzidx])
    }
    if (!identical(type, NA))
        type(ans) <- type
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The SparseArray() constructor
###

### Preference is given to the SVT_SparseArray representation.
SparseArray <- function(x, type=NA)
{
    if (!identical(type, NA))
        type <- S4Arrays:::normarg_array_type(type, "the requested type")

    if (missing(x))
        return(.new_empty_SVT_SparseArray(type))

    if (is(x, "SparseArray")) {
        if (is(x, "SVT_SparseArray"))
            check_svt_version(x)
        if (!identical(type, NA))
            type(x) <- type
        return(x)
    }

    ## Calling SparseArray() on a sparseMatrix derivative (e.g. on a
    ## CsparseMatrix or TsparseMatrix object) produces an SVT_SparseMatrix
    ## object. RsparseMatrix is the exception.
    if (is(x, "RsparseMatrix")) {
        ans <- as(x, "COO_SparseMatrix")
        if (!identical(type, NA))
            type(ans) <- type
        return(ans)
    }

    .SVT_SparseArray(x, type=type)
}

