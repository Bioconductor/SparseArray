### =========================================================================
### NaArray objects
### -------------------------------------------------------------------------
###
### Like SVT_SparseArray objects but the background value is NA instead of
### zero.
###

setClass("NaArray",
    contains="Array",
    representation(
        dim="integer",
        dimnames="list",
        type="character",
        NaSVT="NULL_OR_list",  # NULL or na-Sparse Vector Tree (NaSVT)
        .svt_version="integer"
    ),
    prototype(
        dim=0L,
        dimnames=list(NULL),
        type="logical",
        .svt_version=SVT_VERSION
    )
)

### Extending RectangularData gives us a few things for free (e.g. validity
### method for RectangularData objects, head(), tail(), etc...). Note
### that even though NaMatrix already extends Array (via NaArray),
### we need to make it a *direct* child of Array, and to list Array *before*
### RectangularData in the 'contains' field below. This will ensure that
### method dispatch will always choose the method for Array in case a generic
### has methods defined for both, Array and RectangularData.
### Note that the fact that we need this "hack" is a hint that we could
### achieve a cleaner class hierarchy by inserting a Matrix class in it.
### Matrix would contain Array and RectangularData (in that order). Then
### NaMatrix would contain NaArray and Matrix (in that order).
### Unfortunately the Matrix package already defines a Matrix class so
### we would need to use a different name.
setClass("NaMatrix",
    contains=c("NaArray", "Array", "RectangularData"),
    prototype=prototype(
        dim=c(0L, 0L),
        dimnames=list(NULL, NULL)
    )
)

.validate_NaMatrix <- function(x)
{
    if (length(x@dim) != 2L)
        return("'dim' slot must be an integer vector of length 2")
    TRUE
}
setValidity2("NaMatrix", .validate_NaMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going back and forth between NaArray and NaMatrix
###

### --- From NaArray to NaMatrix ---

setAs("NaArray", "NaMatrix",
    function(from) new("NaMatrix", from)
)

### --- From NaMatrix to NaArray ---

setAs("NaMatrix", "NaArray", function(from) from)  # no-op

setMethod("coerce", c("NaMatrix", "NaArray"),
    function(from, to, strict=TRUE) from  # no-op
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dim(), dimnames(), and `dimnames<-`()
###

setMethod("dim", "NaArray", function(x) x@dim)

setMethod("dimnames", "NaArray",
    function(x) S4Arrays:::simplify_NULL_dimnames(x@dimnames)
)

setReplaceMethod("dimnames", "NaArray",
    function(x, value)
    {
        x@dimnames <- S4Arrays:::normarg_dimnames(value, dim(x))
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### type() getter and setter
###

setMethod("type", "NaArray", function(x) x@type)

.set_NaArray_type <- function(x, value)
{
    stopifnot(is(x, "NaArray"))
    check_svt_version(x)

    value <- S4Arrays:::normarg_array_type(value, "the supplied type")
    x_type <- type(x)
    if (value == x_type)
        return(x)

    new_NaSVT <- SparseArray.Call("C_set_SVT_SparseArray_type",
                                  x@dim, x@type, x@NaSVT, TRUE, value)
    BiocGenerics:::replaceSlots(x, type=value, NaSVT=new_NaSVT, check=FALSE)
}

setReplaceMethod("type", "NaArray", .set_NaArray_type)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### nnacount() and nnawhich()
###

### Returns the number of non-NA array elements in 'x'.
setGeneric("nnacount", function(x) standardGeneric("nnacount"))

### Note that like for the length of atomic vectors in base R, the "non-NA
### count" will be returned as a double if it's > .Machine$integer.max
.get_NaArray_nnacount <- function(x)
{
    stopifnot(is(x, "NaArray"))
    check_svt_version(x)
    SparseArray.Call("C_nzcount_SVT_SparseArray", x@dim, x@NaSVT)
}

setMethod("nnacount", "NaArray", .get_NaArray_nnacount)

### Returns the indices of the non-NA array elements in 'x', either as
### an L-index (if 'arr.ind=FALSE') or as an M-index (if 'arr.ind=TRUE').
setGeneric("nnawhich", signature="x",
    function(x, arr.ind=FALSE) standardGeneric("nnawhich")
)

### Works on any vector-like or array-like object that supports is.na().
.default_nnawhich <- function(x, arr.ind=FALSE)
{
    if (!isTRUEorFALSE(arr.ind))
        stop(wmsg("'arr.ind' must be TRUE or FALSE"))
    which(!is.na(x), arr.ind=arr.ind, useNames=FALSE)
}
setMethod("nnawhich", "ANY", .default_nnawhich)

### Returns an integer vector of length nnacount(x) if 'arr.ind=FALSE', or
### a matrix with nnacount(x) rows if 'arr.ind=TRUE'.
.nnawhich_NaArray <- function(x, arr.ind=FALSE)
{
    stopifnot(is(x, "NaArray"))
    check_svt_version(x)
    if (!isTRUEorFALSE(arr.ind))
        stop(wmsg("'arr.ind' must be TRUE or FALSE"))
    SparseArray.Call("C_nzwhich_SVT_SparseArray", x@dim, x@NaSVT, arr.ind)
}

setMethod("nnawhich", "NaArray", .nnawhich_NaArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level constructor
###

new_NaArray <- function(dim, dimnames=NULL,
                         type="logical", NaSVT=NULL, check=TRUE)
{
    stopifnot(is.integer(dim))
    if (length(dim) == 2L) {
        ans_class <- "NaMatrix"
    } else {
        ans_class <- "NaArray"
    }
    dimnames <- S4Arrays:::normarg_dimnames(dimnames, dim)
    new2(ans_class, dim=dim, dimnames=dimnames,
                    type=type, NaSVT=NaSVT, check=check)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going back and forth between NaArray objects and ordinary arrays
###

.from_NaArray_to_array <- function(from)
{
    stopifnot(is(from, "NaArray"))
    check_svt_version(from)
    SparseArray.Call("C_from_SVT_SparseArray_to_Rarray",
                     from@dim, dimnames(from), from@type, from@NaSVT, TRUE)
}

### S3/S4 combo for as.array.NaArray
as.array.NaArray <- function(x, ...) .from_NaArray_to_array(x)
setMethod("as.array", "NaArray", as.array.NaArray)

.build_NaArray_from_array <- function(x, dimnames=NULL, type=NA)
{
    stopifnot(is.array(x))
    if (is.null(dimnames)) {
        ans_dimnames <- dimnames(x)
    } else {
        ans_dimnames <- S4Arrays:::normarg_dimnames(dimnames, dim(x))
    }
    if (identical(type, NA))
        type <- type(x)
    ans_NaSVT <- SparseArray.Call("C_build_SVT_from_Rarray", x, type, TRUE)
    new_NaArray(dim(x), ans_dimnames, type, ans_NaSVT, check=FALSE)
}

setAs("array", "NaArray",
    function(from) .build_NaArray_from_array(from)
)
setAs("matrix", "NaMatrix",
    function(from) .build_NaArray_from_array(from)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### NaArray() constructor
###

.new_empty_NaArray <- function(type=NA)
{
    if (identical(type, NA))
        type <- "logical"
    new2("NaArray", type=type, check=FALSE)
}

.NaArray <- function(x, dimnames=NULL, type=NA)
{
    if (is.array(x))
        return(.build_NaArray_from_array(x,
                                      dimnames=dimnames, type=type))

    ans <- as(x, "NaArray")
    ans <- S4Arrays:::set_dimnames(ans, dimnames)
    if (!identical(type, NA))
        type(ans) <- type
    ans
}

NaArray <- function(x, dim=NULL, dimnames=NULL, type=NA)
{
    if (!identical(type, NA))
        type <- S4Arrays:::normarg_array_type(type, "the requested type")

    if (is.null(dim)) {
        if (missing(x))
            return(.new_empty_NaArray(type))
        return(.NaArray(x, dimnames=dimnames, type=type))
    }

    dim <- S4Arrays:::normarg_dim(dim)
    ans <- new_NaArray(dim, dimnames=dimnames, check=FALSE)
    if (!missing(x)) {
        nnaidx <- nnawhich(x)
        ans[nnaidx] <- as.vector(x[nnaidx])
    }
    if (!identical(type, NA))
        type(ans) <- type
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show()
###

.show_nnacount <- function(x)
{
    x_nnacount <- nnacount(x)
    x_density <- x_nnacount / length(x)
    sprintf("[nnacount=%s (%s%%)]", format(x_nnacount),
                                    signif(100 * x_density, digits=2))
}

setMethod("show", "NaArray",
    function(object)
    {
        ## Only reason we print the headline in 2 steps is because we
        ## want to make sure to print at least something (part1) even
        ## when printing part2 is going to fail. This will happen for
        ## example if the call to nnacount() in .show_nnacount() fails.
        cat(show_headline_part1(object))
        cat(.show_nnacount(object))
        if (any(dim(object) == 0L)) {
            cat("\n")
            return()
        }
        cat(":\n", sep="")
        S4Arrays:::print_some_array_elements(object)
    }
)

