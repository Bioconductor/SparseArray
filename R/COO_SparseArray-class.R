### =========================================================================
### COO_SparseArray objects
### -------------------------------------------------------------------------
###
### Use COO layout to store the sparse data.
###
### Same as SparseArraySeed objects in the DelayedArray package.
### Extends the Coordinate List (COO) layout used for sparse matrices to
### multiple dimensions.
### See https://en.wikipedia.org/wiki/Sparse_matrix#Coordinate_list_(COO)
### This layout is also used by https://sparse.pydata.org/
###
### The COO_SparseArray API:
### - The SparseArray API (see SparseArray-class.R)
### - Getters nzcoo() and nzdata()
### - Coercion from array to COO_SparseArray
### - Back and forth coercion between COO_SparseArray and [d|l]g[C|R]Matrix
###   objects from the Matrix package
###

setClass("COO_SparseArray",
    contains="SparseArray",
    representation(
        nzcoo="matrix",  # M-index containing the coordinates of the
                         # nonzero elements.
        nzdata="vector"  # A vector (atomic or list) of length 'nrow(nzcoo)'
                         # containing the nonzero elements.
    ),
    prototype(
        nzcoo=matrix(integer(0), ncol=1L),
        nzdata=logical(0)
    )
)

setClass("COO_SparseMatrix",
    contains=c("COO_SparseArray", "SparseMatrix"),
    prototype=prototype(
        dim=c(0L, 0L),
        dimnames=list(NULL, NULL),
        nzcoo=matrix(integer(0), ncol=2L)
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going back and forth between COO_SparseArray and COO_SparseMatrix
###

### --- From COO_SparseArray to COO_SparseMatrix ---

### The user should NOT be able to promote a COO_SparseArray object to
### COO_SparseMatrix. Problem is that the automatic coercion method from
### COO_SparseArray to COO_SparseMatrix silently returns a broken object
### (unfortunately these dummy automatic coercion methods don't bother to
### validate the object they return). So we overwrite it with a method that
### will fail (as expected) thanks to the validity method for SparseMatrix
### objects.
setAs("COO_SparseArray", "COO_SparseMatrix",
    function(from) new("COO_SparseMatrix", from)
)

### --- From COO_SparseMatrix to COO_SparseArray ---

### The user should NOT be able to demote a COO_SparseMatrix object to
### COO_SparseArray, so 'as(x, "COO_SparseArray")' and 'as(x, "SparseArray")'
### should fail or do nothing when 'x' is a COO_SparseMatrix object, even
### when called with 'strict=TRUE'. Making these coercions behave like no-ops
### seems to be the easiest (and safest) way to go.

setAs("COO_SparseMatrix", "COO_SparseArray", function(from) from)  # no-op

### Do NOT use setAs() here! setAs() does really bad things if used to define
### this coercion method e.g. for some reason it calls setIs() internally to
### make COO_SparseMatrix a **direct** extension of SparseArray, thus
### altering (and breaking) our class hierarchy. This is not only conceptually
### wrong but it also seems to break dispatch e.g. calling 'show(x)' on
### COO_SparseMatrix object 'x' does not find the method for SparseArray
### objects despite 'is(x, "SparseArray")' being TRUE.
### Worst part is that this seems to be a "feature" (apparently setAs() tries
### to be really smart here!) but it's just a big mess.
setMethod("coerce", c("COO_SparseMatrix", "SparseArray"),
    function(from, to, strict=TRUE) from  # no-op
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.validate_nzcoo_slot <- function(x)
{
    x_nzcoo <- x@nzcoo
    if (!(is.matrix(x_nzcoo) && typeof(x_nzcoo) == "integer"))
        return("'nzcoo' slot must be an integer matrix")
    x_dim <- x@dim
    if (ncol(x_nzcoo) != length(x_dim))
        return(paste0("'nzcoo' slot must be a matrix with ",
                      "one column per dimension"))
    for (along in seq_along(x_dim)) {
        not_ok <- S4Vectors:::anyMissingOrOutside(x_nzcoo[ , along],
                                                  1L, x_dim[[along]])
        if (not_ok)
            return(paste0("'nzcoo' slot must contain valid indices, ",
                          "that is, indices that are not NA and are ",
                          ">= 1 and <= their corresponding dimension"))
    }
    TRUE
}

.validate_nzdata_slot <- function(x)
{
    x_nzdata <- x@nzdata
    if (!(is.vector(x_nzdata) && length(x_nzdata) == nrow(x@nzcoo)))
        return(paste0("'nzdata' slot must be a vector of length ",
                      "the number of rows in the 'nzcoo' slot"))
    TRUE
}

.validate_COO_SparseArray <- function(x)
{
    msg <- .validate_nzcoo_slot(x)
    if (!isTRUE(msg))
        return(msg)
    msg <- .validate_nzdata_slot(x)
    if (!isTRUE(msg))
        return(msg)
    TRUE
}
setValidity2("COO_SparseArray", .validate_COO_SparseArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setMethod("type", "COO_SparseArray", function(x) type(x@nzdata))

setGeneric("nzcoo", function(x) standardGeneric("nzcoo"))
setMethod("nzcoo", "COO_SparseArray", function(x) x@nzcoo)

setGeneric("nzdata", function(x) standardGeneric("nzdata"))
setMethod("nzdata", "COO_SparseArray", function(x) x@nzdata)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### type() setter
###

.set_COO_SparseArray_type <- function(x, value)
{
    stopifnot(is(x, "COO_SparseArray"))

    value <- S4Arrays:::normarg_array_type(value, "the supplied type")
    x_type <- type(x)
    if (value == x_type)
        return(x)

    new_nzdata <- x@nzdata
    storage.mode(new_nzdata) <- value
    nzidx <- default_nzwhich(new_nzdata)
    new_nzcoo <- x@nzcoo[nzidx, , drop=FALSE]
    new_nzdata <- new_nzdata[nzidx]
    BiocGenerics:::replaceSlots(x, nzcoo=new_nzcoo,
                                   nzdata=new_nzdata,
                                   check=FALSE)
}

setReplaceMethod("type", "COO_SparseArray", .set_COO_SparseArray_type)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .normalize_COO_SparseArray()
###

.remove_zeros_from_nzdata_slot <- function(x)
{
    zero <- vector(type(x@nzdata), length=1L)
    idx0 <- which(x@nzdata == zero)
    if (length(idx0) == 0L)
        return(x)
    new_nzcoo <- x@nzcoo[-idx0, , drop=FALSE]
    new_nzdata <- x@nzdata[-idx0]
    BiocGenerics:::replaceSlots(x, nzcoo=new_nzcoo,
                                   nzdata=new_nzdata,
                                   check=FALSE)
}

.order_nzcoo_slot <- function(x)
{
    oo <- S4Arrays:::Mindex_order(x@nzcoo)
    if (!is.unsorted(oo))
        return(x)
    new_nzcoo <- x@nzcoo[oo, , drop=FALSE]
    new_nzdata <- x@nzdata[oo]
    BiocGenerics:::replaceSlots(x, nzcoo=new_nzcoo,
                                   nzdata=new_nzdata,
                                   check=FALSE)
}

.normalize_COO_SparseArray <- function(x)
{
    stopifnot(is(x, "COO_SparseArray"))
    .order_nzcoo_slot(.remove_zeros_from_nzdata_slot(x))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### nzcount(), nzwhich(), nzvals(), and `nzvals<-`() methods
###

### length(nzdata(x)) and nrow(nzcoo(x)) are guaranteed to be the same but
### the former should be slightly more efficient.
setMethod("nzcount", "COO_SparseArray", function(x) length(nzdata(x)))


### Returns an integer vector of length nzcount(x) if 'arr.ind=FALSE', or
### a matrix with nzcount(x) rows if 'arr.ind=TRUE'.
.nzwhich_COO_SparseArray <- function(x, arr.ind=FALSE)
{
    if (!isTRUEorFALSE(arr.ind))
        stop(wmsg("'arr.ind' must be TRUE or FALSE"))
    ans <- .normalize_COO_SparseArray(x)@nzcoo
    if (arr.ind)
        return(ans)
    Mindex2Lindex(ans, dim=dim(x))
}

setMethod("nzwhich", "COO_SparseArray", .nzwhich_COO_SparseArray)

setMethod("nzvals", "COO_SparseArray",
    function(x) .normalize_COO_SparseArray(x)@nzdata
)

### As a side effect, the returned COO_SparseArray object is normalized.
setReplaceMethod("nzvals", "COO_SparseArray",
    function(x, value)
    {
        if (!is.vector(value))
            stop(wmsg("replacement value must be a vector"))
        x <- .normalize_COO_SparseArray(x)
        x@nzdata[] <- value
        .remove_zeros_from_nzdata_slot(x)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level constructor
###

new_COO_SparseArray <- function(dim, dimnames=NULL,
                                nzcoo=NULL, nzdata=NULL, check=TRUE)
{
    stopifnot(is.integer(dim))
    if (length(dim) == 2L) {
        ans_class <- "COO_SparseMatrix"
    } else {
        ans_class <- "COO_SparseArray"
    }
    dimnames <- S4Arrays:::normarg_dimnames(dimnames, dim)
    new2(ans_class, dim=dim, dimnames=dimnames,
                    nzcoo=nzcoo, nzdata=nzdata, check=check)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

.normarg_nzdata <- function(nzdata, length.out)
{
    if (is.null(nzdata))
        stop(wmsg("'nzdata' cannot be NULL when 'nzcoo' is not NULL"))
    if (!is.vector(nzdata))
        stop(wmsg("'nzdata' must be a vector"))
    ## Same logic as S4Vectors:::V_recycle().
    nzdata_len <- length(nzdata)
    if (nzdata_len == length.out)
        return(nzdata)
    if (nzdata_len > length.out && nzdata_len != 1L)
        stop(wmsg("'length(nzdata)' is greater than 'nrow(nzcoo)'"))
    if (nzdata_len == 0L)
        stop(wmsg("'length(nzdata)' is 0 but 'nrow(nzcoo)' is not"))
    if (length.out %% nzdata_len != 0L)
        warning(wmsg("'nrow(nzcoo)' is not a multiple of 'length(nzdata)'"))
    rep(nzdata, length.out=length.out)
}

COO_SparseArray <- function(dim, nzcoo=NULL, nzdata=NULL, dimnames=NULL,
                                 check=TRUE)
{
    if (!is.numeric(dim))
        stop(wmsg("'dim' must be an integer vector"))
    if (!is.integer(dim))
        dim <- as.integer(dim)
    if (is.null(nzcoo)) {
        if (is.null(nzdata)) {
            nzdata <- logical(0)  # vector()
        } else if (!(is.vector(nzdata) && length(nzdata) == 0L)) {
            stop(wmsg("'nzdata' must be NULL or a vector of length 0 ",
                      "when 'nzcoo' is NULL"))
        }
        nzcoo <- matrix(integer(0), ncol=length(dim))
    } else {
        if (!is.matrix(nzcoo))
            stop(wmsg("'nzcoo' must be a matrix"))
        if (storage.mode(nzcoo) == "double")
            storage.mode(nzcoo) <- "integer"
        if (!is.null(dimnames(nzcoo)))
            dimnames(nzcoo) <- NULL
        nzdata <- .normarg_nzdata(nzdata, nrow(nzcoo))
    }
    new_COO_SparseArray(dim, dimnames, nzcoo, nzdata, check=check)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .dense2sparse() and .sparse2dense()
###

### Works on any array-like object 'x' that supports nzwhich(x) and nzvals(x).
### Returns a COO_SparseArray object.
.dense2sparse <- function(x)
{
    x_dim <- dim(x)
    if (is.null(x_dim))
        stop(wmsg("'x' must be an array-like object"))
    ans_nzcoo <- nzwhich(x, arr.ind=TRUE)  # M-index
    ans_nzdata <- nzvals(x)
    COO_SparseArray(x_dim, ans_nzcoo, ans_nzdata, dimnames(x), check=FALSE)
}

### 'coo' must be a COO_SparseArray object.
### Return an ordinary array.
.sparse2dense <- function(coo)
{
    if (!is(coo, "COO_SparseArray"))
        stop(wmsg("'coo' must be a COO_SparseArray object"))
    coo_nzdata <- nzdata(coo)
    zero <- vector(typeof(coo_nzdata), length=1L)
    ans <- array(zero, dim=dim(coo))
    ans[nzcoo(coo)] <- coo_nzdata
    S4Arrays:::set_dimnames(ans, dimnames(coo))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion to/from COO_SparseArray
###

### S3/S4 combo for as.array.COO_SparseArray
as.array.COO_SparseArray <- function(x, ...) .sparse2dense(x)
setMethod("as.array", "COO_SparseArray", as.array.COO_SparseArray)

setAs("ANY", "COO_SparseArray", function(from) .dense2sparse(from))
setAs("ANY", "COO_SparseMatrix",
    function(from) as(.dense2sparse(from), "COO_SparseMatrix")
)

### Going back and forth between COO_SparseMatrix and [d|l]g[C|R]Matrix objects
### from the Matrix package:

.make_sparseMatrix_from_COO_SparseMatrix <- function(from, to_type,
                                                     orientation)
{
    stopifnot(is(from, "COO_SparseMatrix"))

    ## Late type switching tends to be slightly more memory efficient.
    ## However, switching to a smaller type (e.g. from "complex" to "double"
    ## or from "integer" to "logical") can introduce zeros. In this case,
    ## we must switch the type early. Otherwise we will end up with zeros
    ## in the "x" slot of the resulting dgCMatrix or lgCMatrix object.
    switch_type_early <- coercion_can_introduce_zeros(type(from), to_type)
    if (switch_type_early)
        type(from) <- to_type  # early type switching

    i <- from@nzcoo[ , 1L]
    j <- from@nzcoo[ , 2L]
    nzdata <- from@nzdata

    ## This type switching is safe only if it does not introduce zeros.
    if (!switch_type_early)
        storage.mode(nzdata) <- to_type  # late type switching

    if (orientation == "C") {
        CsparseMatrix(dim(from), i, j, nzdata, dimnames=dimnames(from))
    } else {
        RsparseMatrix(dim(from), i, j, nzdata, dimnames=dimnames(from))
    }
}

.from_COO_SparseMatrix_to_dgCMatrix <- function(from)
    .make_sparseMatrix_from_COO_SparseMatrix(from, "double", "C")

.from_COO_SparseMatrix_to_lgCMatrix <- function(from)
    .make_sparseMatrix_from_COO_SparseMatrix(from, "logical", "C")

setAs("COO_SparseMatrix", "dgCMatrix", .from_COO_SparseMatrix_to_dgCMatrix)
setAs("COO_SparseMatrix", "lgCMatrix", .from_COO_SparseMatrix_to_lgCMatrix)

.from_COO_SparseMatrix_to_dgRMatrix <- function(from)
    .make_sparseMatrix_from_COO_SparseMatrix(from, "double", "R")

.from_COO_SparseMatrix_to_lgRMatrix <- function(from)
    .make_sparseMatrix_from_COO_SparseMatrix(from, "logical", "R")

setAs("COO_SparseMatrix", "dgRMatrix", .from_COO_SparseMatrix_to_dgRMatrix)
setAs("COO_SparseMatrix", "lgRMatrix", .from_COO_SparseMatrix_to_lgRMatrix)

.make_COO_SparseMatrix_from_dgCMatrix_or_lgCMatrix <-
    function(from, use.dimnames=TRUE)
{
    i <- from@i + 1L
    j <- rep.int(seq_len(ncol(from)), diff(from@p))
    ans_nzcoo <- cbind(i, j, deparse.level=0L)
    ans_dimnames <- if (use.dimnames) dimnames(from) else NULL
    COO_SparseArray(dim(from), ans_nzcoo, from@x, ans_dimnames, check=FALSE)
}

.make_COO_SparseMatrix_from_dgRMatrix_or_lgRMatrix <-
    function(from, use.dimnames=TRUE)
{
    i <- rep.int(seq_len(nrow(from)), diff(from@p))
    j <- from@j + 1L
    ans_nzcoo <- cbind(i, j, deparse.level=0L)
    ans_dimnames <- if (use.dimnames) dimnames(from) else NULL
    COO_SparseArray(dim(from), ans_nzcoo, from@x, ans_dimnames, check=FALSE)
}

setAs("dgCMatrix", "COO_SparseMatrix",
    function(from) .make_COO_SparseMatrix_from_dgCMatrix_or_lgCMatrix(from)
)
setAs("lgCMatrix", "COO_SparseMatrix",
    function(from) .make_COO_SparseMatrix_from_dgCMatrix_or_lgCMatrix(from)
)
setAs("dgRMatrix", "COO_SparseMatrix",
    function(from) .make_COO_SparseMatrix_from_dgRMatrix_or_lgRMatrix(from)
)
setAs("lgRMatrix", "COO_SparseMatrix",
    function(from) .make_COO_SparseMatrix_from_dgRMatrix_or_lgRMatrix(from)
)

setAs("Matrix", "COO_SparseArray", function(from) as(from, "COO_SparseMatrix"))

