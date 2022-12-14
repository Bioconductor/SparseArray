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
### - Getters nzcoo() and nzvals()
### - Coercion from array to COO_SparseArray
### - Back and forth coercion between COO_SparseArray and [d|l]g[C|R]Matrix
###   objects from the Matrix package
###

setClass("COO_SparseArray",
    contains="SparseArray",
    representation(
        nzcoo="matrix",  # M-index containing the coordinates of the
                         # nonzero values.
        nzvals="vector"  # A vector (atomic or list) of length
                         # 'nrow(nzcoo)' containing the nonzero values.
    ),
    prototype(
        nzcoo=matrix(integer(0), ncol=1L),
        nzvals=logical(0)
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

### Automatic coercion method from COO_SparseArray to COO_SparseMatrix silently
### returns a broken object (unfortunately these dummy automatic coercion
### methods don't bother to validate the object they return). So we overwrite
### it with a method that will fail (as expected) thanks to the validity
### method for SparseMatrix objects.
setAs("COO_SparseArray", "COO_SparseMatrix",
    function(from) new("COO_SparseMatrix", from)
)

### Also the user should not be able to degrade a COO_SparseMatrix object to a
### COO_SparseArray object so 'as(x, "COO_SparseArray", strict=TRUE)' should
### fail or be a no-op when 'x' is a COO_SparseMatrix object. Making this
### coercion a no-op seems to be the easiest (and safest) way to go.
setAs("COO_SparseMatrix", "COO_SparseArray", function(from) from)  # no-op


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

.validate_nzvals_slot <- function(x)
{
    x_nzvals <- x@nzvals
    if (!(is.vector(x_nzvals) && length(x_nzvals) == nrow(x@nzcoo)))
        return(paste0("'nzvals' slot must be a vector of length ",
                      "the number of rows in the 'nzcoo' slot"))
    TRUE
}

.validate_COO_SparseArray <- function(x)
{
    msg <- .validate_nzcoo_slot(x)
    if (!isTRUE(msg))
        return(msg)
    msg <- .validate_nzvals_slot(x)
    if (!isTRUE(msg))
        return(msg)
    TRUE
}
setValidity2("COO_SparseArray", .validate_COO_SparseArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setMethod("type", "COO_SparseArray", function(x) type(x@nzvals))

setGeneric("nzcoo", function(x) standardGeneric("nzcoo"))
setMethod("nzcoo", "COO_SparseArray", function(x) x@nzcoo)

setGeneric("nzvals", function(x) standardGeneric("nzvals"))
setMethod("nzvals", "COO_SparseArray", function(x) x@nzvals)


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

    new_nzvals <- x@nzvals
    storage.mode(new_nzvals) <- value
    nzidx <- which_is_nonzero(new_nzvals)
    new_nzcoo <- x@nzcoo[nzidx, , drop=FALSE]
    new_nzvals <- new_nzvals[nzidx]
    BiocGenerics:::replaceSlots(x, nzcoo=new_nzcoo,
                                   nzvals=new_nzvals,
                                   check=FALSE)
}

setReplaceMethod("type", "COO_SparseArray", .set_COO_SparseArray_type)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### nzcount()
###

setMethod("nzcount", "COO_SparseArray", function(x) length(nzvals(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level constructor
###

new_COO_SparseArray <- function(dim, dimnames=NULL,
                                nzcoo=NULL, nzvals=NULL, check=TRUE)
{
    stopifnot(is.integer(dim))
    if (length(dim) == 2L) {
        ans_class <- "COO_SparseMatrix"
    } else {
        ans_class <- "COO_SparseArray"
    }
    dimnames <- S4Arrays:::normarg_dimnames(dimnames, dim)
    new2(ans_class, dim=dim, dimnames=dimnames,
                    nzcoo=nzcoo, nzvals=nzvals, check=check)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

.normarg_nzvals <- function(nzvals, length.out)
{
    if (is.null(nzvals))
        stop(wmsg("'nzvals' cannot be NULL when 'nzcoo' is not NULL"))
    if (!is.vector(nzvals))
        stop(wmsg("'nzvals' must be a vector"))
    ## Same logic as S4Vectors:::V_recycle().
    nzvals_len <- length(nzvals)
    if (nzvals_len == length.out)
        return(nzvals)
    if (nzvals_len > length.out && nzvals_len != 1L)
        stop(wmsg("'length(nzvals)' is greater than 'nrow(nzcoo)'"))
    if (nzvals_len == 0L)
        stop(wmsg("'length(nzvals)' is 0 but 'nrow(nzcoo)' is not"))
    if (length.out %% nzvals_len != 0L)
        warning(wmsg("'nrow(nzcoo)' is not a multiple of 'length(nzvals)'"))
    rep(nzvals, length.out=length.out)
}

COO_SparseArray <- function(dim, nzcoo=NULL, nzvals=NULL, dimnames=NULL,
                                 check=TRUE)
{
    if (!is.numeric(dim))
        stop(wmsg("'dim' must be an integer vector"))
    if (!is.integer(dim))
        dim <- as.integer(dim)
    if (is.null(nzcoo)) {
        if (is.null(nzvals)) {
            nzvals <- logical(0)  # vector()
        } else if (!(is.vector(nzvals) && length(nzvals) == 0L)) {
            stop(wmsg("'nzvals' must be NULL or a vector of length 0 ",
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
        nzvals <- .normarg_nzvals(nzvals, nrow(nzcoo))
    }
    new_COO_SparseArray(dim, dimnames, nzcoo, nzvals, check=check)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .dense2sparse() and .sparse2dense()
###

### 'x' must be an array-like object that supports 'type()' and subsetting
### by an M-index subscript.
### Returns a COO_SparseArray object.
.dense2sparse <- function(x)
{
    x_dim <- dim(x)
    if (is.null(x_dim))
        stop(wmsg("'x' must be an array-like object"))
    ans_nzcoo <- which_is_nonzero(x, arr.ind=TRUE)  # M-index
    ans_nzvals <- x[ans_nzcoo]
    ## Work around bug in base::`[`
    if (length(x_dim) == 1L)
        ans_nzvals <- as.vector(ans_nzvals)
    COO_SparseArray(x_dim, ans_nzcoo, ans_nzvals, dimnames(x), check=FALSE)
}

### 'coo' must be a COO_SparseArray object.
### Return an ordinary array.
.sparse2dense <- function(coo)
{
    if (!is(coo, "COO_SparseArray"))
        stop(wmsg("'coo' must be a COO_SparseArray object"))
    coo_nzvals <- nzvals(coo)
    zero <- vector(typeof(coo_nzvals), length=1L)
    ans <- array(zero, dim=dim(coo))
    ans[nzcoo(coo)] <- coo_nzvals
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

.make_sparseMatrix_from_COO_SparseMatrix <- function(from, to_type, form)
{
    stopifnot(is(from, "COO_SparseMatrix"))

    ## Coercion to dg[C|R]Matrix (i.e. to_type="double"):
    ## If 'type(from)' is "logical", "integer", or "raw", we'll coerce 'nzvals'
    ## to "double" right before passing it to [C|R]sparseMatrix() below. This
    ## is ok because it won't introduce zeros in 'nzvals'. Also it should be
    ## slightly more efficient than switching the type of 'from' now.
    ## However, if the coercion to "double" can potentially introduce zeros
    ## (e.g. if 'type(from)' is "complex"), then we need to switch the type
    ## now. Otherwise we will end up with zeros in the "x" slot of the
    ## resulting dg[C|R]Matrix object.

    ## Coercion to lg[C|R]Matrix (i.e. to_type="logical"):
    ## If 'type(from)' is "integer", "double", "complex", or "raw", we'll
    ## coerce 'nzvals' to "logical" right before passing it to
    ## [C|R]sparseMatrix() below. This is ok because it won't introduce
    ## logical zeros (i.e. FALSEs) in 'nzvals'. Also it should be slightly
    ## more efficient than switching the type of 'from' now.
    ## However, if the coercion to "logical" can potentially introduce zeros
    ## (e.g. if 'type(from)' is "character"), then we need to switch the type
    ## now. Otherwise we will end up with zeros in the "x" slot of the
    ## resulting lg[C|R]Matrix object.

    postpone <- coercion_can_introduce_zeros(type(from), to_type)
    if (!postpone)
        type(from) <- to_type  # early type switching

    i <- from@nzcoo[ , 1L]
    j <- from@nzcoo[ , 2L]
    nzvals <- from@nzvals

    ## This type switching is safe only if it does not introduce zeros.
    if (postpone)
        storage.mode(nzvals) <- to_type  # late type switching

    if (form == "C") {
        CsparseMatrix(dim(from), i, j, nzvals, dimnames=dimnames(from))
    } else {
        RsparseMatrix(dim(from), i, j, nzvals, dimnames=dimnames(from))
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### aperm()
###
### Like aperm2() in S4Arrays, extend base::aperm() by allowing dropping
### and/or adding ineffective dimensions.
###

.aperm.COO_SparseArray <- function(a, perm)
{
    a_dim <- dim(a)
    perm <- S4Arrays:::normarg_perm(perm, a_dim)
    msg <- S4Arrays:::validate_perm(perm, a_dim)
    if (!isTRUE(msg))
        stop(wmsg(msg))
    ans_dim <- a_dim[perm]
    ans_dim[is.na(perm)] <- 1L
    ans_nzcoo <- a@nzcoo[ , perm, drop=FALSE]
    ans_nzcoo[ , is.na(perm)] <- 1L
    ans_dimnames <- a@dimnames[perm]
    new_COO_SparseArray(ans_dim, ans_dimnames,
                        ans_nzcoo, a@nzvals, check=FALSE)
}

### S3/S4 combo for aperm.COO_SparseArray
aperm.COO_SparseArray <-
    function(a, perm, ...) .aperm.COO_SparseArray(a, perm, ...)
setMethod("aperm", "COO_SparseArray", aperm.COO_SparseArray)

