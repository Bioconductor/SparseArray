### =========================================================================
### The nz*() functions
### -------------------------------------------------------------------------
###
### A set of generic functions for direct manipulation of the nonzero
### elements of an array-like object.
###


### Works on any vector-like or array-like object 'x' that supports 'type(x)'
### and comparison with zero ('x != zero'). One notable exception are
### ngRMatrix objects from the Matrix package because 'x != FALSE' is
### broken on these objects.
.is_nonzero <- function(x)
{
    ## Make sure to use 'type()' and not 'typeof()'.
    zero <- vector(type(x), length=1L)
    is_nonzero <- x != zero  # broken on ngRMatrix objects!
    is_nonzero | is.na(is_nonzero)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### nzcount()
###
### Returns the number of nonzero array elements in 'x'.

setGeneric("nzcount", function(x) standardGeneric("nzcount"))

setMethod("nzcount", "ANY", function(x) sum(.is_nonzero(x)))

### Not 100% reliable on [d|l]gCMatrix objects because these objects are
### allowed to have zeros in their @x slot!
### See src/SVT_SparseArray_class.c for an example.
setMethod("nzcount", "CsparseMatrix", function(x) length(x@i))
setMethod("nzcount", "RsparseMatrix", function(x) length(x@j))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### nzwhich()
###
### Returns the indices of the nonzero array elements in 'x', either as
### an L-index (if 'arr.ind=FALSE') or as an M-index (if 'arr.ind=TRUE').
### The indices must be returned sorted in strictly ascending order.
### Note that using 'arr.ind=TRUE' won't work if 'nzcount(x)' is >= 2^31.

setGeneric("nzwhich", signature="x",
    function(x, arr.ind=FALSE) standardGeneric("nzwhich")
)

default_nzwhich <- function(x, arr.ind=FALSE)
{
    if (!isTRUEorFALSE(arr.ind))
        stop(wmsg("'arr.ind' must be TRUE or FALSE"))
    which(.is_nonzero(x), arr.ind=arr.ind, useNames=FALSE)
}

setMethod("nzwhich", "ANY", default_nzwhich)

### default_nzwhich() above works fine on [C|R]sparseMatrix derivatives
### (except on ngRMatrix objects) but nzwhich_CsparseMatrix() and
### nzwhich_RsparseMatrix() below are more efficient:
### - nzwhich_CsparseMatrix() is typically 10x to 50x faster than
###   default_nzwhich() on big CsparseMatrix objects;
### - nzwhich_RsparseMatrix() is only slightly faster than default_nzwhich()
###   but is provided mostly to make nzwhich() work on ngRMatrix objects
###   (default_nzwhich() doesn't work on these objects, see above).
### IMPORTANT NOTE: nzwhich_CsparseMatrix() and nzwhich_RsparseMatrix()
### use a shortcut that is **NOT** 100% reliable on [d|l]gCMatrix or
### [d|l]gRMatrix objects because these objects sometimes have zeros in
### their @x slot (see src/SVT_SparseArray_class.c for examples of such
### objects). So in this case the functions will produce a result that
### contains some false positives.
nzwhich_CsparseMatrix <- function(x, arr.ind=FALSE)
{
    if (!isTRUEorFALSE(arr.ind))
        stop(wmsg("'arr.ind' must be TRUE or FALSE"))
    x_nrow <- nrow(x)
    x_ncol <- ncol(x)
    ## If 'x' is a "long object" (i.e. length(x) >= 2^31) then we'll get
    ## an integer overflow in 'x_nrow * (seq_len(x_ncol) - 1L)' below.
    ## Coercing 'x_nrow' to double avoids that.
    if (is.double(length(x)))
        x_nrow <- as.double(x_nrow)
    offsets <- rep.int(x_nrow * (seq_len(x_ncol) - 1L), diff(x@p))
    ans <- x@i + offsets + 1L
    if (!arr.ind)
        return(ans)
    Lindex2Mindex(ans, dim(x))
}

nzwhich_RsparseMatrix <- function(x, arr.ind=FALSE)
{
    if (!isTRUEorFALSE(arr.ind))
        stop(wmsg("'arr.ind' must be TRUE or FALSE"))
    x_nrow <- nrow(x)
    x_ncol <- ncol(x)
    if (is.double(length(x)))
        x_ncol <- as.double(x_ncol)  # see nzwhich_CsparseMatrix() above
    offsets <- rep.int(x_ncol * (seq_len(x_nrow) - 1L), diff(x@p))
    transposed_nzwhich <- x@j + offsets + 1L
    transposed_arr_ind <- Lindex2Mindex(transposed_nzwhich, rev(dim(x)))
    transposed_arr_ind1 <- transposed_arr_ind[ , 1L]
    transposed_arr_ind2 <- transposed_arr_ind[ , 2L]
    oo <- order(transposed_arr_ind1, transposed_arr_ind2)
    arr_ind1 <- transposed_arr_ind2[oo]
    arr_ind2 <- transposed_arr_ind1[oo]
    arr_ind <- cbind(arr_ind1, arr_ind2, deparse.level=0)
    if (arr.ind)
        return(arr_ind)
    Mindex2Lindex(arr_ind, dim(x))
}

setMethod("nzwhich", "CsparseMatrix", nzwhich_CsparseMatrix)
setMethod("nzwhich", "RsparseMatrix", nzwhich_RsparseMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### nzvals()
###
### Extracts the nonzero values from 'x'. The nonzero values are returned
### in a vector of the same type() as 'x' and parallel to nzwhich(x).
### Equivalent to 'x[nzwhich(x)]' (and that's what the default method
### below does). However, specialized methods have the potential to make
### this dramatically faster.

setGeneric("nzvals", function(x) standardGeneric("nzvals"))

### Assumes that array-like object 'x' supports subsetting by a linear
### numeric subscript (L-index).
### Notes:
### - Subsetting by an M-index i.e. subsetting by 'nzwhich(x, arr.ind=TRUE)'
###   instead of 'nzwhich(x)' would work too but only when 'nzcount(x)'
###   is < 2^31.
### - The call to as.vector() should not be necessary since the result of
###   subsetting 'x' by an L- or M-index should already be a vector (atomic
###   or list). However, in the case of a 1D ordinary array, it's not!
###   For example, 'array(11:15)[matrix(3:1)]' is a 1D array.
###   Consider this a bug in base::`[`.
setMethod("nzvals", "ANY", function(x) as.vector(x[nzwhich(x)]))

### Not 100% reliable because [d|l]gCMatrix objects can have zeros in
### their @x slot (see comment for nzwhich_CsparseMatrix() above), but
### consistent with nzwhich_CsparseMatrix().
setMethod("nzvals", "dgCMatrix", function(x) x@x)
setMethod("nzvals", "lgCMatrix", function(x) x@x)

setMethod("nzvals", "ngCMatrix", function(x) rep.int(TRUE, length(x@i)))
setMethod("nzvals", "ngRMatrix", function(x) rep.int(TRUE, length(x@j)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### `nzvals<-`()
###
### Replaces the values of the nonzero array elements in 'x'.
### Equivalent to 'x[nzwhich(x)] <- value' (and that's what the default
### method below does). However specialized methods have the potential to
### make this dramatically faster.
### Notes:
### - 'nzvals(x) <- nzvals(x)' should **always** be a no-op.
### - 'nzvals(x) <- 0L' will wipe out all nonzero values from 'x'.

setGeneric("nzvals<-", signature="x",
    function(x, value) standardGeneric("nzvals<-")
)

setReplaceMethod("nzvals", "ANY",
    function(x, value)
    {
        if (!is.vector(value))
            stop(wmsg("replacement value must be a vector"))
        x[nzwhich(x)] <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sparsity()
###

sparsity <- function(x) { 1 - nzcount(x) / length(x) }

