### =========================================================================
### The nna*() functions
### -------------------------------------------------------------------------
###
### A set of generic functions for direct manipulation of the non-NA
### elements of an array-like object.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### is_nonna()
###
### Returns an array-like object of type "logical".

setGeneric("is_nonna", function(x) standardGeneric("is_nonna"))

### Works on any vector-like or array-like object that supports is.na().
setMethod("is_nonna", "ANY", function(x) !is.na(x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### nnacount()
###
### Returns the number of non-NA array elements in 'x'.

setGeneric("nnacount", function(x) standardGeneric("nnacount"))

setMethod("nnacount", "ANY", function(x) sum(is_nonna(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### nnawhich()
###
### Returns the indices of the non-NA array elements in 'x', either as
### an L-index (if 'arr.ind=FALSE') or as an M-index (if 'arr.ind=TRUE').
### The indices must be returned sorted in strictly ascending order.
### Note that using 'arr.ind=TRUE' won't work if 'nnacount(x)' is >= 2^31.

setGeneric("nnawhich", signature="x",
    function(x, arr.ind=FALSE) standardGeneric("nnawhich")
)

.default_nnawhich <- function(x, arr.ind=FALSE)
{
    if (!isTRUEorFALSE(arr.ind))
        stop(wmsg("'arr.ind' must be TRUE or FALSE"))
    which(is_nonna(x), arr.ind=arr.ind, useNames=FALSE)
}

setMethod("nnawhich", "ANY", .default_nnawhich)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### nnavals()
###
### Extract the non-NA values from 'x'. The non-NA values are returned
### in a vector of the same type() as 'x' and parallel to nnawhich(x).
### Equivalent to 'x[nnawhich(x)]' (and that's what the default method
### below does). However, specialized methods have the potential to make
### this dramatically faster.

setGeneric("nnavals", function(x) standardGeneric("nnavals"))

### Assumes that array-like object 'x' supports subsetting by a linear
### numeric subscript (L-index).
### Notes:
### - Subsetting by an M-index i.e. subsetting by 'nnawhich(x, arr.ind=TRUE)'
###   instead of 'nnawhich(x)' would work too but only when 'nnacount(x)'
###   is < 2^31.
### - The call to as.vector() should not be necessary since the result of
###   subsetting 'x' by an L- or M-index should already be a vector (atomic
###   or list). However, in the case of a 1D ordinary array, it's not!
###   For example, 'array(11:15)[matrix(3:1)]' is a 1D array.
###   Consider this a bug in base::`[`.
### TODO: Maybe change this to 'x[is_nonna(x)]'? But do some timings first.
setMethod("nnavals", "ANY", function(x) as.vector(x[nnawhich(x)]))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### `nnavals<-`()
###
### Replace the values of the non-NA array elements in 'x'.
### Equivalent to 'x[nnawhich(x)] <- value' (and that's what the default
### method below does). However specialized methods have the potential to
### make this dramatically faster.
### Notes:
### - 'nnavals(x) <- nnavals(x)' should **always** be a no-op.
### - 'nnavals(x) <- NA' will wipe out all non-NA values from 'x'.

setGeneric("nnavals<-", signature="x",
    function(x, value) standardGeneric("nnavals<-")
)

### TODO: Maybe change this to 'x[is_nonna(x)] <- value'? But do some
### timings first.
setReplaceMethod("nnavals", "ANY",
    function(x, value)
    {
        if (!is.vector(value))
            stop(wmsg("replacement value must be a vector"))
        x[nnawhich(x)] <- value
        x
    }
)

