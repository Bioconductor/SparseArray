### =========================================================================
### read_block_as_sparse()
### -------------------------------------------------------------------------
###

### Like read_block_as_dense() in S4Arrays, read_block_as_sparse() is
### not meant to be called directly by the end user. They should call
### higher-level user-facing S4Arrays::read_block() function instead,
### with the 'as.sparse' argument set to TRUE.
### Also, like extract_sparse_array() that it is based on,
### read_block_as_sparse() should **always** be called on an array-like
### object 'x' for which 'is_sparse(x)' is TRUE. This is considered
### responsibility of the caller and we trust the caller to do the right
### thing, so, for the sake of efficiency, individual read_block_as_sparse()
### methods don't need to check this again. See extract_sparse_array.R
### for more information.
### Must also return a SparseArray object, like extract_sparse_array().
### Note that the S4Arrays::read_block() wrapper will take care of
### propagating the dimnames, so, for the sake of efficiency, individual
### methods should not try to do it.

setGeneric("read_block_as_sparse", signature="x",
    function(x, viewport) standardGeneric("read_block_as_sparse")
)

### Does NOT propagate the dimnames.
setMethod("read_block_as_sparse", "ANY",
    function(x, viewport)
    {
        Nindex <- makeNindexFromArrayViewport(viewport, expand.RangeNSBS=TRUE)
        extract_sparse_array(x, Nindex)
    }
)

