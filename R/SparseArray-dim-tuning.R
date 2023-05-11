### =========================================================================
### Dim tuning of a SparseArray object
### -------------------------------------------------------------------------
###
### Dim tuning is the act of adding and/or dropping ineffective dimensions
### to/from an array-like object. The exact actions to perform on the
### dimensions of the object are described via the 'dim_tuner' argument.
### See src/dim_tuning_utils.c in the S4Arrays package and
### src/SparseArray_dim_tuning.c in this package for more information.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .tune_SVT_SparseArray_dims()
###
### Workhorse behind drop() method and dim() setter for SVT_SparseArray
### objects.
###
### Unlike with S4Arrays:::tune_dims() and S4Arrays:::tune_dimnames(),
### the 'dim_tuner' vector passed to .tune_SVT_SparseArray_dims() must
### be normalized. See src/SparseArray_dim_tuning.c for more information.
###
### To revert a dim tuning, simply tune again with '-dim_tuner' (i.e. minus
### 'dim_tuner'). More precisely, for .tune_SVT_SparseArray_dims() this
### means that 'svt2' will always be identical to 'svt' here:
###
###   tuned_svt <- SparseArray:::.tune_SVT_SparseArray_dims(svt, dim_tuner)
###   svt2 <- SparseArray:::.tune_SVT_SparseArray_dims(tuned_svt, -dim_tuner)
###
### This should be TRUE for any SVT_SparseArray object 'svt' (with no
### dimnames on its ineffective dimensions) and any 'dim_tuner' vector
### compatible with 'dim(svt)'.

.tune_SVT_SparseArray_dims <- function(x, dim_tuner)
{
    stopifnot(is(x, "SVT_SparseArray"),
              is.integer(dim_tuner))

    ans_SVT <- .Call2("C_tune_SVT_dims",
                      x@dim, x@type, x@SVT, dim_tuner,
                      PACKAGE="SparseArray")
    ans_dim <- S4Arrays:::tune_dims(x@dim, dim_tuner)
    ans_dimnames <- S4Arrays:::tune_dimnames(x@dimnames, dim_tuner)

    new_SVT_SparseArray(ans_dim, ans_dimnames, x@type, ans_SVT, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### drop()
###

### Always returns an SVT_SparseArray object (endomorphism).
.drop_SVT_SparseArray <- function(x)
{
    stopifnot(is(x, "SVT_SparseArray"))
    dim_tuner <- -as.integer(dim(x) == 1L)
    .tune_SVT_SparseArray_dims(x, dim_tuner)
}

### Returns an SVT_SparseArray object or an ordinary vector of type 'type(x)'.
### The dimnames are propagated except when 'length(x)' is 1 (i.e.
### 'all(dim(x) == 1)' is TRUE).
setMethod("drop", "SVT_SparseArray",
    function(x)
    {
        x <- .drop_SVT_SparseArray(x)
        if (length(dim(x)) != 1L)
            return(x)     # SVT_SparseArray object
        a <- as.array(x)  # 1d ordinary array
        ans <- as.vector(a)  # unfortunately, this drops the names
        ## Restore the names.
        a_dimnames <- dimnames(a)  # NULL or list of length 1
        if (!is.null(a_dimnames))
            names(ans) <- a_dimnames[[1L]]
        ans
    }
)

