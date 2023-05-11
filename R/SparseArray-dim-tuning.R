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
### This is the workhorse behind the drop() method and dim() setter for
### SVT_SparseArray objects.
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

### TODO: Maybe introduce tune_array_dims() generic (define it in
### S4Arrays/R/dim-tuning-utils.R) and make .tune_SVT_SparseArray_dims()
### the tune_array_dims() method for SVT_SparseArray objects.
### Then also define the drop() method and dim() replace method for Array
### objects (also in S4Arrays/R/dim-tuning-utils.R), like below, but have
### them use tune_array_dims() instead of .tune_SVT_SparseArray_dims().
### Advantage: other Array derivatives (e.g. DelayedArray objects) will
### only need to implement a tune_array_dims() method and this will give
### them drop() and the dim() setter for free.

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
    is_ineffective <- dim(x) == 1L
    ## We cannot drop all dimensions so we (arbitrarily) keep the first one.
    if (all(is_ineffective))
        is_ineffective[[1L]] <- FALSE
    dim_tuner <- -as.integer(is_ineffective)
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dim() setter
###

.set_SVT_SparseArray_dim <- function(x, value)
{
    stopifnot(is(x, "SVT_SparseArray"))
    x_dim <- dim(x)
    value <- S4Arrays:::normalize_dim_replacement_value(value, x_dim)
    dim_tuner <-
        S4Arrays:::make_dim_tuner_from_old2new_dims(x_dim, value, class(x))
    ans <- .tune_SVT_SparseArray_dims(x, dim_tuner)
    stopifnot(identical(dim(ans), value))  # sanity check
    ans
}

setReplaceMethod("dim", "SVT_SparseArray", .set_SVT_SparseArray_dim)

