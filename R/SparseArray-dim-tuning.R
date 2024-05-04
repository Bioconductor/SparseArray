### =========================================================================
### Dim tuning of a SparseArray object
### -------------------------------------------------------------------------
###
### See R/dim-tuning-utils.R in the S4Arrays package for more information
### about "dim tuning" and the tune_Array_dims() generic.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### tune_Array_dims() method for SVT_SparseArray objects
###
### This is the workhorse behind drop() and dim<-() on SVT_SparseArray
### objects.
###
### Unlike with S4Arrays:::tune_dims() and S4Arrays:::tune_dimnames(),
### the 'dim_tuner' vector passed to .tune_SVT_SparseArray_dims() must
### be normalized. See src/SparseArray_dim_tuning.c for more information.

.tune_SVT_SparseArray_dims <- function(x, dim_tuner)
{
    stopifnot(is(x, "SVT_SparseArray"),
              is.integer(dim_tuner))
    check_svt_version(x)

    ans_SVT <- SparseArray.Call("C_tune_SVT_dims",
                                x@dim, x@type, x@SVT, dim_tuner)
    ans_dim <- S4Arrays:::tune_dims(x@dim, dim_tuner)
    ans_dimnames <- S4Arrays:::tune_dimnames(x@dimnames, dim_tuner)

    new_SVT_SparseArray(ans_dim, ans_dimnames, x@type, ans_SVT, check=FALSE)
}

setMethod("tune_Array_dims", "SVT_SparseArray", .tune_SVT_SparseArray_dims)

