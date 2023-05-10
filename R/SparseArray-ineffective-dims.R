### =========================================================================
### Drop/add ineffective dims from/to a SparseArray object
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .tune_dims() and .tune_dimnames()
###
### Unlike .tune_SVT_SparseArray_dims() below in this file, .tune_dims()
### and .tune_dimnames() both accept a 'dim_tuner' that is not normalized.
### See src/SparseArray_ineffective_dims.c for more information.

.tune_dims <- function(dim, dim_tuner)
{
    stopifnot(is.integer(dim),
              is.logical(dim_tuner))
    .Call2("C_tune_dims", dim, dim_tuner, PACKAGE="SparseArray")
}

.tune_dimnames <- function(dimnames, dim_tuner)
{
    stopifnot(is.null(dimnames) || is.list(dimnames),
              is.logical(dim_tuner))
    .Call2("C_tune_dimnames", dimnames, dim_tuner, PACKAGE="SparseArray")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .tune_SVT_SparseArray_dims()
###
### Workhorse behind dim() setter for SVT_SparseArray objects.

.tune_SVT_SparseArray_dims <- function(x, dim_tuner)
{
    stopifnot(is(x, "SVT_SparseArray"),
              is.logical(dim_tuner))

    ans_SVT <- .Call2("C_tune_SVT_dims",
                      x@dim, x@type, x@SVT, dim_tuner,
                      PACKAGE="SparseArray")
    ans_dim <- .tune_dims(x@dim, dim_tuner)
    ans_dimnames <- .tune_dimnames(x@dimnames, x@dim, dim_tuner)

    new_SVT_SparseArray(ans_dim, ans_dimnames, x@type, ans_SVT, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### drop()
###

### Always returns an SVT_SparseArray object (endomorphism).
.drop_SVT_SparseArray <- function(x)
{
    stopifnot(is(x, "SVT_SparseArray"))
    ## Returns 'ans_dim', 'ans_dimnames', and 'ans_SVT', in a list of length 3.
    C_ans <- .Call2("C_drop_SVT_ineffective_dims",
                    x@dim, x@dimnames, x@type, x@SVT, PACKAGE="SparseArray")
    ans_dim <- C_ans[[1L]]
    ans_dimnames <- C_ans[[2L]]
    ans_SVT <- C_ans[[3L]]
    new_SVT_SparseArray(ans_dim, ans_dimnames, x@type, ans_SVT, check=FALSE)
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

