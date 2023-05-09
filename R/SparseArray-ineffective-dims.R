### =========================================================================
### Drop/add ineffective dims from/to a SparseArray object
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### drop()
###

### Always returns an SVT_SparseArray object (endomorphism).
.drop_SVT_SparseArray <- function(x)
{
    stopifnot(is(x, "SVT_SparseArray"))
    ## Returns 'ans_dim', 'ans_dimnames', and 'ans_SVT', in a list of length 3.
    C_ans <- .Call2("C_drop_SVT_SparseArray_ineffective_dims",
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

