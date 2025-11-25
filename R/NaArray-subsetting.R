### =========================================================================
### NaArray subsetting
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### tune_Array_dims() method NaArray objects
###
### This is the workhorse behind drop() and dim<-() on NaArray objects.
###
### Unlike with S4Arrays:::tune_dims() and S4Arrays:::tune_dimnames(),
### the 'dim_tuner' vector passed to .tune_NaArray_dims() must be
### normalized. See src/SparseArray_dim_tuning.c for more information.

.tune_NaArray_dims <- function(x, dim_tuner)
{
    stopifnot(is(x, "NaArray"), is.integer(dim_tuner))
    check_svt_version(x)

    ans_NaSVT <- SparseArray.Call("C_tune_SVT_dims",
                                  x@dim, x@type, x@NaSVT, dim_tuner)
    ans_dim <- S4Arrays:::tune_dims(x@dim, dim_tuner)
    ans_dimnames <- S4Arrays:::tune_dimnames(x@dimnames, dim_tuner)

    new_NaArray(ans_dim, ans_dimnames, x@type, ans_NaSVT, check=FALSE)
}

setMethod("tune_Array_dims", "NaArray", .tune_NaArray_dims)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .subset_NaSVT_by_Lindex()
### .subset_NaSVT_by_Mindex()
###
### Both return a vector (atomic or list) of the same type() as 'x'.
###

### 'Lindex' must be a numeric vector (integer or double), possibly a long one.
### NA indices are accepted.
.subset_NaSVT_by_Lindex <- function(x, Lindex)
{
    stopifnot(is(x, "NaArray"))
    check_svt_version(x)
    stopifnot(is.vector(Lindex), is.numeric(Lindex))
    on.exit(free_global_OPBufTree())
    ans <- SparseArray.Call("C_subset_SVT_by_Lindex",
                            x@dim, x@type, x@NaSVT, TRUE, Lindex)
    propagate_names_if_1D(ans, dimnames(x), Lindex)
}

setMethod("subset_Array_by_Lindex", "NaArray", .subset_NaSVT_by_Lindex)

### Alright, '.subset_NaSVT_by_Mindex(x, Mindex)' could just have done:
###
###     .subset_NaSVT_by_Lindex(x, Mindex2Lindex(Mindex, dim(x)))
###
### However, the C code in C_subset_NaSVT_by_Mindex() avoids the Mindex2Lindex()
### step and so should be slightly more efficient, at least in theory. But is
### it? Some quick testing suggests that there's actually no significant
### difference!
### TODO: Investigate this more.
.subset_NaSVT_by_Mindex <- function(x, Mindex)
{
    stopifnot(is(x, "NaArray"))
    check_svt_version(x)
    stopifnot(is.matrix(Mindex))
    x_dimnames <- dimnames(x)
    if (!is.numeric(Mindex)) {
        if (!is.character(Mindex))
            stop(wmsg("invalid matrix subscript type \"", type(Mindex), "\""))
        if (is.null(x_dimnames))
            stop(wmsg("NaArray object to subset has no dimnames"))
        ## Subsetting an ordinary array with dimnames on it by a character
        ## matrix is supported in base R but we don't support this yet for
        ## NaArray objects.
        stop("subsetting an NaArray object by a character matrix ",
             "is not supported at the moment")
    }
    on.exit(free_global_OPBufTree())
    ans <- SparseArray.Call("C_subset_SVT_by_Mindex",
                            x@dim, x@type, x@NaSVT, TRUE, Mindex)
    propagate_names_if_1D(ans, x_dimnames, Mindex)
}

setMethod("subset_Array_by_Mindex", "NaArray", .subset_NaSVT_by_Mindex)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .subset_NaSVT_as_Rarray()
###
### Equivalent to 'as.array(.subset_NaSVT_as_NaSVT(x, Nindex))' but slightly
### more efficient.
###
### TODO: Add unit tests for this.

.subset_NaSVT_as_Rarray <- function(x, Nindex, ignore.dimnames=FALSE)
{
    stopifnot(is(x, "NaArray"), isTRUEorFALSE(ignore.dimnames))
    check_svt_version(x)

    new_dim <- S4Arrays:::get_Nindex_lengths(Nindex, x@dim)
    Noffs <- Nindex2Noffs(Nindex)
    ans <- SparseArray.Call("C_subset_SVT_as_Rarray",
                            x@dim, x@type, x@NaSVT, TRUE, Noffs)
    if (!ignore.dimnames)
        dimnames(ans) <- S4Arrays:::subset_dimnames_by_Nindex(x@dimnames,
                                                              Nindex)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .subset_NaSVT_as_NaSVT()
###
### Returns an NaArray object of the same type() as 'x' (endomorphism).

.subset_NaSVT_as_NaSVT <- function(x, Nindex, ignore.dimnames=FALSE)
{
    stopifnot(is(x, "NaArray"), isTRUEorFALSE(ignore.dimnames))
    check_svt_version(x)

    new_dim <- S4Arrays:::get_Nindex_lengths(Nindex, x@dim)
    Noffs <- Nindex2Noffs(Nindex)
    new_NaSVT <- SparseArray.Call("C_subset_SVT_as_SVT",
                                  x@dim, x@type, x@NaSVT, Noffs)

    ## Compute 'new_dimnames'.
    if (is.null(dimnames(x)) || ignore.dimnames) {
        new_dimnames <- vector("list", length(x@dim))
    } else {
        new_dimnames <- S4Arrays:::subset_dimnames_by_Nindex(x@dimnames, Nindex)
    }
    BiocGenerics:::replaceSlots(x, dim=new_dim,
                                   dimnames=new_dimnames,
                                   NaSVT=new_NaSVT,
                                   check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subset_Array_by_Nindex(), extract_array(), and extract_na_array() methods
### for NaArray objects
###

setMethod("subset_Array_by_Nindex", "NaArray",
    function(x, Nindex, drop=TRUE)
    {
        ans_dim <- S4Arrays:::get_Nindex_lengths(Nindex, x@dim)
        if (drop && sum(ans_dim != 1L) <= 1L) {
            ans <- .subset_NaSVT_as_Rarray(x, Nindex)
            return(S4Arrays:::drop_even_if_1D(ans))
        }
        ans <- .subset_NaSVT_as_NaSVT(x, Nindex)
        if (drop)
            ans <- drop(ans)
        ans
    }
)

setMethod("extract_array", "NaArray",
    function(x, index) .subset_NaSVT_as_Rarray(x, index, ignore.dimnames=TRUE)
)

setGeneric("extract_na_array", signature="x",
    function(x, index) standardGeneric("extract_na_array")
)

### No need to propagate the dimnames.
setMethod("extract_na_array", "NaArray",
    function(x, index) .subset_NaSVT_as_NaSVT(x, index, ignore.dimnames=TRUE)
)

