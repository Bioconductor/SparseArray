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
### .subset_NaSVT_by_logical_array()
###
### Returns a vector (atomic or list) of the same type() as 'x'.
###

.subset_NaSVT_by_logical_array <- function(x, y)
    stop("subsetting operation not supported yet")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .subset_NaSVT_by_Lindex()
### .subset_NaSVT_by_Mindex()
###
### Both return a vector (atomic or list) of the same type() as 'x'.
###

.propagate_names_if_1D <- function(ans, x_dimnames, index)
{
    if (length(x_dimnames) != 1L)
        return(ans)
    stopifnot(is.list(x_dimnames))
    x_names <- x_dimnames[[1L]]
    if (is.null(x_names))
        return(ans)
    stopifnot(is.character(x_names),
              identical(length(ans), length(index)))
    setNames(ans, x_names[index])
}

### 'Lindex' must be a numeric vector (integer or double), possibly a long one.
### NA indices are accepted.
.subset_NaSVT_by_Lindex <- function(x, Lindex)
{
    stopifnot(is(x, "NaArray"))
    check_svt_version(x)
    stopifnot(is.vector(Lindex), is.numeric(Lindex))
    on.exit(free_global_OPBufTree())
    ans <- SparseArray.Call("C_subset_NaSVT_by_Lindex",
                            x@dim, x@type, x@NaSVT, Lindex)
    .propagate_names_if_1D(ans, dimnames(x), Lindex)
}

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
    ans <- SparseArray.Call("C_subset_NaSVT_by_Mindex",
                            x@dim, x@type, x@NaSVT, Mindex)
    .propagate_names_if_1D(ans, x_dimnames, Mindex)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subset_NaSVT_by_Nindex()
###
### In addition to being one of the workhorses behind `[` on an
### NaArray object (see below), this is **the** workhorse behind the
### extract_na_array() and extract_array() methods for NaArray objects.
###
### 'Nindex' must be an N-index, that is, a list of numeric vectors (or NULLs),
### one along each dimension in the array to subset. Note that, strictly
### speaking, the vectors in an N-index are expected to be integer vectors,
### but subset_NaSVT_by_Nindex() can handle subscripts of type "double".
### This differs from the 'index' argument in 'extract_array()' where the
### subscripts **must** be integer vectors.
###
### Returns an NaArray object of the same type() as 'x' (endomorphism).

subset_NaSVT_by_Nindex <- function(x, Nindex, ignore.dimnames=FALSE)
{
    stopifnot(is(x, "NaArray"),
              is.list(Nindex),
              length(Nindex) == length(x@dim),
              isTRUEorFALSE(ignore.dimnames))
    check_svt_version(x)

    ## Returns 'new_dim' and 'new_NaSVT' in a list of length 2.
    C_ans <- SparseArray.Call("C_subset_SVT_by_Nindex",
                              x@dim, x@type, x@NaSVT, Nindex)
    new_dim <- C_ans[[1L]]
    new_NaSVT <- C_ans[[2L]]

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
### Single-bracket subsetting method (`[`) for NaArray objects
###

.subset_NaArray <- function(x, i, j, ..., drop=TRUE)
{
    if (missing(x))
        stop(wmsg("'x' is missing"))
    if (!isTRUEorFALSE(drop))
        stop(wmsg("'drop' must be TRUE or FALSE"))
    Nindex <- S4Arrays:::extract_Nindex_from_syscall(sys.call(), parent.frame())
    nsubscript <- length(Nindex)
    if (nsubscript == 0L)
        return(x)  # no-op
    x_dim <- dim(x)
    if (nsubscript == 1L && drop) {
        i <- Nindex[[1L]]
        if (type(i) == "logical" && identical(x_dim, dim(i)))
            return(.subset_NaSVT_by_logical_array(x, i))
        if (is.matrix(i))
            return(.subset_NaSVT_by_Mindex(x, i))
        if (is.numeric(i))
            return(.subset_NaSVT_by_Lindex(x, i))
    }
    if (nsubscript != length(x_dim))
        stop(wmsg("incorrect number of subscripts"))
    ## Note that this normalization will coerce the numeric subscripts
    ## in 'Nindex' to integer. This no longer necessary because now
    ## subset_NaSVT_by_Nindex() can handle subscripts of type "double" at
    ## the C level.
    ## TODO: Consider using a normalization process here that preserves
    ## the numeric subscripts.
    Nindex <- S4Arrays:::normalize_Nindex(Nindex, x)
    ans <- subset_NaSVT_by_Nindex(x, Nindex)
    if (drop)
        ans <- drop(ans)
    ans
}

setMethod("[", "NaArray", .subset_NaArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_na_array() and extract_array() methods for NaArray objects
###

setGeneric("extract_na_array", signature="x",
    function(x, index) standardGeneric("extract_na_array")
)

### No need to propagate the dimnames.
setMethod("extract_na_array", "NaArray",
    function(x, index) subset_NaSVT_by_Nindex(x, index, ignore.dimnames=TRUE)
)

### Note that the default extract_array() method would do the job but it
### relies on single-bracket subsetting so would needlessly go thru the
### complex .subset_NaArray() machinery above to finally call
### subset_NaSVT_by_Nindex(). It would also propagate the dimnames which
### extract_array() does not need to do. The method below completely bypasses
### all this complexity by calling subset_NaSVT_by_Nindex() directly.
setMethod("extract_array", "NaArray",
    function(x, index)
        as.array(subset_NaSVT_by_Nindex(x, index, ignore.dimnames=TRUE))
)

