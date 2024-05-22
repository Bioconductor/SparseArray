### =========================================================================
### SparseArray subsetting
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .subset_SVT_by_logical_array()
###
### Returns a vector (atomic or list) of the same type() as 'x'.
###

.subset_SVT_by_logical_array <- function(x, y)
    stop("subsetting operation not supported yet")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .subset_SVT_by_Lindex()
### .subset_SVT_by_Mindex()
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
.subset_SVT_by_Lindex <- function(x, Lindex)
{
    stopifnot(is(x, "SVT_SparseArray"))
    check_svt_version(x)
    stopifnot(is.vector(Lindex), is.numeric(Lindex))
    on.exit(free_global_OPBufTree())
    ans <- SparseArray.Call("C_subset_SVT_by_Lindex",
                            x@dim, x@type, x@SVT, Lindex)
    .propagate_names_if_1D(ans, dimnames(x), Lindex)
}

### Alright, '.subset_SVT_by_Mindex(x, Mindex)' could just have done:
###
###     .subset_SVT_by_Lindex(x, Mindex2Lindex(Mindex, dim(x)))
###
### However, the C code in C_subset_SVT_by_Mindex() avoids the Mindex2Lindex()
### step and so should be slightly more efficient, at least in theory. But is
### it? Some quick testing suggests that there's actually no significant
### difference!
### TODO: Investigate this more.
.subset_SVT_by_Mindex <- function(x, Mindex)
{
    stopifnot(is(x, "SVT_SparseArray"))
    check_svt_version(x)
    stopifnot(is.matrix(Mindex))
    x_dimnames <- dimnames(x)
    if (!is.numeric(Mindex)) {
        if (!is.character(Mindex))
            stop(wmsg("invalid matrix subscript type \"", type(Mindex), "\""))
        if (is.null(x_dimnames))
            stop(wmsg("SparseArray object to subset has no dimnames"))
        ## Subsetting an ordinary array with dimnames on it by a character
        ## matrix is supported in R base but we don't support this yet for
        ## SparseArray objects.
        stop("subsetting a SparseArray object by a character matrix ",
             "is not supported at the moment")
    }
    on.exit(free_global_OPBufTree())
    ans <- SparseArray.Call("C_subset_SVT_by_Mindex",
                            x@dim, x@type, x@SVT, Mindex)
    .propagate_names_if_1D(ans, x_dimnames, Mindex)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subset_SVT_by_Nindex()
###
### In addition to being one of the workhorses behind `[` on an
### SVT_SparseArray object (see below), this is **the** workhorse behind the
### extract_sparse_array() and extract_array() methods for SVT_SparseArray
### objects.
###
### 'Nindex' must be an N-index, that is, a list of numeric vectors (or NULLs),
### one along each dimension in the array to subset. Note that, strictly
### speaking, the vectors in an N-index are expected to be integer vectors,
### but subset_SVT_by_Nindex() can handle subscripts of type "double".
### This differs from the 'index' argument in 'extract_array()' where the
### subscripts **must** be integer vectors.
###
### Returns an SVT_SparseArray object of the same type() as 'x' (endomorphism).

subset_SVT_by_Nindex <- function(x, Nindex, ignore.dimnames=FALSE)
{
    stopifnot(is(x, "SVT_SparseArray"),
              is.list(Nindex),
              length(Nindex) == length(x@dim),
              isTRUEorFALSE(ignore.dimnames))
    check_svt_version(x)

    ## Returns 'new_dim' and 'new_SVT' in a list of length 2.
    C_ans <- SparseArray.Call("C_subset_SVT_by_Nindex",
                              x@dim, x@type, x@SVT, Nindex)
    new_dim <- C_ans[[1L]]
    new_SVT <- C_ans[[2L]]

    ## Compute 'new_dimnames'.
    if (is.null(dimnames(x)) || ignore.dimnames) {
        new_dimnames <- vector("list", length(x@dim))
    } else {
        new_dimnames <- S4Arrays:::subset_dimnames_by_Nindex(x@dimnames, Nindex)
    }
    BiocGenerics:::replaceSlots(x, dim=new_dim,
                                   dimnames=new_dimnames,
                                   SVT=new_SVT,
                                   check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Single-bracket subsetting method (`[`) for SVT_SparseArray objects
###

.subset_SVT_SparseArray <- function(x, i, j, ..., drop=TRUE)
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
            return(.subset_SVT_by_logical_array(x, i))
        if (is.matrix(i))
            return(.subset_SVT_by_Mindex(x, i))
        if (is.numeric(i))
            return(.subset_SVT_by_Lindex(x, i))
    }
    if (nsubscript != length(x_dim))
        stop(wmsg("incorrect number of subscripts"))
    ## Note that this normalization will coerce the numeric subscripts
    ## in 'Nindex' to integer. This no longer necessary because now
    ## subset_SVT_by_Nindex() can handle subscripts of type "double" at
    ## the C level.
    ## TODO: Consider using a normalization process here that preserves
    ## the numeric subscripts.
    Nindex <- S4Arrays:::normalize_Nindex(Nindex, x)
    ans <- subset_SVT_by_Nindex(x, Nindex)
    if (drop)
        ans <- drop(ans)
    ans
}

setMethod("[", "SVT_SparseArray", .subset_SVT_SparseArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_sparse_array() and extract_array() methods for SVT_SparseArray
### objects
###

### No need to propagate the dimnames.
setMethod("extract_sparse_array", "SVT_SparseArray",
    function(x, index) subset_SVT_by_Nindex(x, index, ignore.dimnames=TRUE)
)

### Note that the default extract_array() method would do the job but it
### relies on single-bracket subsetting so would needlessly go thru the
### complex .subset_SVT_SparseArray() machinery above to finally call
### subset_SVT_by_Nindex(). It would also propagate the dimnames which
### extract_array() does not need to do. The method below completely bypasses
### all this complexity by calling subset_SVT_by_Nindex() directly.
setMethod("extract_array", "SVT_SparseArray",
    function(x, index)
        as.array(subset_SVT_by_Nindex(x, index, ignore.dimnames=TRUE))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_sparse_array() and extract_array() methods for COO_SparseArray
### objects
###

### IMPORTANT NOTE: The returned COO_SparseArray object is guaranteed to be
### **correct** ONLY if the subscripts in 'index' do NOT contain duplicates!
### If they contain duplicates, the correct COO_SparseArray object to return
### should contain repeated nonzero data. However, in order to keep it as
### efficient as possible, the code below does NOT repeat the nonzero data
### that corresponds to duplicates subscripts. It does not check for
### duplicates in 'index' either because this check could have a
### significant cost.
### All this is OK because .extract_COO_SparseArray_subset() should
### always be used in a context where 'index' does NOT contain duplicates.
### The only situation where 'index' CAN contain duplicates is when
### .extract_COO_SparseArray_subset() is called by
### .extract_array_from_COO_SparseArray(), in which case the
### missing nonzero data are added later.
.extract_COO_SparseArray_subset <- function(x, index)
{
    stopifnot(is(x, "COO_SparseArray"))
    ans_dim <- S4Arrays:::get_Nindex_lengths(index, dim(x))
    x_nzcoo <- x@nzcoo
    for (along in seq_along(ans_dim)) {
        i <- index[[along]]
        if (is.null(i))
            next
        x_nzcoo[ , along] <- match(x_nzcoo[ , along], i)
    }
    ## Note that calling rowAnyNAs() on ordinary matrix 'x_nzcoo' would
    ## also work as it would call the rowAnyNAs() S4 generic defined in
    ## MatrixGenerics, and the latter would eventually dispatch on
    ## matrixStats::rowAnyNAs(). However, calling matrixStats::rowAnyNAs()
    ## should be slightly more efficient. Also note that this call is the
    ## only reason why we list matrixStats in the Imports field.
    keep_idx <- which(!matrixStats::rowAnyNAs(x_nzcoo))
    ans_nzcoo <- x_nzcoo[keep_idx, , drop=FALSE]
    ans_nzdata <- x@nzdata[keep_idx]
    COO_SparseArray(ans_dim, ans_nzcoo, ans_nzdata, check=FALSE)
}
setMethod("extract_sparse_array", "COO_SparseArray",
    .extract_COO_SparseArray_subset
)

.extract_array_from_COO_SparseArray <- function(x, index)
{
    coo0 <- .extract_COO_SparseArray_subset(x, index)
    ## If the subscripts in 'index' contain duplicates, 'coo0' is
    ## "incomplete" in the sense that it does not contain the nonzero data
    ## that should have been repeated according to the duplicates in the
    ## subscripts (see IMPORTANT NOTE above).
    ans0 <- as.array(coo0)
    ## We "complete" 'ans0' by repeating the nonzero data according to the
    ## duplicates present in 'index'. Note that this is easy and cheap to
    ## do now because 'ans0' uses a dense representation (it's an ordinary
    ## array). This would be harder to do **natively** on the
    ## COO_SparseArray form (i.e. without converting to dense first
    ## then back to sparse).
    sm_index <- lapply(index,
        function(i) {
            if (is.null(i))
                return(NULL)
            sm <- match(i, i)
            if (isSequence(sm))
                return(NULL)
            sm
        })
    if (all(S4Vectors:::sapply_isNULL(sm_index)))
        return(ans0)
    S4Arrays:::subset_by_Nindex(ans0, sm_index)
}
setMethod("extract_array", "COO_SparseArray",
    .extract_array_from_COO_SparseArray
)

