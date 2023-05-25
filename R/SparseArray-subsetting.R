### =========================================================================
### SparseArray subsetting
### -------------------------------------------------------------------------


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
    ans_nzvals <- x@nzvals[keep_idx]
    COO_SparseArray(ans_dim, ans_nzcoo, ans_nzvals, check=FALSE)
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_sparse_array() method for SVT_SparseArray objects
###

### This is the workhorse behind the extract_sparse_array(), extract_array(),
### and `[` methods for SVT_SparseArray objects.
### Like for 'extract_array()', the supplied 'index' must be a list with
### one list element per dimension in 'x'. Each list element must be an
### integer vector of valid indices along the corresponding dimension
### in 'x', or a NULL.
subset_SVT_SparseArray <- function(x, index, ignore.dimnames=FALSE)
{
    stopifnot(is(x, "SVT_SparseArray"),
              is.list(index),
              length(index) == length(x@dim),
              isTRUEorFALSE(ignore.dimnames))

    ## Returns 'new_dim' and 'new_SVT' in a list of length 2.
    C_ans <- .Call2("C_subset_SVT_SparseArray",
                    x@dim, x@type, x@SVT, index, PACKAGE="SparseArray")
    new_dim <- C_ans[[1L]]
    new_SVT <- C_ans[[2L]]

    ## Compute 'new_dimnames'.
    if (is.null(dimnames(x)) || ignore.dimnames) {
        new_dimnames <- vector("list", length(x@dim))
    } else {
        new_dimnames <- S4Arrays:::subset_dimnames_by_Nindex(x@dimnames, index)
    }
    BiocGenerics:::replaceSlots(x, dim=new_dim,
                                   dimnames=new_dimnames,
                                   SVT=new_SVT,
                                   check=FALSE)
}

### No need to propagate the dimnames.
setMethod("extract_sparse_array", "SVT_SparseArray",
    function(x, index) subset_SVT_SparseArray(x, index, ignore.dimnames=TRUE)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SVT_SparseArray single-bracket subsetting
###

.subset_SVT_SparseArray_by_logical_array <- function(x, y, drop=TRUE)
    stop("subsetting operation not supported yet")

.subset_SVT_SparseArray_by_Mindex <- function(x, Mindex, drop=FALSE)
    stop("subsetting operation not supported yet")

.subset_SVT_SparseArray_by_Lindex <- function(x, Lindex)
    stop("subsetting operation not supported yet")

.single_bracket_SVT_SparseArray <- function(x, i, j, ..., drop=TRUE)
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
    x_ndim <- length(x_dim)
    if (nsubscript == 1L) {
        i <- Nindex[[1L]]
        if (type(i) == "logical" && identical(x_dim, dim(i)))
            return(.subset_SVT_SparseArray_by_logical_array(x, i, drop=drop))
        if (is.matrix(i) && is.numeric(i))
            return(.subset_SVT_SparseArray_by_Mindex(x, i, drop=drop))
        ## Linear single bracket subsetting e.g. x[5:2].
        ## If 'x' is monodimensional and 'drop' is FALSE, we fallback
        ## to "multi-dimensional single bracket subsetting" which is an
        ## endomorphism.
        if (x_ndim != 1L || drop)
            return(.subset_SVT_SparseArray_by_Lindex(x, i))
    }
    if (nsubscript != x_ndim)
        stop(wmsg("incorrect number of subscripts"))
    index <- S4Arrays:::normalize_Nindex(Nindex, x)
    ans <- subset_SVT_SparseArray(x, index)
    if (drop)
        ans <- drop(ans)
    ans
}

setMethod("[", "SVT_SparseArray", .single_bracket_SVT_SparseArray)

### Note that the default extract_array() method would do the job but it
### relies on single-bracket subsetting so would needlessly go thru the
### complex .single_bracket_SVT_SparseArray() machinery above. It would
### also propagate the dimnames which extract_array() does not need to do.
### The method below completely bypasses all this complexity.
setMethod("extract_array", "SVT_SparseArray",
    function(x, index)
        as.array(subset_SVT_SparseArray(x, index, ignore.dimnames=TRUE))
)

