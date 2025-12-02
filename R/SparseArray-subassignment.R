### =========================================================================
### SparseArray subassignment
### -------------------------------------------------------------------------
###


adjust_left_type <- function(x, value)
{
    stopifnot(is(x, "SVT_SparseArray") || is(x, "NaArray"))
    check_svt_version(x)
    if (!is.vector(value))
        stop(wmsg("the supplied value must be a vector for this form ",
                  "of subassignment to an SVT_SparseArray object"))
    ## Change 'x' type if necessary.
    new_type <- type(c(vector(type(x)), vector(type(value))))
    type(x) <- new_type
    x
}

### Adjust the type of 'value' and recycle it to the length of the
### subassignment M/L-index.
normalize_right_value <- function(value, left_type, index_len)
{
    if (length(value) == 0L)
        stop(wmsg("right value has length zero"))
    storage.mode(value) <- left_type
    S4Vectors:::recycleVector(value, index_len)
}

adjust_right_array_dim <- function(right_array, selection_dim)
{
    right_dim <- unname(dim(right_array))
    if (identical(selection_dim, right_dim))
        return(right_array)
    effdim_idx1 <- which(selection_dim != 1L)
    effdim_idx2 <- which(right_dim != 1L)
    if (!identical(selection_dim[effdim_idx1], right_dim[effdim_idx2]))
        stop(wmsg("dimensions of right array don't ",
                  "match dimensions of array selection"))
    dim(right_array) <- selection_dim
    right_array
}

### Same as 'array(data, selection_dim)' but:
### - returns an error if 'data' is longer than the array to construct
###   (strangely array() truncates the data in this case);
### - issues a warning if the length of the array to construct (which is
###   the length of the array selection in the context where array2()
###   is used) is not a multiple of 'length(data)'.
### Note that the error and warning messages are intentionally worded to
### make the most sense in the context where array2() is used.
array2 <- function(data, selection_dim)
{
    stopifnot(is.vector(data), is.integer(selection_dim))
    if (length(data) > prod(selection_dim))
        stop(wmsg("right value is longer than array selection"))
    ans <- array(vector(typeof(data), 1L), dim=selection_dim)
    ## Will issue "number of items to replace is not a multiple of
    ## replacement length" warning if 'prod(selection_dim)' is not a
    ## multiple of 'length(data)'.
    ans[] <- data
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subassign_Array_by_Lindex() and subassign_Array_by_Mindex() methods for
### SVT_SparseArray objects
###

.subassign_SVT_by_Lindex <- function(x, Lindex, value)
{
    x <- adjust_left_type(x, value)
    stopifnot(is.vector(Lindex), is.numeric(Lindex))

    ## No-op (except for type adjustment above) if array selection is empty.
    if (length(Lindex) == 0L)
        return(x)

    value <- normalize_right_value(value, type(x), length(Lindex))
    new_SVT <- SparseArray.Call("C_subassign_SVT_by_Lindex",
                                x@dim, x@type, x@SVT, FALSE, Lindex, value)
    BiocGenerics:::replaceSlots(x, SVT=new_SVT, check=FALSE)
}

setMethod("subassign_Array_by_Lindex", "SVT_SparseArray",
    .subassign_SVT_by_Lindex
)

.subassign_SVT_by_Mindex <- function(x, Mindex, value)
{
    x <- adjust_left_type(x, value)
    stopifnot(is.matrix(Mindex), is.numeric(Mindex))

    ## No-op (except for type adjustment above) if array selection is empty.
    if (nrow(Mindex) == 0L)
        return(x)

    value <- normalize_right_value(value, type(x), nrow(Mindex))
    new_SVT <- SparseArray.Call("C_subassign_SVT_by_Mindex",
                                x@dim, x@type, x@SVT, FALSE, Mindex, value)
    BiocGenerics:::replaceSlots(x, SVT=new_SVT, check=FALSE)
}

setMethod("subassign_Array_by_Mindex", "SVT_SparseArray",
    .subassign_SVT_by_Mindex
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subassign_Array_by_Nindex() method for SVT_SparseArray objects
###
### Like the 'index' argument in 'extract_array()', the 'Nindex' argument in
### all the functions below must be a **normalized** N-index, that is, a list
### with one list element per dimension in 'x' where each list element is
### either NULL or an integer vector of valid indices along the corresponding
### dimension in 'x'.

### 'Rvector' is considered "short" if it can be cleanly recycled along the
### first (a.k.a. leftmost or innermost) dimension of the array selection.
### This is a requirement of .subassign_SVT_with_short_Rvector().
### Note that 'Rvector' is never actually recycled. Instead the C code
### behind .subassign_SVT_with_short_Rvector() will cycle over its elements.
### 'selection_dim' is assumed to hold the dimensions of a non-empty
### array selection (in other words it cannot contain zeros).
### Returns TRUE or FALSE indicating whether 'Rvector' is considered "short"
### or not.
is_short <- function(Rvector_len, selection_dim)
{
    stopifnot(isSingleInteger(Rvector_len), is.integer(selection_dim))
    if (Rvector_len == 0L)
        stop(wmsg("right value has length zero"))
    selection_dim[[1L]] %% Rvector_len == 0L
}

### This handles subassignment by an N-index and with a right value that is
### a "short vector" that gets recycled along the first (a.k.a. leftmost
### or innermost) dimension of 'x', like in:
###     x[ , 1:2] <- 0
### or:
###     x[1:12, ] <- c(0.6, 0, 2.5)
### We want to support this in the most efficient way possible so we
### don't actually recycle the right value at the R level. Instead we
### will **virtually** recycle it at the C level.
### See is_short() above for more information.
.subassign_SVT_with_short_Rvector <- function(x, Nindex, Rvector)
{
    x <- adjust_left_type(x, Rvector)

    ## No-op (except for type change above) if array selection is empty.
    selection_dim <- S4Arrays:::get_Nindex_lengths(Nindex, x@dim)
    if (any(selection_dim == 0L))
        return(x)

    Rvector_len <- length(Rvector)
    stopifnot(is_short(Rvector_len, selection_dim))

    ## Prepare 'Noffs' and 'Rvector'.
    Norder <- S4Arrays:::get_Nindex_order(Nindex)
    Nindex <- S4Arrays:::subset_Nindex_by_Nindex(Nindex, Norder)
    Noffs <- Nindex2Noffs(Nindex)
    Norder1 <- Norder[[1L]]
    if (!is.null(Norder1))
        Rvector <- Rvector[((Norder1 - 1L) %% Rvector_len) + 1L]
    storage.mode(Rvector) <- type(x)

    new_SVT <- SparseArray.Call("C_subassign_SVT_with_short_Rvector",
                                x@dim, x@type, x@SVT, FALSE, Noffs, Rvector)
    BiocGenerics:::replaceSlots(x, SVT=new_SVT, check=FALSE)
}

.subassign_SVT_by_Noffs_with_Rarray <- function(x, Noffs, Rarray)
{
    new_SVT <- SparseArray.Call("C_subassign_SVT_with_Rarray",
                                x@dim, x@type, x@SVT, FALSE, Noffs, Rarray)
    BiocGenerics:::replaceSlots(x, SVT=new_SVT, check=FALSE)
}

.subassign_SVT_with_Rarray <- function(x, Nindex, Rarray)
{
    stopifnot(is(x, "SVT_SparseArray"), is.list(Nindex))
    check_svt_version(x)
    stopifnot(is.array(Rarray))

    ## Change 'x' type if necessary.
    new_type <- type(c(vector(type(x)), vector(type(Rarray))))
    type(x) <- new_type

    ## No-op (except for type change above) if array selection is empty.
    selection_dim <- S4Arrays:::get_Nindex_lengths(Nindex, x@dim)
    Rarray <- adjust_right_array_dim(Rarray, selection_dim)
    if (any(selection_dim == 0L))
        return(x)

    ## Prepare 'Noffs' and 'Rarray'.
    Norder <- S4Arrays:::get_Nindex_order(Nindex)
    Nindex <- S4Arrays:::subset_Nindex_by_Nindex(Nindex, Norder)
    Noffs <- Nindex2Noffs(Nindex)
    Rarray <- S4Arrays:::subset_by_Nindex(Rarray, Norder)
    storage.mode(Rarray) <- new_type

    .subassign_SVT_by_Noffs_with_Rarray(x, Noffs, Rarray)
}

.subassign_SVT_with_SVT <- function(x, Nindex, y)
{
    stopifnot(is(x, "SVT_SparseArray"), is.list(Nindex))
    check_svt_version(x)
    stopifnot(is(y, "SVT_SparseArray"))
    check_svt_version(y)

    ## Change 'x' type if necessary.
    new_type <- type(c(vector(type(x)), vector(type(y))))
    type(x) <- new_type

    ## No-op (except for type change above) if array selection is empty.
    selection_dim <- S4Arrays:::get_Nindex_lengths(Nindex, x@dim)
    y <- adjust_right_array_dim(y, selection_dim)
    if (any(selection_dim == 0L))
        return(x)

    ## Prepare 'Noffs' and 'y'.
    Norder <- S4Arrays:::get_Nindex_order(Nindex)
    Nindex <- S4Arrays:::subset_Nindex_by_Nindex(Nindex, Norder)
    Noffs <- Nindex2Noffs(Nindex)
    y <- S4Arrays:::subset_by_Nindex(y, Norder)
    type(y) <- new_type

    new_SVT <- SparseArray.Call("C_subassign_SVT_with_SVT",
                                x@dim, x@type, x@SVT, FALSE, Noffs,
                                y@dim, y@type, y@SVT, FALSE)
    BiocGenerics:::replaceSlots(x, SVT=new_SVT, check=FALSE)
}

.subassign_SVT_by_Nindex <- function(x, Nindex, value)
{
    stopifnot(is(x, "SVT_SparseArray"), is.list(Nindex))
    check_svt_version(x)
    if (is.vector(value)) {
        ## Change 'x' type if necessary.
        new_type <- type(c(vector(type(x)), vector(type(value))))
        type(x) <- new_type

        ## No-op (except for type change above) if array selection is empty.
        selection_dim <- S4Arrays:::get_Nindex_lengths(Nindex, x@dim)
        if (any(selection_dim == 0L))
            return(x)

        if (is_short(length(value), selection_dim))
            return(.subassign_SVT_with_short_Rvector(x, Nindex, value))

        ## Turn 'value' into an ordinary array of same dimensions as
        ## the array selection, with recycling if necessary.
        value <- array2(value, selection_dim)
    }
    if (is.array(value))
        return(.subassign_SVT_with_Rarray(x, Nindex, value))
    if (is(value, "SVT_SparseArray"))
        return(.subassign_SVT_with_SVT(x, Nindex, value))
    stop(wmsg("the right value must be an ordinary vector or array, ",
              "or an SVT_SparseArray object, for this subassignment"))
}

setMethod("subassign_Array_by_Nindex", "SVT_SparseArray",
    .subassign_SVT_by_Nindex
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### An alternate implementation of coercion from ordinary array to
### SVT_SparseArray that uses subassignment
###

### Not used at the moment.
### The workhorse behind coercion from array to SVT_SparseArray is
### .build_SVT_SparseArray_from_array() defined in R/SVT_SparseArray-class.R.
### build_SVT_SparseArray_from_array2() is an alternative to that.
### It starts by creating an all-zero SVT_SparseArray object 'svt0' of
### the same dimension as 'x' (the array to coerce), then it does:
###
###     svt0[ , , ] <- x  # 3D case
###
### Note that an important feature of .build_SVT_SparseArray_from_array()
### is its ability to delay the switch to the requested 'type' as much as
### possible i.e. until the leaves of the SVT_SparseArray object to return
### get created at the C level. build_SVT_SparseArray_from_array2() doesn't
### do that: it switches the type of the returned SVT_SparseArray, which is
### not as memory efficient.
build_SVT_SparseArray_from_array2 <- function(x, dimnames=NULL, type=NA)
{
    stopifnot(is.array(x))
    if (is.null(dimnames)) {
        ans_dimnames <- dimnames(x)
    } else {
        ans_dimnames <- S4Arrays:::normarg_dimnames(dimnames, dim(x))
    }
    svt0 <- SVT_SparseArray(dim=dim(x), dimnames=dimnames(x), type=type(x))
    Noffs <- vector(mode="list", length=length(dim(svt0)))
    ans <- .subassign_SVT_by_Noffs_with_Rarray(svt0, Noffs, x)
    if (!identical(type, NA))
        type(ans) <- type
    ans
}

