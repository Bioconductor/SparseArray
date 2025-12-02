### =========================================================================
### NaArray subassignment
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subassign_Array_by_Lindex() and subassign_Array_by_Mindex() methods for
### NaArray objects
###

.subassign_NaSVT_by_Lindex <- function(x, Lindex, value)
{
    x <- adjust_left_type(x, value)
    stopifnot(is.vector(Lindex), is.numeric(Lindex))

    ## No-op (except for type adjustment above) if selection is empty.
    if (length(Lindex) == 0L)
        return(x)

    value <- normalize_right_value(value, type(x), length(Lindex))
    new_NaSVT <- SparseArray.Call("C_subassign_SVT_by_Lindex",
                                  x@dim, x@type, x@NaSVT, TRUE, Lindex, value)
    BiocGenerics:::replaceSlots(x, NaSVT=new_NaSVT, check=FALSE)
}

setMethod("subassign_Array_by_Lindex", "NaArray",
    function(x, Lindex, value) .subassign_NaSVT_by_Lindex(x, Lindex, value)
)

.subassign_NaSVT_by_Mindex <- function(x, Mindex, value)
{
    x <- adjust_left_type(x, value)
    stopifnot(is.matrix(Mindex), is.numeric(Mindex))

    ## No-op (except for type adjustment above) if array selection is empty.
    if (nrow(Mindex) == 0L)
        return(x)

    value <- normalize_right_value(value, type(x), nrow(Mindex))
    new_NaSVT <- SparseArray.Call("C_subassign_SVT_by_Mindex",
                                  x@dim, x@type, x@NaSVT, TRUE, Mindex, value)
    BiocGenerics:::replaceSlots(x, NaSVT=new_NaSVT, check=FALSE)
}

setMethod("subassign_Array_by_Mindex", "NaArray",
    .subassign_NaSVT_by_Mindex
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subassign_Array_by_Nindex() method for NaArray objects
###

.subassign_NaSVT_with_short_Rvector <- function(x, Nindex, Rvector)
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

    new_NaSVT <- SparseArray.Call("C_subassign_SVT_with_short_Rvector",
                                  x@dim, x@type, x@NaSVT, TRUE, Noffs, Rvector)
    BiocGenerics:::replaceSlots(x, NaSVT=new_NaSVT, check=FALSE)
}

.subassign_NaSVT_by_Noffs_with_Rarray <- function(x, Noffs, Rarray)
{
    new_NaSVT <- SparseArray.Call("C_subassign_SVT_with_Rarray",
                                  x@dim, x@type, x@NaSVT, TRUE, Noffs, Rarray)
    BiocGenerics:::replaceSlots(x, NaSVT=new_NaSVT, check=FALSE)
}

.subassign_NaSVT_with_Rarray <- function(x, Nindex, Rarray)
{
    stopifnot(is(x, "NaArray"), is.list(Nindex))
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

    .subassign_NaSVT_by_Noffs_with_Rarray(x, Noffs, Rarray)
}

.subassign_NaSVT_with_NaSVT <- function(x, Nindex, y)
{
    stopifnot(is(x, "NaArray"), is.list(Nindex))
    check_svt_version(x)
    stopifnot(is(y, "NaArray"))
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

    new_NaSVT <- SparseArray.Call("C_subassign_SVT_with_SVT",
                                  x@dim, x@type, x@NaSVT, TRUE, Noffs,
                                  y@dim, y@type, y@NaSVT, TRUE)
    BiocGenerics:::replaceSlots(x, NaSVT=new_NaSVT, check=FALSE)
}

.subassign_NaSVT_by_Nindex <- function(x, Nindex, value)
{
    stopifnot(is(x, "NaArray"), is.list(Nindex))
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
            return(.subassign_NaSVT_with_short_Rvector(x, Nindex, value))

        ## Turn 'value' into an ordinary array of same dimensions as
        ## the array selection, with recycling if necessary.
        value <- array2(value, selection_dim)
    }
    if (is.array(value))
        return(.subassign_NaSVT_with_Rarray(x, Nindex, value))
    if (is(value, "NaArray"))
        return(.subassign_NaSVT_with_NaSVT(x, Nindex, value))
    stop(wmsg("the right value must be an ordinary vector or array, ",
              "or an NaArray object, for this subassignment"))
}

setMethod("subassign_Array_by_Nindex", "NaArray",
    .subassign_NaSVT_by_Nindex
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### An alternate implementation of coercion from ordinary array to
### NaArray that uses subassignment
###
### Not used at the moment. See comment for build_SVT_SparseArray_from_array2()
### in R/SparseArray-subassignment.R
###

build_NaArray_from_array2 <- function(x, dimnames=NULL, type=NA)
{
    stopifnot(is.array(x))
    if (is.null(dimnames)) {
        ans_dimnames <- dimnames(x)
    } else {
        ans_dimnames <- S4Arrays:::normarg_dimnames(dimnames, dim(x))
    }
    naa0 <- NaArray(dim=dim(x), dimnames=dimnames(x), type=type(x))
    Noffs <- vector(mode="list", length=length(dim(naa0)))
    ans <- .subassign_NaSVT_by_Noffs_with_Rarray(naa0, Noffs, x)
    if (!identical(type, NA))
        type(ans) <- type
    ans
}

