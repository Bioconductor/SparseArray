### =========================================================================
### SparseArray subassignment
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .subassign_SVT_by_logical_array()
###

.subassign_SVT_by_logical_array <- function(x, y, value)
    stop("subassignment operation not supported yet")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .subassign_SVT_by_Mindex()
### .subassign_SVT_by_Lindex()
###

.adjust_left_type <- function(x, value)
{
    stopifnot(is(x, "SVT_SparseArray"))
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
.normalize_right_value <- function(value, left_type, index_len)
{
    if (length(value) == 0L)
        stop(wmsg("replacement has length zero"))
    storage.mode(value) <- left_type
    S4Vectors:::recycleVector(value, index_len)
}

.subassign_SVT_by_Mindex <- function(x, Mindex, value)
{
    x <- .adjust_left_type(x, value)
    stopifnot(is.matrix(Mindex), is.numeric(Mindex))

    ## No-op (except for type adjustment above) if selection if empty.
    if (nrow(Mindex) == 0L)
        return(x)

    value <- .normalize_right_value(value, type(x), nrow(Mindex))

    if (storage.mode(Mindex) != "integer")
        storage.mode(Mindex) <- "integer"
    new_SVT <- SparseArray.Call("C_subassign_SVT_by_Mindex",
                                x@dim, x@type, x@SVT, Mindex, value)
    BiocGenerics:::replaceSlots(x, SVT=new_SVT, check=FALSE)
}

.subassign_SVT_by_Lindex <- function(x, Lindex, value)
{
    x <- .adjust_left_type(x, value)
    stopifnot(is.vector(Lindex), is.numeric(Lindex))

    ## No-op (except for type adjustment above) if selection if empty.
    if (length(Lindex) == 0L)
        return(x)

    value <- .normalize_right_value(value, type(x), length(Lindex))

    new_SVT <- SparseArray.Call("C_subassign_SVT_by_Lindex",
                                x@dim, x@type, x@SVT, Lindex, value)
    BiocGenerics:::replaceSlots(x, SVT=new_SVT, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .subassign_SVT_by_Nindex()
###
### Like the 'index' argument in 'extract_array()', the 'Nindex' argument in
### all the functions below must be an N-index, that is, a list with one list
### element per dimension in 'x'. Each list element must be an integer vector
### of valid indices along the corresponding dimension in 'x', or a NULL.

.subassign_SVT_with_short_Rvector <- function(x, Nindex, Rvector)
{
    stopifnot(is.vector(Rvector))
    SparseArray.Call("C_subassign_SVT_with_short_Rvector",
                     x@dim, x@type, x@SVT, Nindex, Rvector)
}

.subassign_SVT_with_Rarray <- function(x, Nindex, Rarray)
{
    stopifnot(is.array(Rarray))
    SparseArray.Call("C_subassign_SVT_with_Rarray",
                     x@dim, x@type, x@SVT, Nindex, Rarray)
}

.subassign_SVT_with_SVT <- function(x, Nindex, v)
{
    stopifnot(is(v, "SVT_SparseArray"))
    check_svt_version(v)
    SparseArray.Call("C_subassign_SVT_with_SVT",
                     x@dim, x@type, x@SVT, Nindex, v@dim, v@type, v@SVT)
}

.subassign_SVT_by_Nindex <- function(x, Nindex, value)
{
    stopifnot(is(x, "SVT_SparseArray"), is.list(Nindex))
    check_svt_version(x)
    if (!is.vector(value) && !is.array(value) && !is(value, "SVT_SparseArray"))
        stop(wmsg("the supplied value must be an ordinary vector or array, ",
                  "or an SVT_SparseArray object, for this subassignment"))

    ## Change 'x' type if necessary. */
    new_type <- type(c(vector(type(x)), vector(type(value))))
    type(x) <- new_type

    ## No-op (except for type change above) if selection if empty. */
    selection_dim <- S4Arrays:::get_Nindex_lengths(Nindex, x@dim)
    if (any(selection_dim == 0L))
        return(x)

    if (is.vector(value)) {
        value_len <- length(value)
        if (value_len == 0L)
            stop(wmsg("replacement has length zero"))
        selection_len <- prod(selection_dim)
        if (value_len > selection_len)
            stop(wmsg("the supplied value is longer than the selection"))
        storage.mode(value) <- new_type
        if (value_len <= selection_dim[[1L]] &&
            selection_dim[[1L]] %% value_len == 0L)
        {
            ## We want to support things like 'x[ , 1:2] <- 0'
            ## or 'x[1:12, ] <- c(0.6, 0, 2.5)' in the most efficient
            ## way so no recycling of 'value' at the R level.
            new_SVT <- .subassign_SVT_with_short_Rvector(x, Nindex, value)
        } else {
            ## Turn 'value' into an ordinary array of the same dimensions
            ## as the selection, with recycling if necessary.
            a <- array(vector(typeof(value), 1L), dim=selection_dim)
            a[] <- value
            new_SVT <- .subassign_SVT_with_Rarray(x, Nindex, value)
        }
    } else {
        if (!identical(selection_dim, unname(dim(value))))
            stop(wmsg("the selection and supplied value must have ",
                      "the same dimensions"))
        if (is.array(value)) {
            storage.mode(value) <- new_type
            new_SVT <- .subassign_SVT_with_Rarray(x, Nindex, value)
        } else {
            type(value) <- new_type
            new_SVT <- .subassign_SVT_with_SVT(x, Nindex, value)
        }
    }
    BiocGenerics:::replaceSlots(x, SVT=new_SVT, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subassignment method (`[<-`) for SVT_SparseArray objects
###

.subassign_SVT_SparseArray <- function(x, i, j, ..., value)
{
    if (missing(x))
        stop(wmsg("'x' is missing"))
    Nindex <- S4Arrays:::extract_Nindex_from_syscall(sys.call(), parent.frame())
    nsubscript <- length(Nindex)
    x_dim <- dim(x)
    x_ndim <- length(x_dim)
    if (nsubscript == 1L) {
        i <- Nindex[[1L]]
        if (type(i) == "logical" && identical(x_dim, dim(i)))
            return(.subassign_SVT_by_logical_array(x, i, value))
        if (is.matrix(i) && is.numeric(i))
            return(.subassign_SVT_by_Mindex(x, i, value))
        ## Linear single bracket subassignment e.g. x[5:2] <- 4.
        return(.subassign_SVT_by_Lindex(x, i, value))
    }
    if (nsubscript != x_ndim)
        stop(wmsg("incorrect number of subscripts"))
    Nindex <- S4Arrays:::normalize_Nindex(Nindex, x)
    .subassign_SVT_by_Nindex(x, Nindex, value)
}

setReplaceMethod("[", "SVT_SparseArray", .subassign_SVT_SparseArray)

