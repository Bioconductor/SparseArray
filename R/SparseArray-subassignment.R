### =========================================================================
### SparseArray subassignment
### -------------------------------------------------------------------------
###


.subassign_SVT_SparseArray_by_logical_array <- function(x, y, value)
    stop("subassignment operation not supported yet")

.subassign_SVT_SparseArray_by_Mindex <- function(x, Mindex, value)
{
    stopifnot(is(x, "SVT_SparseArray"),
              is.matrix(Mindex), is.numeric(Mindex))
    if (!is.vector(value))
        stop(wmsg("the supplied value must be a vector for this form ",
                  "of subassignment to an SVT_SparseArray object"))

    ## Change 'x' type if necessary. */
    new_type <- type(c(vector(type(x)), vector(type(value))))
    type(x) <- new_type

    ## No-op (except for type change above) if selection if empty. */
    if (nrow(Mindex) == 0L)
        return(x)

    ## Normalize 'value'. */
    if (length(value) == 0L)
        stop(wmsg("replacement has length zero"))
    storage.mode(value) <- new_type
    value <- S4Vectors:::recycleVector(value, nrow(Mindex))

    if (storage.mode(Mindex) != "integer")
        storage.mode(Mindex) <- "integer"
    new_SVT <- .Call2("C_subassign_SVT_by_Mindex",
                      x@dim, x@type, x@SVT, Mindex, value,
                      PACKAGE="SparseArray")
    BiocGenerics:::replaceSlots(x, SVT=new_SVT, check=FALSE)
}

.subassign_SVT_SparseArray_by_Lindex <- function(x, Lindex, value)
{
    stopifnot(is(x, "SVT_SparseArray"),
              is.vector(Lindex), is.numeric(Lindex))
    if (!is.vector(value))
        stop(wmsg("the supplied value must be a vector for this form ",
                  "of subassignment to an SVT_SparseArray object"))

    ## Change 'x' type if necessary. */
    new_type <- type(c(vector(type(x)), vector(type(value))))
    type(x) <- new_type

    ## No-op (except for type change above) if selection if empty. */
    if (length(Lindex) == 0L)
        return(x)

    ## Normalize 'value'. */
    if (length(value) == 0L)
        stop(wmsg("replacement has length zero"))
    storage.mode(value) <- new_type
    value <- S4Vectors:::recycleVector(value, length(Lindex))

    new_SVT <- .Call2("C_subassign_SVT_by_Lindex",
                      x@dim, x@type, x@SVT, Lindex, value,
                      PACKAGE="SparseArray")
    BiocGenerics:::replaceSlots(x, SVT=new_SVT, check=FALSE)
}


.subassign_SVT_with_short_Rvector <- function(x, index, Rvector)
{
    stopifnot(is.vector(Rvector))
    .Call2("C_subassign_SVT_with_short_Rvector",
           x@dim, x@type, x@SVT, index, Rvector, PACKAGE="SparseArray")
}

.subassign_SVT_with_Rarray <- function(x, index, Rarray)
{
    stopifnot(is.array(Rarray))
    .Call2("C_subassign_SVT_with_Rarray", x@dim, x@type, x@SVT,
           index, Rarray, PACKAGE="SparseArray")
}

.subassign_SVT_with_SVT <- function(x, index, v)
{
    stopifnot(is(v, "SVT_SparseArray"))
    .Call2("C_subassign_SVT_with_SVT",
           x@dim, x@type, x@SVT, index, v@dim, v@type, v@SVT,
           PACKAGE="SparseArray")
}

### Like for 'extract_array()', the supplied 'index' must be a list with
### one list element per dimension in 'x'. Each list element must be an
### integer vector of valid indices along the corresponding dimension
### in 'x', or a NULL.
.subassign_SVT_SparseArray_by_index <- function(x, index, value)
{
    stopifnot(is(x, "SVT_SparseArray"), is.list(index))
    if (!is.vector(value) && !is.array(value) && !is(value, "SVT_SparseArray"))
        stop(wmsg("the supplied value must be an ordinary vector or array, ",
                  "or an SVT_SparseArray object, for this subassignment"))

    ## Change 'x' type if necessary. */
    new_type <- type(c(vector(type(x)), vector(type(value))))
    type(x) <- new_type

    ## No-op (except for type change above) if selection if empty. */
    selection_dim <- S4Arrays:::get_Nindex_lengths(index, x@dim)
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
            new_SVT <- .subassign_SVT_with_short_Rvector(x, index, value)
        } else {
            ## Turn 'value' into an ordinary array of the same dimensions
            ## as the selection, with recycling if necessary.
            a <- array(vector(typeof(value), 1L), dim=selection_dim)
            a[] <- value
            new_SVT <- .subassign_SVT_with_Rarray(x, index, value)
        }
    } else {
        if (!identical(selection_dim, unname(dim(value))))
            stop(wmsg("the selection and supplied value must have ",
                      "the same dimensions"))
        if (is.array(value)) {
            storage.mode(value) <- new_type
            new_SVT <- .subassign_SVT_with_Rarray(x, index, value)
        } else {
            type(value) <- new_type
            new_SVT <- .subassign_SVT_with_SVT(x, index, value)
        }
    }
    BiocGenerics:::replaceSlots(x, SVT=new_SVT, check=FALSE)
}

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
            return(.subassign_SVT_SparseArray_by_logical_array(x, i, value))
        if (is.matrix(i) && is.numeric(i))
            return(.subassign_SVT_SparseArray_by_Mindex(x, i, value))
        ## Linear single bracket subassignment e.g. x[5:2] <- 4.
        return(.subassign_SVT_SparseArray_by_Lindex(x, i, value))
    }
    if (nsubscript != x_ndim)
        stop(wmsg("incorrect number of subscripts"))
    index <- S4Arrays:::normalize_Nindex(Nindex, x)
    .subassign_SVT_SparseArray_by_index(x, index, value)
}

setReplaceMethod("[", "SVT_SparseArray", .subassign_SVT_SparseArray)

