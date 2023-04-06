### =========================================================================
### Combining SparseArray objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Binding()
###
### We only support binding SparseArray objects along the rows or the cols
### at the moment. No binding along an arbitrary dimension yet! (i.e. no
### "abind" method yet)
###

### Similar to simple_abind() (see abind.R) but works on COO_SparseArray
### objects.
.abind_COO_SparseArray_objects <- function(objects, along)
{
    stopifnot(is.list(objects))
    if (length(objects) == 0L)
        return(NULL)

    ## Check dim compatibility.
    dims <- S4Arrays:::get_dims_to_bind(objects, along)
    if (is.character(dims))
        stop(wmsg(dims))
    if (length(objects) == 1L)
        return(objects[[1L]])

    ## Compute 'ans_dim' and 'ans_dimnames'.
    ans_dim <- S4Arrays:::combine_dims_along(dims, along)
    ans_dimnames <- S4Arrays:::combine_dimnames_along(objects, dims, along)

    ## Combine the "nzcoo" slots.
    offsets <- cumsum(dims[along, -ncol(dims)])
    nzcoo_list <- lapply(seq_along(objects),
        function(i) {
            object <- objects[[i]]
            nzcoo <- object@nzcoo
            if (i >= 2L)
                nzcoo[ , along] <- nzcoo[ , along, drop=FALSE] +
                                   offsets[[i - 1L]]
            nzcoo
        }
    )
    ans_nzcoo <- do.call(rbind, nzcoo_list)

    ## Combine the "nzvals" slots.
    ans_nzvals <- unlist(lapply(objects, slot, "nzvals"), use.names=FALSE)

    COO_SparseArray(ans_dim, ans_nzcoo, ans_nzvals, ans_dimnames, check=FALSE)
}

setMethod("arbind", "COO_SparseArray",
    function(...) .abind_COO_SparseArray_objects(list(...), along=1L)
)

setMethod("acbind", "COO_SparseArray",
    function(...) .abind_COO_SparseArray_objects(list(...), along=2L)
)

### Similar to simple_abind() (see abind.R) but works on SVT_SparseArray
### objects.
.abind_SVT_SparseArray_objects <- function(objects, along)
{
    stopifnot(is.list(objects))
    if (length(objects) == 0L)
        return(NULL)

    ## Check dim compatibility.
    dims <- S4Arrays:::get_dims_to_bind(objects, along)
    if (is.character(dims))
        stop(wmsg(dims))
    if (length(objects) == 1L)
        return(objects[[1L]])

    ## Compute 'ans_dimnames'.
    ans_dimnames <- S4Arrays:::combine_dimnames_along(objects, dims, along)

    ## Compute 'ans_type'.
    ans_type <- type(unlist(
        lapply(objects, function(object) vector(type(object)))
    ))
    objects <- lapply(objects, `type<-`, ans_type)

    ## Returns 'ans_dim' and 'ans_SVT' in a list of length 2.
    C_ans <- .Call2("C_abind_SVT_SparseArray_objects",
                    objects, along, ans_type, PACKAGE="SparseArray")
    ans_dim <- C_ans[[1L]]
    ans_SVT <- C_ans[[2L]]

    new_SVT_SparseArray(ans_dim, ans_dimnames, ans_type, ans_SVT, check=FALSE)
}

setMethod("arbind", "SVT_SparseArray",
    function(...) .abind_SVT_SparseArray_objects(list(...), along=1L)
)

setMethod("acbind", "SVT_SparseArray",
    function(...) .abind_SVT_SparseArray_objects(list(...), along=2L)
)

### TODO: The methods below are defined for SparseArray objects but it seems
### that they could as well be defined more generally for Array objects.

### The generics have the 'deparse.level' argument. We ignore it.
setMethod("rbind", "SparseArray", function(...) arbind(...))
setMethod("cbind", "SparseArray", function(...) acbind(...))

### Arguments 'use.names', 'ignore.mcols', and 'check' are ignored.
setMethod("bindROWS", "SparseArray",
    function(x, objects=list(), use.names=TRUE, ignore.mcols=FALSE, check=TRUE)
    {
        args <- c(list(x), unname(objects))
        do.call(rbind, args)
    }
)

