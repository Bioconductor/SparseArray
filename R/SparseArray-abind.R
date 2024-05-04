### =========================================================================
### Combining multidimensional SparseArray objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### abind()
###

.abind_COO_SparseArray_objects <- function(objects, dims, along, ans_dimnames)
{
    ## Compute 'ans_dim'.
    ans_dim <- S4Arrays:::combine_dims_along(dims, along)

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

    ## Combine the @nzdata slots.
    ans_nzdata <- unlist(lapply(objects, slot, "nzdata"), use.names=FALSE)

    COO_SparseArray(ans_dim, ans_nzcoo, ans_nzdata, ans_dimnames, check=FALSE)
}

.abind_SVT_SparseArray_objects <- function(objects, along, ans_dimnames)
{
    ## Compute 'ans_type'.
    ans_type <- type(unlist(
        lapply(objects, function(object) vector(type(object)))
    ))
    objects <- lapply(objects, `type<-`, ans_type)

    ## Returns 'ans_dim' and 'ans_SVT' in a list of length 2.
    C_ans <- SparseArray.Call("C_abind_SVT_SparseArray_objects",
                              objects, along, ans_type)
    ans_dim <- C_ans[[1L]]
    ans_SVT <- C_ans[[2L]]

    new_SVT_SparseArray(ans_dim, ans_dimnames, ans_type, ans_SVT, check=FALSE)
}

.abind_SparseArray_objects <- function(..., along=NULL, rev.along=NULL)
{
    objects <- S4Vectors:::delete_NULLs(list(...))
    if (length(objects) == 0L)
        return(NULL)

    ndims <- vapply(objects, function(object) length(dim(object)), integer(1))
    N <- max(ndims)
    along <- S4Arrays:::get_along(N, along=along, rev.along=rev.along)
    ans_ndim <- max(N, along)
    objects <- S4Arrays:::add_missing_dims(objects, ans_ndim)

    ## Check dim compatibility.
    dims <- S4Arrays:::get_dims_to_bind(objects, along)
    if (is.character(dims))
        stop(wmsg(dims))
    x <- objects[[1L]]
    if (length(objects) == 1L)
        return(x)

    ## Compute 'ans_dimnames'.
    ans_dimnames <- S4Arrays:::combine_dimnames_along(objects, dims, along)

    if (is(x, "COO_SparseArray"))
        return(.abind_COO_SparseArray_objects(objects, dims, along,
                                              ans_dimnames))
    if (is(x, "SVT_SparseArray")) {
        check_svt_version(x)
        return(.abind_SVT_SparseArray_objects(objects, along,
                                              ans_dimnames))
    }
    stop(wmsg(class(x)[[1L]], " objects are not supported"))
}

setMethod("abind", "SparseArray", .abind_SparseArray_objects)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rbind(), cbind()
###

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

