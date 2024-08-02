### =========================================================================
### Combining multidimensional NaArray objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### abind()
###

.abind_NaArray_objects <- function(..., along=NULL, rev.along=NULL)
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

    ## Compute 'ans_type'.
    ans_type <- type(unlist(
        lapply(objects, function(object) vector(type(object)))
    ))
    objects <- lapply(objects, `type<-`, ans_type)

    ## Returns 'ans_dim' and 'ans_NaSVT' in a list of length 2.
    C_ans <- SparseArray.Call("C_abind_SVT_SparseArray_objects",
                              objects, "NaSVT", along, ans_type)
    ans_dim <- C_ans[[1L]]
    ans_NaSVT <- C_ans[[2L]]

    new_NaArray(ans_dim, ans_dimnames, ans_type, ans_NaSVT, check=FALSE)
}

setMethod("abind", "NaArray", .abind_NaArray_objects)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rbind(), cbind()
###

### TODO: The methods below are defined for NaArray objects but it seems
### that they could as well be defined more generally for Array objects.

### The generics have the 'deparse.level' argument. We ignore it.
setMethod("rbind", "NaArray", function(...) arbind(...))
setMethod("cbind", "NaArray", function(...) acbind(...))

### Arguments 'use.names', 'ignore.mcols', and 'check' are ignored.
setMethod("bindROWS", "NaArray",
    function(x, objects=list(), use.names=TRUE, ignore.mcols=FALSE, check=TRUE)
    {
        args <- c(list(x), unname(objects))
        do.call(rbind, args)
    }
)

