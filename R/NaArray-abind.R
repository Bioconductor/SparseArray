### =========================================================================
### Combining multidimensional NaArray objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### abind()
###
### Note that making abind() work on NaArray objects also makes arbind(),
### acbind(), rbind(), cbind(), and bindROWS() work on these objects.
###

.abind_NaArray_objects <- function(..., along=NULL, rev.along=NULL)
{
    objects <- S4Arrays:::homogenize_abind_input_objects(list(...),
                                           along=along, rev.along=rev.along)
    if (length(objects) == 0L)
        return(NULL)

    x <- objects[[1L]]
    ## Note that unary abind() is not necessarily a no-op because
    ## homogenize_abind_input_objects() can add dimensions to the input object.
    if (length(objects) == 1L)
        return(x)

    ## Compute 'ans_dimnames'.
    dims <- attr(objects, "dims")
    along <- attr(objects, "along")
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

