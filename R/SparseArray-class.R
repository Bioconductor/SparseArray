### =========================================================================
### SparseArray objects
### -------------------------------------------------------------------------


### SparseArray is virtual class with 2 concrete subclasses: COO_SparseArray
### and SVT_SparseArray.
###
### The SparseArray API:
### 1) Implemented in this file:
###    - Getters dim(), length(), dimnames(), type()
###    - Setters `dimnames<-`() and `type<-`()
###    - An is_sparse() method that returns TRUE
###    - nzcount(), nzwhich(), and nzvals() generics
###    - sparsity()
### 2) Implemented elsewhere:
###    - nzcount(), nzwhich(), and nzvals() methods
###    - as.array()
###    - extract_array() and extract_sparse_array()
###    - Subsetting (`[`) and subassignment (`[<-`)
###    - read_block_as_dense() and read_block_as_sparse()
###    - abind(), arbind(), acbind()
###    - aperm()

setClass("SparseArray",
    contains="Array",
    representation(
        "VIRTUAL",
        dim="integer",
        dimnames="list"    # List with one list element per dimension. Each
                           # list element must be NULL or a character vector.
    ),
    prototype(
        dim=0L,
        dimnames=list(NULL)
    )
)

.validate_SparseArray <- function(x)
{
    msg <- S4Arrays:::validate_dim_slot(x, "dim")
    if (!isTRUE(msg))
        return(msg)
    msg <- S4Arrays:::validate_dimnames_slot(x, x@dim)
    if (!isTRUE(msg))
        return(msg)
    TRUE
}
setValidity2("SparseArray", .validate_SparseArray)

### Extending RectangularData gives us a few things for free (e.g. validity
### method for RectangularData objects, head(), tail(), etc...). Note
### that even though SparseMatrix already extends Array (via SparseArray),
### we need to make it a *direct* child of Array, and to list Array *before*
### RectangularData in the 'contains' field below. This will ensure that
### method dispatch will always choose the method for Array in case a generic
### has methods defined for both, Array and RectangularData.
### Note that the fact that we need this "hack" is a hint that we could
### achieve a cleaner class hierarchy by inserting a Matrix class in it.
### Matrix would contain Array and RectangularData (in that order). Then
### SparseMatrix would contain SparseArray and Matrix (in that order).
### Unfortunately the Matrix package already defines a Matrix class so
### we would need to use a different name.
setClass("SparseMatrix",
    contains=c("SparseArray", "Array", "RectangularData"),
    representation("VIRTUAL"),
    prototype(
        dim=c(0L, 0L),
        dimnames=list(NULL, NULL)
    )
)

.validate_SparseMatrix <- function(x)
{
    if (length(x@dim) != 2L)
        return("'dim' slot must be an integer vector of length 2")
    TRUE
}
setValidity2("SparseMatrix", .validate_SparseMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dim(), dimnames(), and `dimnames<-`()
###

setMethod("dim", "SparseArray", function(x) x@dim)

setMethod("dimnames", "SparseArray",
    function(x) S4Arrays:::simplify_NULL_dimnames(x@dimnames)
)

setReplaceMethod("dimnames", "SparseArray",
    function(x, value)
    {
        x@dimnames <- S4Arrays:::normarg_dimnames(value, dim(x))
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### coercion_can_introduce_zeros()
###

coercion_can_introduce_zeros <- function(from_type, to_type)
{
    if (!isSingleString(from_type))
        stop(wmsg("'from_type' must be a single string"))
    if (!isSingleString(to_type))
        stop(wmsg("'to_type' must be a single string"))
    if (!(to_type %in% c("double", "logical")))
        stop(wmsg("'to_type' must be \"double\" or \"logical\""))
    SparseArray.Call("C_coercion_can_introduce_zeros", from_type, to_type)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### is_sparse() method
###

setMethod("is_sparse", "SparseArray", function(x) TRUE)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### nzcount(), nzwhich(), and nzvals() generics + sparsity()
###

### Returns the number of nonzero array elements in 'x'.
setGeneric("nzcount", function(x) standardGeneric("nzcount"))

### Not 100% reliable because [d|l]gCMatrix objects are allowed to have
### zeros in their @x slot! See src/SVT_SparseArray_class.c for an example.
setMethod("nzcount", "CsparseMatrix", function(x) length(x@i))
setMethod("nzcount", "RsparseMatrix", function(x) length(x@j))

### Returns the indices of the nonzero array elements in 'x', either as
### an L-index (if 'arr.ind=FALSE') or as an M-index (if 'arr.ind=TRUE').
setGeneric("nzwhich", signature="x",
    function(x, arr.ind=FALSE) standardGeneric("nzwhich")
)

default_nzwhich <- function(x, arr.ind=FALSE)
{
    if (!isTRUEorFALSE(arr.ind))
        stop(wmsg("'arr.ind' must be TRUE or FALSE"))
    ## Make sure to use 'type()' and not 'typeof()'.
    zero <- vector(type(x), length=1L)
    is_nonzero <- x != zero
    which(is_nonzero | is.na(is_nonzero), arr.ind=arr.ind, useNames=FALSE)
}
setMethod("nzwhich", "ANY", default_nzwhich)

### default_nzwhich() above works on a CsparseMatrix derivative but
### nzwhich_CsparseMatrix() will typically be 50x or 100x faster, or more!
### However, this is **NOT** 100% reliable because [d|l]gCMatrix objects are
### allowed to have zeros in their @x slot! See src/SVT_SparseArray_class.c
### for an example.
nzwhich_CsparseMatrix <- function(x, arr.ind=FALSE)
{
    if (!isTRUEorFALSE(arr.ind))
        stop(wmsg("'arr.ind' must be TRUE or FALSE"))
    x_nrow <- nrow(x)
    x_ncol <- ncol(x)
    offsets <- rep.int(x_nrow * seq_len(x_ncol) - (x_nrow - 1L), diff(x@p))
    ans <- x@i + offsets
    if (!arr.ind)
        return(ans)
    Lindex2Mindex(ans, dim(x))
}
setMethod("nzwhich", "CsparseMatrix", nzwhich_CsparseMatrix)

### Returns the values of the nonzero array elements in a vector of the
### same type() as 'x' and parallel to nzwhich(x).
### Equivalent to x[nzwhich(x)] (and that's what the default method below
### does). However specialized methods can make this dramatically faster.
setGeneric("nzvals", function(x) standardGeneric("nzvals"))

setMethod("nzvals", "ANY", function(x) x[nzwhich(x)])

### Not 100% reliable because [d|l]gCMatrix objects are allowed to have
### zeros in their @x slot! See src/SVT_SparseArray_class.c for an example.
setMethod("nzvals", "dgCMatrix", function(x) x@x)
setMethod("nzvals", "lgCMatrix", function(x) x@x)

sparsity <- function(x) { 1 - nzcount(x) / length(x) }


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show()
###

setMethod("classNameForDisplay", "SparseArray", function(x) "SparseArray")
setMethod("classNameForDisplay", "SparseMatrix", function(x) "SparseMatrix")

.show_headline_part1 <- function(x)
{
    sprintf("<%s %s> of type \"%s\" ", paste0(dim(x), collapse=" x "),
                                       classNameForDisplay(x), type(x))
}

.show_headline_part2 <- function(x)
{
    ## Calling nzcount(x) will fail if 'x' is an SVT_SparseArray object
    ## that uses version 0 of the SVT internal layout.
    x_nzcount <- nzcount(x)
    x_density <- x_nzcount / length(x)
    sprintf("[nzcount=%s (%s%%)]", format(x_nzcount),
                                   signif(100 * x_density, digits=2))
}

setMethod("show", "SparseArray",
    function(object)
    {
        ## Only reason we print the headline in 2 steps is because we
        ## want to make sure to print at least something (part1) even
        ## when printing part2 is going to fail. This will happen for
        ## example if the call to nzcount() in .show_headline_part2() fails.
        cat(.show_headline_part1(object))
        cat(.show_headline_part2(object))
        if (any(dim(object) == 0L)) {
            cat("\n")
            return()
        }
        cat(":\n", sep="")
        S4Arrays:::print_some_array_elements(object)
    }
)

