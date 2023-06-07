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
###    - nzcount() generic
###    - sparsity()
### 2) Implemented elsewhere:
###    - nzcount() methods
###    - which() methods
###    - as.array()
###    - extract_array() and extract_sparse_array()
###    - Subsetting (`[`) and subassignment (`[<-`)
###    - read_block_as_dense() and read_block_as_sparse()
###    - arbind() and acbind()
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
### which_is_nonzero() and coercion_can_introduce_zeros()
###

which_is_nonzero <- function(x, arr.ind=FALSE)
{
    ## Make sure to use 'type()' and not 'typeof()'.
    zero <- vector(type(x), length=1L)
    is_nonzero <- x != zero
    which(is_nonzero | is.na(is_nonzero), arr.ind=arr.ind)
}

coercion_can_introduce_zeros <- function(from_type, to_type)
{
    if (!isSingleString(from_type))
        stop(wmsg("'from_type' must be a single string"))
    if (!isSingleString(to_type))
        stop(wmsg("'to_type' must be a single string"))
    if (!(to_type %in% c("double", "logical")))
        stop(wmsg("'to_type' must be \"double\" or \"logical\""))
    .Call2("C_coercion_can_introduce_zeros", from_type, to_type,
           PACKAGE="SparseArray")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### is_sparse(), nzcount(), and sparsity()
###

setMethod("is_sparse", "SparseArray", function(x) TRUE)

setGeneric("nzcount", function(x) standardGeneric("nzcount"))

### Not 100% reliable because [d|l]gCMatrix objects are allowed to have
### zeros in their "x" slot! See src/SVT_SparseArray_class.c for an example.
setMethod("nzcount", "CsparseMatrix", function(x) length(x@i))
setMethod("nzcount", "RsparseMatrix", function(x) length(x@j))

sparsity <- function(x) { 1 - nzcount(x) / length(x) }


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show()
###

SparseArray_as_one_line_summary <- function(x)
{
    sprintf("<%s %s> of type \"%s\" (nzcount=%s)",
            paste0(dim(x), collapse=" x "), class(x),
            type(x), format(nzcount(x)))
}

setMethod("show", "SparseArray",
    function(object)
    {
        #grey <- make_style("grey")
        #cat(grey(SparseArray_as_one_line_summary(object)))
        cat(SparseArray_as_one_line_summary(object))
        if (any(dim(object) == 0L)) {
            cat("\n")
            return()
        }
        #cat(grey(":"), "\n", sep="")
        cat(":\n", sep="")
        S4Arrays:::print_some_array_elements(object)
    }
)

