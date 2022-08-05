### =========================================================================
### SparseArray objects
### -------------------------------------------------------------------------


### Virtual class with 2 concrete subclasses: COO_SparseArray and
### SVT_SparseArray.
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
### SparseArray API:
### - Getters: dim(), length(), dimnames(), type().
### - Setters: `dimnames<-`(), `type<-`().
### - is_sparse().

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
    if (identical(to_type, "double"))
        return(from_type %in% c("logical", "integer", "raw"))
    if (identical(to_type, "logical"))
        return(from_type %in% c("integer", "double", "complex", "raw"))
    stop(wmsg("'to_type' must be \"double\" or \"logical\""))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generics is_sparse(), is_sparse<-(), nzcount(), and sparsity()
###

### is_sparse() detects **structural** sparsity which is a global qualitative
### property of array-like object 'x' rather than a quantitative one.
### In other words it doesn't look at the data in 'x' to decide whether 'x'
### should be considered sparse or not. Said otherwise, it is NOT about
### quantitative sparsity measured by sparsity().
### IMPORTANT: Seeds for which is_sparse() returns TRUE **must** support
### extract_sparse_array(). More info about this in SparseArray-subsetting.R
### where the extract_sparse_array() generic is defined.
setGeneric("is_sparse", function(x) standardGeneric("is_sparse"))

setGeneric("is_sparse<-", signature="x",
    function(x, value) standardGeneric("is_sparse<-")
)

### By default, nothing is considered sparse.
setMethod("is_sparse", "ANY", function(x) FALSE)

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

