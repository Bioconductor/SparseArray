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
### 2) Implemented elsewhere:
###    - nzcount(), nzwhich(), nzvals(), and `nzvals<-`()
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
### is_sparse() method
###

setMethod("is_sparse", "SparseArray", function(x) TRUE)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show()
###

setMethod("classNameForDisplay", "SparseArray", function(x) "SparseArray")
setMethod("classNameForDisplay", "SparseMatrix", function(x) "SparseMatrix")

show_headline_part1 <- function(x)
{
    sprintf("<%s %s> of type \"%s\" ", paste0(dim(x), collapse=" x "),
                                       classNameForDisplay(x), type(x))
}

.show_nzcount <- function(x)
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
        ## example if the call to nzcount() in .show_nzcount() fails.
        cat(show_headline_part1(object))
        cat(.show_nzcount(object))
        if (any(dim(object) == 0L)) {
            cat("\n")
            return()
        }
        cat(":\n", sep="")
        S4Arrays:::print_some_array_elements(object)
    }
)

