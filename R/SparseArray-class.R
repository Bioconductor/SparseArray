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
###    - nzcount(), nzwhich(), nzvals(), and `nzvals<-`() generics
###    - sparsity()
### 2) Implemented elsewhere:
###    - nzcount(), nzwhich(), nzvals(), and `nzvals<-`() methods
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
### nzcount(), nzwhich(), nzvals(), and `nzvals<-`() generics
###

### Returns the number of nonzero array elements in 'x'.
setGeneric("nzcount", function(x) standardGeneric("nzcount"))

### Not 100% reliable on [d|l]gCMatrix objects because these objects are
### allowed to have zeros in their @x slot!
### See src/SVT_SparseArray_class.c for an example.
setMethod("nzcount", "CsparseMatrix", function(x) length(x@i))
setMethod("nzcount", "RsparseMatrix", function(x) length(x@j))

### Returns the indices of the nonzero array elements in 'x', either as
### an L-index (if 'arr.ind=FALSE') or as an M-index (if 'arr.ind=TRUE').
setGeneric("nzwhich", signature="x",
    function(x, arr.ind=FALSE) standardGeneric("nzwhich")
)

### Works on any vector-like or array-like object 'x' that supports type(x)
### and comparison with zero (x != zero). One notable exception are ngRMatrix
### objects from the Matrix package ('x != FALSE' is broken on these objects).
default_nzwhich <- function(x, arr.ind=FALSE)
{
    if (!isTRUEorFALSE(arr.ind))
        stop(wmsg("'arr.ind' must be TRUE or FALSE"))
    ## Make sure to use 'type()' and not 'typeof()'.
    zero <- vector(type(x), length=1L)
    is_nonzero <- x != zero  # broken on ngRMatrix objects!
    which(is_nonzero | is.na(is_nonzero), arr.ind=arr.ind, useNames=FALSE)
}
setMethod("nzwhich", "ANY", default_nzwhich)

### default_nzwhich() above works fine on [C|R]sparseMatrix derivatives
### (except on ngRMatrix objects) but nzwhich_CsparseMatrix() and
### nzwhich_RsparseMatrix() below are more efficient:
### - nzwhich_CsparseMatrix() is typically 10x to 50x faster than
###   default_nzwhich() on big CsparseMatrix objects;
### - nzwhich_RsparseMatrix() is only slightly faster than default_nzwhich()
###   but is provided mostly to make nzwhich() work on ngRMatrix objects
###   (default_nzwhich() doesn't work on these objects, see above).
### IMPORTANT NOTE: nzwhich_CsparseMatrix() and nzwhich_RsparseMatrix()
### use a shortcut that is **NOT** 100% reliable on [d|l]gCMatrix or
### [d|l]gRMatrix objects because these objects sometimes have zeros in
### their @x slot (see src/SVT_SparseArray_class.c for examples of such
### objects). So in this case the functions will produce a result that
### contains some false positives.
nzwhich_CsparseMatrix <- function(x, arr.ind=FALSE)
{
    if (!isTRUEorFALSE(arr.ind))
        stop(wmsg("'arr.ind' must be TRUE or FALSE"))
    x_nrow <- nrow(x)
    x_ncol <- ncol(x)
    ## If 'x' is a "long object" (i.e. length(x) >= 2^31) then we'll get
    ## an integer overflow in 'x_nrow * (seq_len(x_ncol) - 1L)' below.
    ## Coercing 'x_nrow' to double avoids that.
    if (is.double(length(x)))
        x_nrow <- as.double(x_nrow)
    offsets <- rep.int(x_nrow * (seq_len(x_ncol) - 1L), diff(x@p))
    ans <- x@i + offsets + 1L
    if (!arr.ind)
        return(ans)
    Lindex2Mindex(ans, dim(x))
}
nzwhich_RsparseMatrix <- function(x, arr.ind=FALSE)
{
    if (!isTRUEorFALSE(arr.ind))
        stop(wmsg("'arr.ind' must be TRUE or FALSE"))
    x_nrow <- nrow(x)
    x_ncol <- ncol(x)
    if (is.double(length(x)))
        x_ncol <- as.double(x_ncol)  # see nzwhich_CsparseMatrix() above
    offsets <- rep.int(x_ncol * (seq_len(x_nrow) - 1L), diff(x@p))
    transposed_nzwhich <- x@j + offsets + 1L
    transposed_arr_ind <- Lindex2Mindex(transposed_nzwhich, rev(dim(x)))
    transposed_arr_ind1 <- transposed_arr_ind[ , 1L]
    transposed_arr_ind2 <- transposed_arr_ind[ , 2L]
    oo <- order(transposed_arr_ind1, transposed_arr_ind2)
    arr_ind1 <- transposed_arr_ind2[oo]
    arr_ind2 <- transposed_arr_ind1[oo]
    arr_ind <- cbind(arr_ind1, arr_ind2, deparse.level=0)
    if (arr.ind)
        return(arr_ind)
    Mindex2Lindex(arr_ind, dim(x))
}

setMethod("nzwhich", "CsparseMatrix", nzwhich_CsparseMatrix)
setMethod("nzwhich", "RsparseMatrix", nzwhich_RsparseMatrix)

### Returns the values of the nonzero array elements in a vector of the
### same type() as 'x' and parallel to nzwhich(x).
### Equivalent to 'x[nzwhich(x)]' (and that's what the default method
### below does). However specialized methods have the potential to make
### this dramatically faster.
setGeneric("nzvals", function(x) standardGeneric("nzvals"))

### Assumes that array-like object 'x' supports subsetting by a linear
### numeric subscript (L-index). An alternative would be to use an M-index
### instead, that is, to subset by 'nzwhich(x, arr.ind=TRUE)' instead
### of 'nzwhich(x)'.
### Note that the call to as.vector() should not be necessary since
### the result of subsetting 'x' by an L- or M-index should already
### be a vector (atomic or list). However, in the case of a 1D ordinary
### array, it's not! For example, 'array(11:15)[matrix(3:1)]' is a 1D array.
### This is a bug in base::`[`.
setMethod("nzvals", "ANY", function(x) as.vector(x[nzwhich(x)]))

### Not 100% reliable because [d|l]gCMatrix objects can have zeros in
### their @x slot (see comment for nzwhich_CsparseMatrix() above), but
### consistent with nzwhich_CsparseMatrix().
setMethod("nzvals", "dgCMatrix", function(x) x@x)
setMethod("nzvals", "lgCMatrix", function(x) x@x)

setMethod("nzvals", "ngCMatrix", function(x) rep.int(TRUE, length(x@i)))
setMethod("nzvals", "ngRMatrix", function(x) rep.int(TRUE, length(x@j)))

### Replace the values of the nonzero array elements in 'x'.
### Equivalent to 'x[nzwhich(x)] <- value' (and that's what the default
### method below does). However specialized methods have the potential to
### make this dramatically faster.
setGeneric("nzvals<-", signature="x",
    function(x, value) standardGeneric("nzvals<-")
)

setReplaceMethod("nzvals", "ANY",
    function(x, value)
    {
        if (!is.vector(value))
            stop(wmsg("replacement value must be a vector"))
        x[nzwhich(x)] <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sparsity()
###

sparsity <- function(x) { 1 - nzcount(x) / length(x) }


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

