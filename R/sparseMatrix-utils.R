### =========================================================================
### Low-level manipulation of sparseMatrix derivatives (from Matrix package)
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level CsparseMatrix and RsparseMatrix constructors
###
### NOT exported
###

### Returns "double" or "logical".
.infer_sparseMatrix_type_from_input_type <- function(input_type)
{
    switch(input_type, 'double'=, 'integer'=, 'raw'="double",
                       'logical'="logical",
           stop(wmsg("unsupported input type: ", input_type)))
}

new_CsparseMatrix <- function(dim, p, i, x, dimnames=NULL)
{
    stopifnot(is.integer(dim), length(dim) == 2L)
    x_type <- typeof(x)
    ans_type <- .infer_sparseMatrix_type_from_input_type(x_type)
    ans_class <- if (ans_type == "double") "dgCMatrix" else "lgCMatrix"
    if (ans_type != x_type)
        storage.mode(x) <- ans_type
    ans_dimnames <- S4Arrays:::normarg_dimnames(dimnames, dim)
    new(ans_class, Dim=dim, p=p, i=i, x=x, Dimnames=ans_dimnames)
}

new_RsparseMatrix <- function(dim, p, j, x, dimnames=NULL)
{
    stopifnot(is.integer(dim), length(dim) == 2L)
    x_type <- typeof(x)
    ans_type <- .infer_sparseMatrix_type_from_input_type(x_type)
    ans_class <- if (ans_type == "double") "dgRMatrix" else "lgRMatrix"
    if (ans_type != x_type)
        storage.mode(x) <- ans_type
    ans_dimnames <- S4Arrays:::normarg_dimnames(dimnames, dim)
    new(ans_class, Dim=dim, p=p, j=j, x=x, Dimnames=ans_dimnames)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### High-level CsparseMatrix and RsparseMatrix constructors
###
### NOT exported
###
### Note that CsparseMatrix() is a simpler version of Matrix::sparseMatrix()
### that is typically 50%-60% faster and more memory efficient.
###
### For both constructors:
### - 'i', 'j', 'nzdata' must be **parallel** atomic vectors.
### - 'i' and 'j' must be integer vectors with no NAs.
### - 'nzdata' must be a double, integer, raw, or logical vector, with no
###   zeros (NAs are ok).

CsparseMatrix <- function(dim, i, j, nzdata, dimnames=NULL)
{
    stopifnot(is.integer(dim), length(dim) == 2L,
              is.integer(i), is.integer(j), is.atomic(nzdata))
    oo <- order(j, i)
    ans_i <- i[oo] - 1L  # CsparseMatrix objects want this zero-based
    ans_p <- c(0L, cumsum(tabulate(j[oo], nbins=dim[[2L]])))
    ans_x <- nzdata[oo]
    new_CsparseMatrix(dim, ans_p, ans_i, ans_x, dimnames=dimnames)
}

RsparseMatrix <- function(dim, i, j, nzdata, dimnames=NULL)
{
    stopifnot(is.integer(dim), length(dim) == 2L,
              is.integer(i), is.integer(j), is.atomic(nzdata))
    oo <- order(i, j)
    ans_j <- j[oo] - 1L  # RsparseMatrix objects want this zero-based
    ans_p <- c(0L, cumsum(tabulate(i[oo], nbins=dim[[1L]])))
    ans_x <- nzdata[oo]
    new_RsparseMatrix(dim, ans_p, ans_j, ans_x, dimnames=dimnames)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion from Array derivative to sparseMatrix derivative
###

### These coercions will work out-of-the-box on any Array derivative that
### supports coercion to SparseArray.
.from_Array_to_sparseMatrix <- function(from, to)
{
    ## If 'from' is a SparseArray derivative, 'as(from, "SparseArray")' will
    ## be a no-op and thus doing 'as(as(from, "SparseArray"), to)' below will
    ## lead to an infinite recursion. We explicitly guard against this.
    if (is(from, "SparseArray"))
        stop(wmsg("coercion from ", class(from), " to ",
                  to, " is not supported"))
    as(as(from, "SparseArray"), to)
}

setAs("Array", "dgCMatrix",
    function(from) .from_Array_to_sparseMatrix(from, "dgCMatrix")
)
setAs("Array", "dgRMatrix",
    function(from) .from_Array_to_sparseMatrix(from, "dgRMatrix")
)
setAs("Array", "lgCMatrix",
    function(from) .from_Array_to_sparseMatrix(from, "lgCMatrix")
)
setAs("Array", "lgRMatrix",
    function(from) .from_Array_to_sparseMatrix(from, "lgRMatrix")
)

### These coercions will work out-of-the-box on any Array derivative that
### supports type() and coercion to [d|l]gCMatrix and to [d|l]gRMatrix.
setAs("Array", "CsparseMatrix",
    function(from)
    {
        ans_type <- .infer_sparseMatrix_type_from_input_type(type(from))
        ans_class <- if (ans_type == "double") "dgCMatrix" else "lgCMatrix"
        as(from, ans_class)
    }
)
setAs("Array", "RsparseMatrix",
    function(from)
    {
        ans_type <- .infer_sparseMatrix_type_from_input_type(type(from))
        ans_class <- if (ans_type == "double") "dgRMatrix" else "lgRMatrix"
        as(from, ans_class)
    }
)

### We give the preference to the CsparseMatrix representation (compressed
### column-oriented form).
setAs("Array", "sparseMatrix",
    function(from) as(from, "CsparseMatrix")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colMins_dgCMatrix()
### colMaxs_dgCMatrix()
### colRanges_dgCMatrix()
### colVars_dgCMatrix()
###
### NOT exported
###
### Don't turn these into formal S4 methods for dgCMatrix objects to avoid
### conflict with the methods defined in the sparseMatrixStats package!
### They do NOT propagate the colnames (the methods defined in matrixStats
### don't either).

colMins_dgCMatrix <- function (x, na.rm=FALSE)
{
    stopifnot(is(x, "dgCMatrix"))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    .Call2("C_colMins_dgCMatrix", x, na.rm, PACKAGE="SparseArray")
}

colMaxs_dgCMatrix <- function (x, na.rm=FALSE)
{
    stopifnot(is(x, "dgCMatrix"))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    .Call2("C_colMaxs_dgCMatrix", x, na.rm, PACKAGE="SparseArray")
}

### About 2x faster than the method for dgCMatrix objects defined
### in sparseMatrixStats.
colRanges_dgCMatrix <- function (x, na.rm=FALSE)
{
    stopifnot(is(x, "dgCMatrix"))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    .Call2("C_colRanges_dgCMatrix", x, na.rm, PACKAGE="SparseArray")
}

### About 2.5x faster than the method for dgCMatrix objects defined
### in sparseMatrixStats.
colVars_dgCMatrix <- function(x, na.rm=FALSE)
{
    stopifnot(is(x, "dgCMatrix"))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    .Call2("C_colVars_dgCMatrix", x, na.rm, PACKAGE="SparseArray")
}

