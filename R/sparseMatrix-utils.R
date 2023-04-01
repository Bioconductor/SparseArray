### =========================================================================
### Low-level manipulation of sparseMatrix derivatives
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
### CsparseMatrix() -- NOT exported
###
### A simpler version of Matrix::sparseMatrix() that is typically 50%-60%
### faster and more memory efficient.

### 'i', 'j', 'nzvals' must be **parallel** atomic vectors.
### 'i' and 'j' must be integer vectors with no NAs.
### 'nzvals' must be a double, integer, raw, or logical vector, with no
### zeros and possibly with NAs.
CsparseMatrix <- function(dim, i, j, nzvals, dimnames=NULL)
{
    stopifnot(is.integer(dim), length(dim) == 2L,
              is.integer(i), is.integer(j), is.atomic(nzvals))
    oo <- order(j, i)
    ans_i <- i[oo] - 1L  # CsparseMatrix objects want this zero-based
    ans_p <- c(0L, cumsum(tabulate(j[oo], nbins=dim[[2L]])))
    ans_x <- nzvals[oo]
    new_CsparseMatrix(dim, ans_p, ans_i, ans_x, dimnames=dimnames)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion from Array to sparseMatrix
###

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
### RsparseMatrix() -- NOT exported
###

RsparseMatrix <- function(dim, i, j, nzvals, dimnames=NULL)
{
    stopifnot(is.integer(dim), length(dim) == 2L,
              is.integer(i), is.integer(j), is.atomic(nzvals))
    oo <- order(i, j)
    ans_j <- j[oo] - 1L  # RsparseMatrix objects want this zero-based
    ans_p <- c(0L, cumsum(tabulate(i[oo], nbins=dim[[1L]])))
    ans_x <- nzvals[oo]
    new_RsparseMatrix(dim, ans_p, ans_j, ans_x, dimnames=dimnames)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rowsum() method for dgCMatrix objects
###

.rowsum_dgCMatrix <- function(x, group, reorder=TRUE, na.rm=FALSE)
{
    stopifnot(is(x, "dgCMatrix"))
    ugroup <- S4Arrays:::compute_ugroup(group, nrow(x), reorder)
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    group <- match(group, ugroup)
    ans <- .Call2("C_rowsum_dgCMatrix", x, group, length(ugroup), na.rm,
                                        PACKAGE="SparseArray")
    dimnames(ans) <- list(as.character(ugroup), colnames(x))
    ans
}

### S3/S4 combo for rowsum.dgCMatrix
rowsum.dgCMatrix <- function(x, group, reorder=TRUE, ...)
    .rowsum_dgCMatrix(x, group, reorder=reorder, ...)
setMethod("rowsum", "dgCMatrix", rowsum.dgCMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colMins_dgCMatrix()
### colMaxs_dgCMatrix()
### colRanges_dgCMatrix()
### colVars_dgCMatrix()
###
### NOT exported.
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

