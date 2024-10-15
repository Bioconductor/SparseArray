### =========================================================================
### Low-level manipulation of sparseMatrix derivatives (from Matrix package)
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level [C|R|T]sparseMatrix constructors
###
### WARNING: They take **0-based** row/column indices!
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

### 'i' must be an integer vector containing 0-based row indices.
### Returns a dgCMatrix, lgCMatrix, or ngCMatrix (depends on 'typeof(x)').
new_CsparseMatrix <- function(dim, p, i, x=NULL, dimnames=NULL)
{
    stopifnot(is.integer(dim), length(dim) == 2L)
    ans_dimnames <- S4Arrays:::normarg_dimnames(dimnames, dim)
    if (is.null(x))
        return(new("ngCMatrix", Dim=dim, p=p, i=i, Dimnames=ans_dimnames))
    x_type <- typeof(x)
    ans_type <- .infer_sparseMatrix_type_from_input_type(x_type)
    ans_class <- if (ans_type == "double") "dgCMatrix" else "lgCMatrix"
    if (ans_type != x_type)
        storage.mode(x) <- ans_type
    new(ans_class, Dim=dim, p=p, i=i, x=x, Dimnames=ans_dimnames)
}

### 'j' must be an integer vector containing 0-based column indices.
### Returns a dgRMatrix, lgRMatrix, or ngRMatrix (depends on 'typeof(x)').
new_RsparseMatrix <- function(dim, p, j, x=NULL, dimnames=NULL)
{
    stopifnot(is.integer(dim), length(dim) == 2L)
    ans_dimnames <- S4Arrays:::normarg_dimnames(dimnames, dim)
    if (is.null(x))
        return(new("ngRMatrix", Dim=dim, p=p, j=j, Dimnames=ans_dimnames))
    x_type <- typeof(x)
    ans_type <- .infer_sparseMatrix_type_from_input_type(x_type)
    ans_class <- if (ans_type == "double") "dgRMatrix" else "lgRMatrix"
    if (ans_type != x_type)
        storage.mode(x) <- ans_type
    new(ans_class, Dim=dim, p=p, j=j, x=x, Dimnames=ans_dimnames)
}

### 'i' and 'j' must be **parallel** integer vectors containing 0-based row
### and column indices, respectively.
### Returns a dgTMatrix, lgTMatrix, or ngTMatrix (depends on 'typeof(x)').
new_TsparseMatrix <- function(dim, i, j, x=NULL, dimnames=NULL)
{
    stopifnot(is.integer(dim), length(dim) == 2L)
    ans_dimnames <- S4Arrays:::normarg_dimnames(dimnames, dim)
    if (is.null(x))
        return(new("ngTMatrix", Dim=dim, i=i, j=j, Dimnames=ans_dimnames))
    x_type <- typeof(x)
    ans_type <- .infer_sparseMatrix_type_from_input_type(x_type)
    ans_class <- if (ans_type == "double") "dgTMatrix" else "lgTMatrix"
    if (ans_type != x_type)
        storage.mode(x) <- ans_type
    new(ans_class, Dim=dim, i=i, j=j, x=x, Dimnames=ans_dimnames)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### High-level [C|R|T]sparseMatrix constructors
###
### For the 3 constructors:
### - 'i' and 'j' must be **parallel** integer vectors that contain 1-based
###   row and column indices, respectively.
### - 'nzdata' must be NULL or an atomic vector **parallel** to 'i' and 'j'.
###   Its type must be "double", "integer", "raw", or "logical". It should
###   not contain zeros (NAs are ok) but this is not checked.
### - Each constructor takes care of removing triplets with duplicated (i,j)
###   coordinates from the supplied (i,j,nzdata) triplets.
###
### NOT exported
###

.check_i_j_nzdata <- function(i, j, nzdata)
{
    stopifnot(is.integer(i), is.integer(j), length(i) == length(j))
    if (!is.null(nzdata))
        stopifnot(is.atomic(nzdata), length(nzdata) == length(i))
}

.drop_triplets_with_dup_i_j_pairs <- function(i, j, nzdata=NULL)
{
    .check_i_j_nzdata(i, j, nzdata)
    dup_idx <- which(duplicatedIntegerPairs(i, j, fromLast=TRUE))
    if (length(dup_idx) != 0L) {
        i <- i[-dup_idx]
        j <- j[-dup_idx]
        if (!is.null(nzdata))
            nzdata <- nzdata[-dup_idx]
    }
    list(i, j, nzdata)
}

### Note that this is a simplified version of Matrix::sparseMatrix() that
### will typically be 50%-60% faster and more memory efficient.
CsparseMatrix <- function(dim, i, j, nzdata=NULL, dimnames=NULL)
{
    stopifnot(is.integer(dim), length(dim) == 2L)
    triplets <- .drop_triplets_with_dup_i_j_pairs(i, j, nzdata=nzdata)
    i <- triplets[[1L]]
    j <- triplets[[2L]]
    nzdata <- triplets[[3L]]
    oo <- orderIntegerPairs(j, i)
    ans_p <- c(0L, cumsum(tabulate(j[oo], nbins=dim[[2L]])))
    ans_i <- i[oo] - 1L  # new_CsparseMatrix() wants this 0-based
    ans_x <- if (is.null(nzdata)) NULL else nzdata[oo]
    new_CsparseMatrix(dim, ans_p, ans_i, x=ans_x, dimnames=dimnames)
}

RsparseMatrix <- function(dim, i, j, nzdata=NULL, dimnames=NULL)
{
    stopifnot(is.integer(dim), length(dim) == 2L)
    triplets <- .drop_triplets_with_dup_i_j_pairs(i, j, nzdata=nzdata)
    i <- triplets[[1L]]
    j <- triplets[[2L]]
    nzdata <- triplets[[3L]]
    oo <- orderIntegerPairs(i, j)
    ans_p <- c(0L, cumsum(tabulate(i[oo], nbins=dim[[1L]])))
    ans_j <- j[oo] - 1L  # new_RsparseMatrix() wants this 0-based
    ans_x <- if (is.null(nzdata)) NULL else nzdata[oo]
    new_RsparseMatrix(dim, ans_p, ans_j, x=ans_x, dimnames=dimnames)
}

### Besides optional removal of the triplets with duplicated (i,j) coordinates,
### the supplied (i,j,nzdata) triplets are stored as-is in the returned
### TsparseMatrix object, that is, their original order is preserved.
TsparseMatrix <- function(dim, i, j, nzdata=NULL, dimnames=NULL, drop.dups=TRUE)
{
    stopifnot(is.integer(dim), length(dim) == 2L, isTRUEorFALSE(drop.dups))
    if (drop.dups) {
        triplets <- .drop_triplets_with_dup_i_j_pairs(i, j, nzdata=nzdata)
        i <- triplets[[1L]]
        j <- triplets[[2L]]
        nzdata <- triplets[[3L]]
    } else {
        .check_i_j_nzdata(i, j, nzdata)
    }
    new_TsparseMatrix(dim, i - 1L, j - 1L, x=nzdata, dimnames=dimnames)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Critically endangered coercions from/to sparseMatrix derivatives
###
### This simply brings back some basic coercion methods originally defined
### in the Matrix package but that the lazy Matrix maintainers have decided
### to eradicate from the surface of Earth in a shocking attempt at making
### their users lives miserable.
###

### --- from ordinary matrix to [d|l|n]gCMatrix and [d|l|n]gRMatrix ---

### Not deprecated yet. Cold feet maybe?
#setAs("matrix", "dgCMatrix",
#    function(from)
#        as(as(as(from, "dMatrix"), "generalMatrix"), "CsparseMatrix")
#)

### Deprecated in Matrix 1.7-0
setAs("matrix", "dgRMatrix",
    function(from)
        as(as(as(from, "dMatrix"), "generalMatrix"), "RsparseMatrix")
)

### Deprecated in Matrix 1.7-0
setAs("matrix", "lgCMatrix",
    function(from)
        as(as(as(from, "lMatrix"), "generalMatrix"), "CsparseMatrix")
)

### Never supported by Matrix?
setAs("matrix", "lgRMatrix",
    function(from)
        as(as(as(from, "lMatrix"), "generalMatrix"), "RsparseMatrix")
)

### Deprecated in Matrix 1.7-0
setAs("matrix", "ngCMatrix",
    function(from)
        as(as(as(from, "nMatrix"), "generalMatrix"), "CsparseMatrix")
)

### Never supported by Matrix?
setAs("matrix", "ngRMatrix",
    function(from)
        as(as(as(from, "nMatrix"), "generalMatrix"), "RsparseMatrix")
)

### --- other useful coercions ---

### Deprecated in Matrix 1.7-0
setAs("dgCMatrix", "ngCMatrix", function(from) as(from, "nMatrix"))

### Never supported by Matrix?
setAs("CsparseMatrix", "ngCMatrix", function(from) as(from, "nMatrix"))
setAs("RsparseMatrix", "ngRMatrix", function(from) as(from, "nMatrix"))
setAs("TsparseMatrix", "ngTMatrix", function(from) as(from, "nMatrix"))


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
setAs("Array", "ngCMatrix",
    function(from) .from_Array_to_sparseMatrix(from, "ngCMatrix")
)
setAs("Array", "ngRMatrix",
    function(from) .from_Array_to_sparseMatrix(from, "ngRMatrix")
)

### These coercions will work out-of-the-box on any Array derivative that
### supports type() and coercion to [d|l]gCMatrix, [d|l]gRMatrix, and
### [d|l]gTMatrix, respectively.
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
setAs("Array", "TsparseMatrix",
    function(from)
    {
        ans_type <- .infer_sparseMatrix_type_from_input_type(type(from))
        ans_class <- if (ans_type == "double") "dgTMatrix" else "lgTMatrix"
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
### IMPORTANT NOTE: The functions below precede the sparseMatrixStats package
### so we probaby don't need them anymore. Do NOT turn them into formal S4
### methods for dgCMatrix objects to avoid conflicts with the methods defined
### in the sparseMatrixStats package!
### They do NOT propagate the colnames. The corresponding methods defined in
### matrixStats didn't either at the time this was implemented.

colMins_dgCMatrix <- function (x, na.rm=FALSE)
{
    stopifnot(is(x, "dgCMatrix"))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    SparseArray.Call("C_colMins_dgCMatrix", x, na.rm)
}

colMaxs_dgCMatrix <- function (x, na.rm=FALSE)
{
    stopifnot(is(x, "dgCMatrix"))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    SparseArray.Call("C_colMaxs_dgCMatrix", x, na.rm)
}

### About 2x faster than the method for dgCMatrix objects defined
### in sparseMatrixStats.
colRanges_dgCMatrix <- function (x, na.rm=FALSE)
{
    stopifnot(is(x, "dgCMatrix"))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    SparseArray.Call("C_colRanges_dgCMatrix", x, na.rm)
}

### About 2.5x faster than the method for dgCMatrix objects defined
### in sparseMatrixStats.
colVars_dgCMatrix <- function(x, na.rm=FALSE)
{
    stopifnot(is(x, "dgCMatrix"))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    SparseArray.Call("C_colVars_dgCMatrix", x, na.rm)
}

