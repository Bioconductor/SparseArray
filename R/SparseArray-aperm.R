### =========================================================================
### Transposition of a SparseArray object
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Transposition
###

### S3/S4 combo for t.SVT_SparseMatrix
t.SVT_SparseMatrix <- function(x)
{
    new_SVT <- .Call2("C_transpose_2D_SVT",
                      x@dim, x@type, x@SVT, PACKAGE="SparseArray")
    BiocGenerics:::replaceSlots(x, dim=rev(x@dim),
                                   dimnames=rev(x@dimnames),
                                   SVT=new_SVT,
                                   check=FALSE)
}
setMethod("t", "SVT_SparseMatrix", t.SVT_SparseMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Multidimensional transposition
###

.transpose_SVT <- function(x)
{
    stopifnot(is(x, "SVT_SparseArray"))
    new_SVT <- .Call2("C_transpose_SVT",
                      x@dim, x@type, x@SVT, PACKAGE="SparseArray")
    BiocGenerics:::replaceSlots(x, dim=rev(x@dim),
                                   dimnames=rev(x@dimnames),
                                   SVT=new_SVT,
                                   check=FALSE)
}

.transpose_SVT_v2 <- function(x)
{
    stopifnot(is(x, "SVT_SparseArray"))
    new_SVT <- .Call2("C_transpose_SVT_v2",
                      x@dim, x@type, x@SVT, PACKAGE="SparseArray")
    BiocGenerics:::replaceSlots(x, dim=rev(x@dim),
                                   dimnames=rev(x@dimnames),
                                   SVT=new_SVT,
                                   check=FALSE)
}

.transpose_SVT_v3 <- function(x)
{
    stopifnot(is(x, "SVT_SparseArray"))
    new_SVT <- .Call2("C_transpose_SVT_v3",
                      x@dim, x@type, x@SVT, PACKAGE="SparseArray")
    BiocGenerics:::replaceSlots(x, dim=rev(x@dim),
                                   dimnames=rev(x@dimnames),
                                   SVT=new_SVT,
                                   check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### aperm()
###

### Like aperm2() in S4Arrays, extend base::aperm() by allowing dropping
### and/or adding ineffective dimensions.
.aperm.COO_SparseArray <- function(a, perm)
{
    a_dim <- dim(a)
    perm <- S4Arrays:::normarg_perm(perm, a_dim)
    msg <- S4Arrays:::validate_perm(perm, a_dim)
    if (!isTRUE(msg))
        stop(wmsg(msg))
    ans_dim <- a_dim[perm]
    ans_dim[is.na(perm)] <- 1L
    ans_nzcoo <- a@nzcoo[ , perm, drop=FALSE]
    ans_nzcoo[ , is.na(perm)] <- 1L
    ans_dimnames <- a@dimnames[perm]
    new_COO_SparseArray(ans_dim, ans_dimnames,
                        ans_nzcoo, a@nzvals, check=FALSE)
}

### S3/S4 combo for aperm.COO_SparseArray
aperm.COO_SparseArray <-
    function(a, perm, ...) .aperm.COO_SparseArray(a, perm, ...)
setMethod("aperm", "COO_SparseArray", aperm.COO_SparseArray)

.aperm.SVT_SparseArray <- function(a, perm)
{
    perm <- S4Arrays:::normarg_perm(perm, a@dim)
    new_SVT <- .Call2("C_aperm_SVT",
                      a@dim, a@type, a@SVT, perm, PACKAGE="SparseArray")
    BiocGenerics:::replaceSlots(a, dim=a@dim[perm],
                                   dimnames=a@dimnames[perm],
                                   SVT=new_SVT,
                                   check=FALSE)
}

### S3/S4 combo for aperm.SVT_SparseArray
aperm.SVT_SparseArray <-
    function(a, perm, ...) .aperm.SVT_SparseArray(a, perm, ...)
setMethod("aperm", "SVT_SparseArray", aperm.SVT_SparseArray)

