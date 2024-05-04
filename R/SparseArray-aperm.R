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
    new_SVT <- SparseArray.Call("C_transpose_2D_SVT", x@dim, x@type, x@SVT)
    BiocGenerics:::replaceSlots(x, dim=rev(x@dim),
                                   dimnames=rev(x@dimnames),
                                   SVT=new_SVT,
                                   check=FALSE)
}
setMethod("t", "SVT_SparseMatrix", t.SVT_SparseMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### aperm()
###

### Like aperm2() in S4Arrays, extend base::aperm() by allowing dropping
### and/or adding ineffective dimensions.
.aperm_COO <- function(x, perm)
{
    stopifnot(is(x, "COO_SparseArray"))
    perm <- S4Arrays:::normarg_perm(perm, x@dim)
    msg <- S4Arrays:::validate_perm(perm, x@dim)
    if (!isTRUE(msg))
        stop(wmsg(msg))
    ans_dim <- x@dim[perm]
    ans_dim[is.na(perm)] <- 1L
    ans_nzcoo <- x@nzcoo[ , perm, drop=FALSE]
    ans_nzcoo[ , is.na(perm)] <- 1L
    ans_dimnames <- x@dimnames[perm]
    new_COO_SparseArray(ans_dim, ans_dimnames,
                        ans_nzcoo, x@nzdata, check=FALSE)
}

.aperm_SVT <- function(x, perm, .NAME=c("C_aperm_SVT", "C_aperm0_SVT"))
{
    stopifnot(is(x, "SVT_SparseArray"))
    check_svt_version(x)

    .NAME <- match.arg(.NAME)
    perm <- S4Arrays:::normarg_perm(perm, x@dim)
    new_SVT <- SparseArray.Call(.NAME, x@dim, x@type, x@SVT, perm)
    BiocGenerics:::replaceSlots(x, dim=x@dim[perm],
                                   dimnames=x@dimnames[perm],
                                   SVT=new_SVT,
                                   check=FALSE)
}

### S3/S4 combo for aperm.COO_SparseArray
aperm.COO_SparseArray <- function(a, perm, ...) .aperm_COO(a, perm, ...)
setMethod("aperm", "COO_SparseArray", aperm.COO_SparseArray)

### S3/S4 combo for aperm.SVT_SparseArray
aperm.SVT_SparseArray <- function(a, perm, ...) .aperm_SVT(a, perm, ...)
setMethod("aperm", "SVT_SparseArray", aperm.SVT_SparseArray)

