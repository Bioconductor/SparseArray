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
### Both .aperm_COO() and .aperm_SVT() support S4Arrays::aperm2() extended
### semantic.
###

.aperm_COO <- function(a, perm)
{
    stopifnot(is(a, "COO_SparseArray"))
    perm <- S4Arrays:::normarg_perm(perm, a@dim)
    msg <- S4Arrays:::validate_perm(perm, a@dim)
    if (!isTRUE(msg))
        stop(wmsg(msg))
    ans_dim <- a@dim[perm]
    ans_dim[is.na(perm)] <- 1L
    ans_nzcoo <- a@nzcoo[ , perm, drop=FALSE]
    ans_nzcoo[ , is.na(perm)] <- 1L
    ans_dimnames <- a@dimnames[perm]
    new_COO_SparseArray(ans_dim, ans_dimnames,
                        ans_nzcoo, a@nzdata, check=FALSE)
}

.aperm_SVT <- function(a, perm, .NAME=c("C_aperm_SVT", "C_aperm0_SVT"))
{
    stopifnot(is(a, "SVT_SparseArray"))
    check_svt_version(a)

    .NAME <- match.arg(.NAME)

    aperm0_SVT <- function(x, perm) {
        new_SVT <- SparseArray.Call(.NAME, x@dim, x@type, x@SVT, perm)
        BiocGenerics:::replaceSlots(x, dim=x@dim[perm],
                                       dimnames=x@dimnames[perm],
                                       SVT=new_SVT,
                                       check=FALSE)
    }
    S4Arrays:::extended_aperm(a, perm, aperm0_SVT)
}

### S3/S4 combo for aperm.COO_SparseArray
aperm.COO_SparseArray <- function(a, perm, ...) .aperm_COO(a, perm, ...)
setMethod("aperm", "COO_SparseArray", aperm.COO_SparseArray)

### S3/S4 combo for aperm.SVT_SparseArray
aperm.SVT_SparseArray <- function(a, perm, ...) .aperm_SVT(a, perm, ...)
setMethod("aperm", "SVT_SparseArray", aperm.SVT_SparseArray)

