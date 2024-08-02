### =========================================================================
### Transposition of an NaArray object
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Transposition
###

### S3/S4 combo for t.NaMatrix
t.NaMatrix <- function(x)
{
    check_svt_version(x)
    new_NaSVT <- SparseArray.Call("C_transpose_2D_SVT", x@dim, x@type, x@NaSVT)
    BiocGenerics:::replaceSlots(x, dim=rev(x@dim),
                                   dimnames=rev(x@dimnames),
                                   NaSVT=new_NaSVT,
                                   check=FALSE)
}
setMethod("t", "NaMatrix", t.NaMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### aperm()
###
### Supports S4Arrays::aperm2() extended semantic.
###

.aperm_NaSVT <- function(a, perm, .NAME=c("C_aperm_SVT", "C_aperm0_SVT"))
{
    stopifnot(is(a, "NaArray"))
    check_svt_version(a)

    .NAME <- match.arg(.NAME)

    aperm0_NaSVT <- function(x, perm) {
        new_NaSVT <- SparseArray.Call(.NAME, x@dim, x@type, x@NaSVT, perm)
        BiocGenerics:::replaceSlots(x, dim=x@dim[perm],
                                       dimnames=x@dimnames[perm],
                                       NaSVT=new_NaSVT,
                                       check=FALSE)
    }
    S4Arrays:::extended_aperm(a, perm, aperm0_NaSVT)
}

### S3/S4 combo for aperm.NaArray
aperm.NaArray <- function(a, perm, ...) .aperm_NaSVT(a, perm, ...)
setMethod("aperm", "NaArray", aperm.NaArray)

