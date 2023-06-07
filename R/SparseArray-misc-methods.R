### =========================================================================
### Operate natively on SparseArray objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Various "unary isometric" array transformations
###
### A "unary isometric" array transformation is a transformation that returns
### an array-like object with the same dimensions as the input and where each
### element is the result of applying a function to the corresponding element
### in the input.
###
### Note that some "unary isometric" transformations preserve sparsity (e.g.
### is.na(), nchar(), round(), sqrt(), log1p(), etc...) and others don't
### (e.g. is.finite(), !, log(), etc..). We only implement the former.
###
### All the "unary isometric" array transformations implemented in this
### section return a COO_SparseArray object of the same dimensions as the
### input COO_SparseArray object.
###
### TODO: Add methods for SVT_SparseArray objects.

.UNARY_ISO_OPS <- c("is.na", "is.infinite", "is.nan", "tolower", "toupper")

for (.Generic in .UNARY_ISO_OPS) {
    setMethod(.Generic, "COO_SparseArray",
        function(x)
        {
            GENERIC <- match.fun(.Generic)
            new_nzvals <- GENERIC(x@nzvals)
            BiocGenerics:::replaceSlots(x, nzvals=new_nzvals, check=FALSE)
        }
    )
}

setMethod("nchar", "COO_SparseArray",
    function(x, type="chars", allowNA=FALSE, keepNA=NA)
    {
        new_nzvals <- nchar(x@nzvals, type=type, allowNA=allowNA, keepNA=keepNA)
        BiocGenerics:::replaceSlots(x, nzvals=new_nzvals, check=FALSE)
    }
)

