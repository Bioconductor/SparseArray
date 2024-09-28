### =========================================================================
### Miscellaneous operations on SparseArray objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Various "unary isometric" array transformations
###
### A "unary isometric" array transformation is a transformation that returns
### an array-like object with the same dimensions as the input and where each
### element is the result of applying a function to the corresponding element
### in the input.
###
### Note that:
### - Some "unary isometric" transformations preserve sparsity (e.g. is.na(),
###   nchar(), round(), sqrt(), log1p(), etc...) and others don't (e.g.
###   is.finite(), !, log(), etc..). SparseArray objects only need to support
###   the former.
### - All operations from the 'Math' and 'Math2' groups are "unary isometric"
###   transformations (see '?S4groupGeneric'). The corresponding methods for
###   SparseArray objects are implemented in R/SparseArray-Math-methods.R
### - All the "unary isometric" methods implemented below return an array-like
###   object of the same class as the input (endomorphism).

### --- Methods for COO_SparseArray objects ---

.isoFUN_COO <- function(isoFUN, x, ...)
{
    GENERIC <- match.fun(isoFUN)
    new_nzdata <- GENERIC(x@nzdata, ...)
    BiocGenerics:::replaceSlots(x, nzdata=new_nzdata, check=FALSE)
}

setMethod("is.na", "COO_SparseArray",
    function(x) .isoFUN_COO("is.na", x)
)
setMethod("is.nan", "COO_SparseArray",
    function(x) .isoFUN_COO("is.nan", x)
)
setMethod("is.infinite", "COO_SparseArray",
    function(x) .isoFUN_COO("is.infinite", x)
)
setMethod("tolower", "COO_SparseArray",
    function(x) .isoFUN_COO("tolower", x)
)
setMethod("toupper", "COO_SparseArray",
    function(x) .isoFUN_COO("toupper", x)
)
setMethod("nchar", "COO_SparseArray",
    function(x, type="chars", allowNA=FALSE, keepNA=NA)
        .isoFUN_COO("nchar", x, type=type, allowNA=allowNA, keepNA=keepNA)
)

### --- Methods for SVT_SparseArray objects ---

### Returns a "logical" SVT_SparseArray object.
.isFUN_SVT <- function(isFUN, x)
{
    stopifnot(is(x, "SVT_SparseArray"))
    check_svt_version(x)
    new_SVT <- SparseArray.Call("C_SVT_apply_isFUN",
                                x@dim, x@type, x@SVT, isFUN)
    BiocGenerics:::replaceSlots(x, type="logical", SVT=new_SVT, check=FALSE)
}

setMethod("is.na", "SVT_SparseArray",
    function(x) .isFUN_SVT("is.na", x)
)
setMethod("is.nan", "SVT_SparseArray",
    function(x) .isFUN_SVT("is.nan", x)
)
setMethod("is.infinite", "SVT_SparseArray",
    function(x) .isFUN_SVT("is.infinite", x)
)

### TODO: Support more methods!


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Various "N-ary isometric" array transformations
###
### An "N-ary isometric" array transformation is a transformation that takes
### one or more array-like objects of the same dimensions (a.k.a. conformable
### arrays) and returns an array-like object of the same dimensions.
###
### Note that:
### - All operations from the 'Ops' group are "N-ary isometric"
###   transformations (see '?S4groupGeneric'). The corresponding
###   methods for SparseArray objects are implemented in files
###   R/SparseArray-[Arith|Compare|Logic]-methods.R.
### - If all the input arrays have the same class then the "N-ary isometric"
###   methods implemented below return an array-like object of that class
###   (endomorphism).

### Binary pmin() and pmax() between two conformable array-like objects.
### We go for a pure R implementation for now that relies on `<` (or `>`),
### nzwhich(), linear subsetting (`[`), linear subassignment (`[<-`), and
### is.na(), so it works on any array-like object that supports these
### operations e.g. on ordinary arrays or dgCMatrix objects. However it
### won't work on COO_SparseArray objects because these objects don't support
### some of the required operations like `<` or linear subsetting. Sparsity
### is preserved all along so it's pretty efficient. Of course nothing would
### beat a C implementation but this is good enough for now.
### About base::pmin() and base::pmax(): These are also pure R implementations
### that would **almost** work on SVT_SparseArray objects except that they
### rely on logical negation (!) of the array which SVT_SparseArray objects
### don't support because it does NOT preserve sparsity. However this what
### pmin()/pmax() use on dgCMatrix objects (there are no dedicated pmin()
### or pmax() methods for dgCMatrix objects as of Matrix 1.7-0), which is
### not very efficient.
.pminmax2 <- function(op, x, y, na.rm=FALSE)
{
    stopifnot(is(x, "SVT_SparseArray"), is(y, "SVT_SparseArray"))
    op <- match.fun(op)
    ans <- x
    replace_idx <- nzwhich(op(y, x))
    ## This subassignment will take care of setting the type of 'ans'
    ## to the "biggest" of 'type(x)' and 'type(y)' so it's a good
    ## idea to do it even if there's nothing to replace (i.e. even
    ## if 'length(replace_idx)' is 0).
    ans[replace_idx] <- y[replace_idx]

    if (na.rm) {
        ## The above replacement propagated NAs from 'y' to 'ans'.
        ## We need to revert that.
        is_na <- is.na(y)
    } else {
        ## The above replacement may have replaced NAs in 'ans' with
        ## non-NA values from 'y'. We need to revert that.
        is_na <- is.na(x)
    }
    restore_idx <- nzwhich(is_na)
    ans[restore_idx] <- x[restore_idx]

    ans_dimnames <- S4Arrays:::get_first_non_NULL_dimnames(list(x, y))
    S4Arrays:::set_dimnames(ans, ans_dimnames)
}

.pmin2 <- function(x, y, na.rm=FALSE) .pminmax2("<", x, y, na.rm=na.rm)
.pmax2 <- function(x, y, na.rm=FALSE) .pminmax2(">", x, y, na.rm=na.rm)

.psummarize <- function(NaryFUN, binaryFUN, ..., na.rm=FALSE)
{
    NaryFUN <- match.fun(NaryFUN)
    binaryFUN <- match.fun(binaryFUN)
    objects <- S4Vectors:::delete_NULLs(list(...))
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))

    if (length(objects) == 0L)
        stop(wmsg("no input"))  # should never happen

    x <- objects[[1L]]
    if (length(objects) == 1L)
        return(x)   # no-op

    if (length(objects) == 2L) {
        y <- objects[[2L]]
    } else {
        ## Recursive.
        y <- do.call(NaryFUN, c(objects[-1L], list(na.rm=na.rm)))
    }
    binaryFUN(x, y, na.rm=na.rm)
}

### Method dispatch will select these methods if and only if **all** the
### objects passed thru the ellipsis (...) are SparseArray derivatives.
### Note that even though the methods are defined for SparseArray objects,
### they will fail on COO_SparseArray objects (see comment for "Binary pmin()
### and pmax() between two conformable array-like objects" above).
setMethod("pmin", "SparseArray",
    function(..., na.rm=FALSE) .psummarize("pmin", .pmin2, ..., na.rm=na.rm)
)
setMethod("pmax", "SparseArray",
    function(..., na.rm=FALSE) .psummarize("pmax", .pmax2, ..., na.rm=na.rm)
)

