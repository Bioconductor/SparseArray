### =========================================================================
### matrixStats methods for SparseMatrix and SparseArray objects
### -------------------------------------------------------------------------
###
### About matrixStats usage in Bioconductor: Based on some quick grep-based
### inspection, the matrixStats operations used by Bioconductor software
### packages are (looking at the col* functions only):
###   (1) Heavily used: colSums, colMeans, colMedians, colVars, colSds,
###       colMaxs, colMins, colMeans2, colSums2
###   (2) Not so heavily used: colRanges, colRanks, colQuantiles, colMads,
###       colIQRs
###   (3) Marginally used: colAlls, colCumsums, colWeightedMeans, colAnyNAs
###
### Notes:
### - colSums() and colMeans() are functions actually defined in the base
###   package but we still count them as part of the matrixStats family.
### - All other matrix col/row summarization operations are from the
###   matrixStats package with corresponding generics defined in the
###   MatrixGenerics package.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers
###

### A silly trick used only to trigger an error if called with arguments.
.check_unused_arguments <- function() NULL

.check_dims <- function(dims, method)
{
    if (!identical(dims, 1))
        stop(wmsg("the ", method, "() method for SVT_SparseArray objects ",
                  "does not support the 'dims' argument"))
}

.check_rows_cols <- function(rows, cols, method)
{
    if (!(is.null(rows) && is.null(cols)))
        stop(wmsg("the ", method, "() method for SVT_SparseArray objects ",
                  "does not support the 'rows' or 'cols' argument"))
}

### Returns TRUE or FALSE.
.normarg_useNames <- function(useNames=NA)
{
    if (!(is.logical(useNames) && length(useNames) == 1L))
        stop(wmsg("'useNames' must be a single logical value"))
    !isFALSE(useNames)
}

.stopifnot_2D_object <- function(x, method)
{
    if (length(dim(x)) != 2L)
        stop(wmsg("the ", method, "() method for SVT_SparseArray objects ",
                  "only supports 2D objects (i.e. SVT_SparseMatrix objects) ",
                  "at the moment"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .colStats_SparseArray(), .rowStats_SparseArray()
###
### Workhorses behind all the matrixStats methods for SVT_SparseArray
### objects, with the exception of the colMedians()/rowMedians() methods
### at the moment.
###

.normarg_dims <- function(dims)
{
    if (!isSingleNumber(dims))
        stop(wmsg("'dims' must be a single integer"))
    if (!is.integer(dims))
        dims <- as.integer(dims)
    dims
}

### Return an ordinary array with 'length(dim(x)) - dims' dimensions.
.colStats_SparseArray <- function(op, x, na.rm=FALSE, center=NULL, dims=1L,
                                  useNames=NA)
{
    stopifnot(isSingleString(op), is(x, "SparseArray"))
    if (is(x, "SVT_SparseArray")) {
        check_svt_version(x)
    } else {
        x <- as(x, "SVT_SparseArray")
    }

    ## Check 'na.rm'.
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))

    ## Check and normalize 'center'.
    if (is.null(center)) {
        center <- NA_real_
    } else {
        if (!isSingleNumberOrNA(center))
            stop(wmsg("'center' must be NULL, or a single number"))
        if (!is.double(center))
            center <- as.double(center)
    }

    ## Check and normalize 'dims'.
    dims <- .normarg_dims(dims)
    if (dims <= 0L || dims > length(x@dim))
        stop(wmsg("'dims' must be a single integer that is ",
                  "> 0 and <= length(dim(x)) for the col*() functions, and ",
                  ">= 0 and < length(dim(x)) for the row*() functions"))

    ## Normalize 'useNames'.
    useNames <- .normarg_useNames(useNames)

    x_dimnames <- if (useNames) x@dimnames else NULL
    SparseArray.Call("C_colStats_SVT",
                     x@dim, x_dimnames, x@type, x@SVT,
                     op, na.rm, center, dims)
}

### WORK-IN-PROGRESS!
### The naive implementation for .rowStats_SparseArray() would be to simply do:
###
###   aperm(.colStats_SparseArray(op, aperm(x), ..., dims=length(dim(x))-dims))
###
### This is semantically correct for any number of dimensions. However,
### it is VERY inefficient when 'x' is a SVT_SparseArray object with more
### than 2 dimensions because multidimensional transposition of 'x' (i.e.
### 'aperm(x)') is VERY expensive in that case. So we use some tricks below
### to avoid this multidimensional transposition.
.rowStats_SparseArray <- function(op, x, na.rm=FALSE, center=NULL, dims=1L,
                                  useNames=NA)
{
    stopifnot(isSingleString(op), is(x, "SparseArray"))

    ## Check and normalize 'dims'.
    dims <- .normarg_dims(dims)
    if (dims < 0L || dims >= length(x@dim))
        stop(wmsg("'dims' must be a single integer that is ",
                  "> 0 and <= length(dim(x)) for the col*() functions, and ",
                  ">= 0 and < length(dim(x)) for the row*() functions"))

    ## Normalize 'useNames'.
    useNames <- .normarg_useNames(useNames)

    x_ndim <- length(x@dim)

    if (is(x, "COO_SparseArray")) {
        ans <- .colStats_SparseArray(op, as(aperm(x), "SVT_SparseArray"),
                                     na.rm=na.rm, center=center,
                                     dims=x_ndim-dims, useNames=useNames)
        if (!is.null(dim(ans)))
            ans <- aperm(ans)
        return(ans)
    }

    check_svt_version(x)
    if (x_ndim <= 2L || dims == 0L || length(x) == 0L) {
        if (x_ndim >= 2L && dims != 0L)
            x <- aperm(x)
        ans <- .colStats_SparseArray(op, x, na.rm=na.rm, center=center,
                                     dims=x_ndim-dims, useNames=useNames)
        if (!is.null(dim(ans)))
            ans <- aperm(ans)
        return(ans)
    }

    extract_j_slice <- function(j) {
        index <- vector("list", x_ndim)
        index[[2L]] <- j
        slice <- subset_SVT_by_Nindex(x, index, ignore.dimnames=TRUE)
        dim(slice) <- dim(x)[-2L]  # 'x_ndim - 1' dimensions
        slice
    }

    if (dims == 1L) {
        x <- aperm(x, perm=c(2:1, 3:x_ndim))
        ## We summarize the individual slices obtained by walking along the
        ## 2nd dimension of 'x'. Each slice has 'x_ndim - 1' dimensions.
        ans <- sapply(seq_len(x@dim[[2L]]),
            function(j) {
                slice <- extract_j_slice(j)  # 'x_ndim - 1' dimensions
                .colStats_SparseArray(op, slice, na.rm=na.rm, center=center,
                                      dims=x_ndim-1L, useNames=FALSE)
            })
        if (useNames)
            names(ans) <- x@dimnames[[2L]]
        return(ans)
    }

    if (dims == 2L) {
        ## We summarize the individual slices obtained by walking along the
        ## 2nd dimension of 'x'. Each slice has 'x_ndim - 1' dimensions.
        ans_cols <- lapply(seq_len(x@dim[[2L]]),
            function(j) {
                slice <- extract_j_slice(j)  # 'x_ndim - 1' dimensions
                .rowStats_SparseArray(op, slice, na.rm=na.rm, center=center,
                                      useNames=FALSE)
            })
        ans <- do.call(cbind, ans_cols)
        if (useNames)
            dimnames(ans) <- dimnames(x)[1:2]
        return(ans)
    }

    stop(wmsg("row*(<SVT_SparseArray>) summarizations ",
              "don't support 'dims' >= 3 yet"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colCountNAs/rowCountNAs
###
### Not part of the matrixStats API!

.colCountNAs_SVT_SparseMatrix <- function(x)
{
    if (is.null(x@SVT)) {
        ans <- integer(ncol(x))
    } else {
        ans <- vapply(x@SVT,
            function(lv) {
                if (is.null(lv))
                    return(0L)
                lv_vals <- lv[[2L]]
                sum(is.na(lv_vals))
            }, integer(1), USE.NAMES=FALSE)
    }
    setNames(ans, colnames(x))
}
#setMethod("colCountNAs", "SVT_SparseMatrix", .colCountNAs_SVT_SparseMatrix)

.rowCountNAs_SVT_SparseMatrix <- function(x)
{
    .colCountNAs_SVT_SparseMatrix(t(x))
}
#setMethod("rowCountNAs", "SVT_SparseMatrix", .rowCountNAs_SVT_SparseMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colAnyNAs/rowAnyNAs
###

.colAnyNAs_SparseArray <-
    function(x, rows=NULL, cols=NULL, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colAnyNAs")
    .colStats_SparseArray("anyNA", x, dims=dims, useNames=useNames)
}
setMethod("colAnyNAs", "SparseArray", .colAnyNAs_SparseArray)

.rowAnyNAs_SparseArray <-
    function(x, rows=NULL, cols=NULL, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "rowAnyNAs")
    .rowStats_SparseArray("anyNA", x, dims=dims, useNames=useNames)
}
setMethod("rowAnyNAs", "SparseArray", .rowAnyNAs_SparseArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colAnys/rowAnys and colAlls/rowAlls
###

.colAnys_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colAnys")
    .colStats_SparseArray("any", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colAnys", "SparseArray", .colAnys_SparseArray)

.rowAnys_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "rowAnys")
    .rowStats_SparseArray("any", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("rowAnys", "SparseArray", .rowAnys_SparseArray)

.colAlls_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colAlls")
    .colStats_SparseArray("all", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colAlls", "SparseArray", .colAlls_SparseArray)

.rowAlls_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "rowAlls")
    .rowStats_SparseArray("all", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("rowAlls", "SparseArray", .rowAlls_SparseArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colMins/rowMins, colMaxs/rowMaxs, and colRanges/rowRanges
###

.colMins_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colMins")
    .colStats_SparseArray("min", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colMins", "SparseArray", .colMins_SparseArray)

.rowMins_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "rowMins")
    .rowStats_SparseArray("min", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("rowMins", "SparseArray", .rowMins_SparseArray)

.colMaxs_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colMaxs")
    .colStats_SparseArray("max", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colMaxs", "SparseArray", .colMaxs_SparseArray)

.rowMaxs_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "rowMaxs")
    .rowStats_SparseArray("max", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("rowMaxs", "SparseArray", .rowMaxs_SparseArray)

.bind_mins_maxs <- function(mins, maxs, just.use.c)
{
    ## Bind 'mins' and 'maxs' together.
    if (just.use.c)
        return(c(mins, maxs))
    if (is.null(dim(mins))) {
        ans <- cbind(mins, maxs, deparse.level=0L)
        dimnames(ans) <- S4Arrays:::simplify_NULL_dimnames(dimnames(ans))
        return(ans)
    }
    ans_dimnames <- dimnames(mins)
    dim(mins) <- c(dim(mins), 1L)
    dim(maxs) <- c(dim(maxs), 1L)
    ans <- S4Arrays:::simple_abind(mins, maxs, along=length(dim(mins)))
    S4Arrays:::set_dimnames(ans, ans_dimnames)
}

.colRanges_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colRanges")
    ## Using two passes at the moment and binding the two results in R.
    ## TODO: Do all this in a single pass by calling
    ## '.colStats_SparseArray("range", ...)' and modifying .Call ENTRY POINT
    ## C_colStats_SVT to perform the binding from the very start at the C level.
    mins <- .colStats_SparseArray("min", x, na.rm=na.rm, dims=dims,
                                         useNames=useNames)
    maxs <- .colStats_SparseArray("max", x, na.rm=na.rm, dims=dims,
                                         useNames=FALSE)
    .bind_mins_maxs(mins, maxs, dims == length(dim(x)))
}
setMethod("colRanges", "SparseArray", .colRanges_SparseArray)

.rowRanges_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "rowRanges")
    ## Using two passes at the moment and binding the two results in R.
    ## TODO: Do all this in a single pass by calling
    ## '.rowStats_SparseArray("range", ...)' and modifying .Call ENTRY POINT
    ## C_colStats_SVT to perform the binding from the very start at the C level.
    mins <- .rowStats_SparseArray("min", x, na.rm=na.rm, dims=dims,
                                         useNames=useNames)
    maxs <- .rowStats_SparseArray("max", x, na.rm=na.rm, dims=dims,
                                         useNames=FALSE)
    .bind_mins_maxs(mins, maxs, dims == 0L)
}
setMethod("rowRanges", "SparseArray", .rowRanges_SparseArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colSums/rowSums, colProds/rowProds, and colMeans/rowMeans
###
### The colSums/rowSums/colMeans/rowMeans functions in base R propagate the
### dimnames so we do the same.

.colSums_SparseArray <- function(x, na.rm=FALSE, dims=1)
{
    .colStats_SparseArray("sum", x, na.rm=na.rm, dims=dims)
}
setMethod("colSums", "SparseArray", .colSums_SparseArray)

.rowSums_SparseArray <- function(x, na.rm=FALSE, dims=1)
{
    .rowStats_SparseArray("sum", x, na.rm=na.rm, dims=dims)
}
setMethod("rowSums", "SparseArray", .rowSums_SparseArray)

.colProds_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colProds")
    .colStats_SparseArray("prod", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colProds", "SparseArray", .colProds_SparseArray)

.rowProds_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "rowProds")
    .rowStats_SparseArray("prod", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("rowProds", "SparseArray", .rowProds_SparseArray)

.colMeans_SparseArray <- function(x, na.rm=FALSE, dims=1)
{
    .colStats_SparseArray("mean", x, na.rm=na.rm, dims=dims)
}
setMethod("colMeans", "SparseArray", .colMeans_SparseArray)

.rowMeans_SparseArray <- function(x, na.rm=FALSE, dims=1)
{
    .rowStats_SparseArray("mean", x, na.rm=na.rm, dims=dims)
}
setMethod("rowMeans", "SparseArray", .rowMeans_SparseArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colSums2/rowSums2 and colMeans2/rowMeans2
###

.colSums2_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colSums2")
    .colStats_SparseArray("sum", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colSums2", "SparseArray", .colSums2_SparseArray)

.rowSums2_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "rowSums2")
    .rowStats_SparseArray("sum", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("rowSums2", "SparseArray", .rowSums2_SparseArray)

.colMeans2_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colMeans2")
    .colStats_SparseArray("mean", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colMeans2", "SparseArray", .colMeans2_SparseArray)

.rowMeans2_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "rowMeans2")
    .rowStats_SparseArray("mean", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("rowMeans2", "SparseArray", .rowMeans2_SparseArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colVars/rowVars and colSds/rowSds
###

### Equivalent to 'var(c(x, integer(padding)), ...)' but doesn't actually
### realize the padding with zeros.
.padded_var <- function(x, padding=0L, na.rm=FALSE, center=NULL)
{
    if (na.rm)
        x <- x[!is.na(x)]
    nvals <- length(x) + padding
    if (nvals <= 1L)
        return(NA_real_)
    if (is.null(center)) {
        center <- sum(x) / nvals
    } else {
        stopifnot(isSingleNumberOrNA(center))
    }
    delta <- x - center
    s <- sum(delta * delta) + center * center * padding
    s / (nvals - 1L)
}

### Returns a numeric vector of length 'ncol(x)'.
.normarg_center <- function(center, x, na.rm=FALSE)
{
    if (is.null(center))
        return(colMeans(x, na.rm=na.rm))
    if (!is.numeric(center))
        stop(wmsg("'center' must be NULL or a numeric vector"))
    x_ncol <- ncol(x)
    if (length(center) != x_ncol) {
        if (length(center) != 1L)
            stop(wmsg("'center' must have one element per row ",
                      "or column in the SparseMatrix object"))
        center <- rep.int(center, x_ncol)
    }
    center
}

### Original "pure R" implementation. Was originally used by the colVars()
### method for SVT_SparseMatrix objects. No longer used!
.colVars_SparseArray_SparseMatrix <-
    function(x, na.rm=FALSE, center=NULL, useNames=NA)
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    useNames <- .normarg_useNames(useNames)
    x_nrow <- nrow(x)
    x_ncol <- ncol(x)
    if (x_nrow <= 1L) {
        ans <- rep.int(NA_real_, x_ncol)
    } else {
        center <- .normarg_center(center, x, na.rm=na.rm)
        ans <- center * center * x_nrow / (x_nrow - 1L)
        if (!is.null(x@SVT)) {
            ans <- vapply(seq_along(x@SVT),
                function(i) {
                    lv <- x@SVT[[i]]
                    if (is.null(lv))
                        return(ans[[i]])
                    lv_vals <- lv[[2L]]
                    padding <- x_nrow - length(lv_vals)
                    .padded_var(lv_vals, padding, na.rm=na.rm,
                                center=center[[i]])
                }, numeric(1), USE.NAMES=FALSE)
        }
    }
    if (useNames)
        names(ans) <- colnames(x)
    ans
}

.colVars_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL,
                dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colVars")
    .colStats_SparseArray("var1", x, na.rm=na.rm, center=center,
                                  dims=dims, useNames=useNames)
}
setMethod("colVars", "SparseArray", .colVars_SparseArray)

.rowVars_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL,
                dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "rowVars")
    .rowStats_SparseArray("var1", x, na.rm=na.rm, center=center,
                                  dims=dims, useNames=useNames)
}
setMethod("rowVars", "SparseArray", .rowVars_SparseArray)

.colSds_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL,
                dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colSds")
    .colStats_SparseArray("sd1", x, na.rm=na.rm, center=center,
                                 dims=dims, useNames=useNames)
}
setMethod("colSds", "SparseArray", .colSds_SparseArray)

.rowSds_SparseArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL,
                dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "rowSds")
    .rowStats_SparseArray("sd1", x, na.rm=na.rm, center=center,
                                 dims=dims, useNames=useNames)
}
setMethod("rowSds", "SparseArray", .rowSds_SparseArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colMedians/rowMedians
###
### TODO: How hard would it be to replace current "pure R" implementation
### with C implementation available thru .Call ENTRY POINT C_colStats_SVT ?

### All values in 'x' are **assumed** to be >= 0 but we don't check this!
### 'padding' is expected to be < length(x).
.positive_padded_median <- function(x, padding=0L)
{
    x_len <- length(x)
    stopifnot(padding < x_len)
    n <- x_len + padding
    if (n %% 2L == 1L) {
        middle <- (n + 1L) %/% 2L
        partial <- middle - padding
        return(sort(x, partial=partial)[partial])
    }
    i1 <- n %/% 2L - padding
    i2 <- i1 + 1L
    mean(sort(x, partial=i2)[i1:i2])
}

### Equivalent to 'median(c(x, integer(padding)), ...)' but doesn't actually
### realize the padding with zeros.
.padded_median <- function(x, padding=0L, na.rm=FALSE)
{
    if (na.rm) {
        x <- x[!is.na(x)]
    } else {
        if (anyNA(x))
            return(NA_real_)
    }
    n <- length(x) + padding
    if (n == 0L)
        return(NA_real_)
    if (padding > length(x))
        return(0)

    ## Handle case where we have more positive values than non-positive values.
    pos_idx <- which(x > 0L)
    pos_count <- length(pos_idx)
    nonpos_count <- n - pos_count
    if (pos_count > nonpos_count) {
        ans <- .positive_padded_median(x[pos_idx], padding=nonpos_count)
        return(ans)
    }

    ## Handle case where we have more negative values than non-negative values.
    neg_count <- length(x) - pos_count
    nonneg_count <- n - neg_count
    if (neg_count > nonneg_count) {
        ans <- - .positive_padded_median(-x[-pos_idx], padding=nonneg_count)
        return(ans)
    }

    if (n %% 2L == 1L)
        return(0)

    half <- n %/% 2L
    if (pos_count == half) {
        right <- min(x[pos_idx])
    } else {
        right <- 0
    }
    if (neg_count == half) {
        left <- max(x[-pos_idx])
    } else {
        left <- 0
    }
    (left + right) * 0.5
}

.colMedians_SVT_SparseMatrix <- function(x, na.rm=FALSE, useNames=NA)
{
    stopifnot(is(x, "SVT_SparseMatrix"))
    check_svt_version(x)
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    useNames <- .normarg_useNames(useNames)
    x_nrow <- nrow(x)
    x_ncol <- ncol(x)
    if (x_nrow == 0L) {
        ans <- rep.int(NA_real_, x_ncol)
    } else {
        ans <- numeric(x_ncol)
        if (!is.null(x@SVT)) {
            ans <- vapply(seq_along(x@SVT),
                function(i) {
                    lv <- x@SVT[[i]]
                    if (is.null(lv))
                        return(ans[[i]])
                    lv_vals <- lv[[2L]]
                    padding <- x_nrow - length(lv_vals)
                    .padded_median(lv_vals, padding, na.rm=na.rm)
                }, numeric(1), USE.NAMES=FALSE)
        }
    }
    if (useNames)
        names(ans) <- colnames(x)
    ans
}

setMethod("colMedians", "SparseArray",
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
    {
        .check_unused_arguments(...)
        .stopifnot_2D_object(x, "colMedians")
        .check_rows_cols(rows, cols, "colMedians")
        if (!is(x, "SVT_SparseArray"))
            x <- as(x, "SVT_SparseArray")
        .colMedians_SVT_SparseMatrix(x, na.rm=na.rm, useNames=useNames)
    }
)

setMethod("rowMedians", "SparseArray",
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
    {
        .check_unused_arguments(...)
        .stopifnot_2D_object(x, "rowMedians")
        .check_rows_cols(rows, cols, "rowMedians")
        tx <- t(x)
        if (!is(tx, "SVT_SparseArray"))
            x <- as(tx, "SVT_SparseArray")
        .colMedians_SVT_SparseMatrix(tx, na.rm=na.rm, useNames=useNames, ...)
    }
)

