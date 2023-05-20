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
### - All other matrix row/column summarization operations are from the
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

.check_useNames <- function(useNames)
{
    if (!(is.logical(useNames) && length(useNames) == 1L))
        stop(wmsg("'useNames' must be a single logical value"))
}

.stopifnot_2D_object <- function(x, method)
{
    if (length(dim(x)) != 2L)
        stop(wmsg("the ", method, "() method for SVT_SparseArray objects ",
                  "only supports 2D objects (i.e. SVT_SparseMatrix objects) ",
                  "at the moment"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Workhorse behind all the matrixStats methods
###

### Return an ordinary array with 'length(dim(x)) - dims' dimensions.
.colStats_SVT <- function(op, x, na.rm=FALSE, center=NULL, dims=1L, useNames=NA)
{
    stopifnot(isSingleString(op), is(x, "SVT_SparseArray"))

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
    if (!isSingleNumber(dims))
        stop(wmsg("'dims' must be a single integer"))
    if (!is.integer(dims))
        dims <- as.integer(dims)

    ## Check 'useNames'.
    .check_useNames(useNames)

    x_dimnames <- if (isFALSE(useNames)) NULL else x@dimnames
    .Call2("C_colStats_SVT", x@dim, x_dimnames, x@type, x@SVT,
                             op, na.rm, center, dims,
                             PACKAGE="SparseArray")
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
### colAnys/rowAnys and colAlls/rowAlls
###

.colAnys_SVT <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colAnys")
    .colStats_SVT("any", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colAnys", "SVT_SparseArray", .colAnys_SVT)

.rowAnys_SVT <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .stopifnot_2D_object(x, "rowAnys")
    .check_rows_cols(rows, cols, "rowAnys")
    .colAnys_SVT(t(x), na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("rowAnys", "SVT_SparseArray", .rowAnys_SVT)

.colAlls_SVT <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colAlls")
    .colStats_SVT("all", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colAlls", "SVT_SparseArray", .colAlls_SVT)

.rowAlls_SVT <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .stopifnot_2D_object(x, "rowAlls")
    .check_rows_cols(rows, cols, "rowAlls")
    .colAlls_SVT(t(x), na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("rowAlls", "SVT_SparseArray", .rowAlls_SVT)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colMins/rowMins, colMaxs/rowMaxs, and colRanges/rowRanges
###

### Original "pure R" implementation. Was originally used by the colMins()
### and colMaxs() methods for SVT_SparseMatrix objects.
### No longer used!
.colMinsMaxs_SVT_SparseMatrix <- function(x, FUN, na.rm=FALSE, useNames=NA)
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    .check_useNames(useNames)
    x_nrow <- nrow(x)
    x_ncol <- ncol(x)
    if (x_ncol == 0L) {
        ans <- vector(type(x), length=0L)
    } else if (x_nrow == 0L) {
        ans <- rep.int(suppressWarnings(FUN()), x_ncol)
    } else if (is.null(x@SVT)) {
        ans <- vector(type(x), length=x_ncol)
    } else {
        zero <- vector(type(x), length=1L)
        ## We use lapply() rather than vapply() because we don't know in
        ## advance the type of the result we're going to get for each column.
        ans <- unlist(lapply(x@SVT,
            function(lv) {
                if (is.null(lv))
                    return(zero)
                lv_vals <- lv[[2L]]
                ## Suppress the warning that min() or max() will issue
                ## if 'na.rm' is TRUE and 'lv_vals' contains no non-NA values.
                res <- suppressWarnings(FUN(lv_vals, na.rm=na.rm))
                if (length(lv_vals) < x_nrow)
                    res <- FUN(res, zero)
                res
            }), use.names=FALSE)
    }
    if (!isFALSE(useNames))
        names(ans) <- colnames(x)
    ans
}

.colMins_SVT <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colMins")
    .colStats_SVT("min", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colMins", "SVT_SparseArray", .colMins_SVT)

.rowMins_SVT <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .stopifnot_2D_object(x, "rowMins")
    .check_rows_cols(rows, cols, "rowMins")
    .colMins_SVT(t(x), na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("rowMins", "SVT_SparseArray", .rowMins_SVT)

.colMaxs_SVT <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colMaxs")
    .colStats_SVT("max", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colMaxs", "SVT_SparseArray", .colMaxs_SVT)

.rowMaxs_SVT <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .stopifnot_2D_object(x, "rowMaxs")
    .check_rows_cols(rows, cols, "rowMaxs")
    .colMaxs_SVT(t(x), na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("rowMaxs", "SVT_SparseArray", .rowMaxs_SVT)

.colRanges_SVT <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colRanges")

    ## Using two passes at the moment and binding the two results in R.
    ## TODO: Do all this in a single pass. Call '.colStats_SVT("range", ...)'
    ## and modify .Call ENTRY POINT C_colStats_SVT to perform the binding
    ## from the very start.
    mins <- .colStats_SVT("min", x, na.rm=na.rm, dims=dims, useNames=useNames)
    maxs <- .colStats_SVT("max", x, na.rm=na.rm, dims=dims, useNames=FALSE)

    ## Bind 'mins' and 'maxs' together.
    if (dims == length(dim(x)))
        return(c(mins, maxs))
    if (is.null(dim(mins)))
        return(cbind(mins, maxs, deparse.level=0L))
    ans_dimnames <- dimnames(mins)
    dim(mins) <- c(dim(mins), 1L)
    dim(maxs) <- c(dim(maxs), 1L)
    ans <- S4Arrays:::simple_abind(mins, maxs, along=length(dim(mins)))
    S4Arrays:::set_dimnames(ans, ans_dimnames)
}
setMethod("colRanges", "SVT_SparseArray", .colRanges_SVT)

.rowRanges_SVT <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .stopifnot_2D_object(x, "rowRanges")
    .check_rows_cols(rows, cols, "rowRanges")
    .colRanges_SVT(t(x), na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("rowRanges", "SVT_SparseArray", .rowRanges_SVT)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colSums/rowSums, colProds/rowProds, and colMeans/rowMeans
###
### The colSums/rowSums/colMeans/rowMeans functions in base R propagate the
### dimnames so we do the same.

.colSums_SVT <- function(x, na.rm=FALSE, dims=1)
{
    .colStats_SVT("sum", x, na.rm=na.rm, dims=dims)
}
setMethod("colSums", "SVT_SparseArray", .colSums_SVT)

.rowSums_SVT <- function(x, na.rm=FALSE, dims=1)
{
    .stopifnot_2D_object(x, "rowSums")
    .colSums_SVT(t(x), na.rm=na.rm, dims=dims)
}
setMethod("rowSums", "SVT_SparseArray", .rowSums_SVT)

.colProds_SVT <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colProds")
    .colStats_SVT("prod", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colProds", "SVT_SparseArray", .colProds_SVT)

.rowProds_SVT <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .stopifnot_2D_object(x, "rowProds")
    .check_rows_cols(rows, cols, "rowProds")
    .colProds_SVT(t(x), na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("rowProds", "SVT_SparseArray", .rowProds_SVT)

.colMeans_SVT <- function(x, na.rm=FALSE, dims=1)
{
    .colStats_SVT("mean", x, na.rm=na.rm, dims=dims)
}
setMethod("colMeans", "SVT_SparseArray", .colMeans_SVT)

.rowMeans_SVT <- function(x, na.rm=FALSE, dims=1)
{
    .stopifnot_2D_object(x, "rowMeans")
    .colMeans_SVT(t(x), na.rm=na.rm, dims=dims)
}
setMethod("rowMeans", "SVT_SparseArray", .rowMeans_SVT)


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
.colVars_SVT_SparseMatrix <- function(x, na.rm=FALSE, center=NULL, useNames=NA)
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    .check_useNames(useNames)
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
    if (isTRUE(useNames))
        names(ans) <- colnames(x)
    ans
}

.colVars_SVT <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL,
                            dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colVars")
    .colStats_SVT("var1", x, na.rm=na.rm, center=center,
                             dims=dims, useNames=useNames)
}
setMethod("colVars", "SVT_SparseArray", .colVars_SVT)

.rowVars_SVT <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL,
                            dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .stopifnot_2D_object(x, "rowVars")
    .check_rows_cols(rows, cols, "rowVars")
    .colVars_SVT(t(x), na.rm=na.rm, center=center,
                       dims=dims, useNames=useNames)
}
setMethod("rowVars", "SVT_SparseArray", .rowVars_SVT)

.colSds_SVT <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL,
                           dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .check_rows_cols(rows, cols, "colSds")
    .colStats_SVT("sd1", x, na.rm=na.rm, center=center,
                            dims=dims, useNames=useNames)
}
setMethod("colSds", "SVT_SparseArray", .colSds_SVT)

.rowSds_SVT <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL,
                           dims=1, ..., useNames=NA)
{
    .check_unused_arguments(...)
    .stopifnot_2D_object(x, "rowSds")
    .check_rows_cols(rows, cols, "rowSds")
    .colSds_SVT(t(x), na.rm=na.rm, center=center,
                      dims=dims, useNames=useNames)
}
setMethod("rowSds", "SVT_SparseArray", .rowSds_SVT)


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
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    .check_useNames(useNames)
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
    if (isTRUE(useNames))
        names(ans) <- colnames(x)
    ans
}

setMethod("colMedians", "SVT_SparseArray",
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
    {
        .check_unused_arguments(...)
        .stopifnot_2D_object(x, "colMedians")
        .check_rows_cols(rows, cols, "colMedians")
        .colMedians_SVT_SparseMatrix(x, na.rm=na.rm, useNames=useNames)
    }
)

setMethod("rowMedians", "SVT_SparseArray",
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
    {
        .check_unused_arguments(...)
        .stopifnot_2D_object(x, "rowMedians")
        .check_rows_cols(rows, cols, "rowMedians")
        .colMedians_SVT_SparseMatrix(t(x), na.rm=na.rm, useNames=useNames, ...)
    }
)

