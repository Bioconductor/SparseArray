### =========================================================================
### matrixStats methods for SparseMatrix objects
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
### - The colSums() and colMeans() generics are defined in the BiocGenerics
###   package.
### - All other matrix row/column summarization operations are from the
###   matrixStats package with corresponding generics defined in the
###   MatrixGenerics package.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers
###

.check_dims <- function(dims, method)
{
    if (!identical(dims, 1))
        stop(wmsg("\"", method, "\" method for SVT_SparseMatrix ",
                  "objects does not support the 'dims' argument"))
}

.check_rows_cols <- function(rows, cols, method)
{
    if (!(is.null(rows) && is.null(cols)))
        stop(wmsg("\"", method, "\" method for SparseMatrix objects ",
                  "does not support arguments 'rows' and 'cols'"))
}

.check_useNames <- function(useNames)
{
    if (!(is.logical(useNames) && length(useNames) == 1L))
        stop(wmsg("'useNames' must be a single logical value"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colSums/rowSums
###
### Like base::colSums() and base::rowSums(), the colSums() and rowSums()
### generics in BiocGenerics have a 'dims' argument. We do NOT support it.
### base::colSums() and base::rowSums() propagate the colnames and rownames.

.colSums_SVT_SparseMatrix <- function(x, na.rm=FALSE, dims=1)
{
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    .check_dims(dims, "colSums")
    if (is.null(x@SVT)) {
        ans <- numeric(ncol(x))
    } else {
        ans <- vapply(x@SVT,
            function(lv) {
                if (is.null(lv))
                    return(0)
                lv_vals <- lv[[2L]]
                sum(lv_vals, na.rm=na.rm)
            }, numeric(1), USE.NAMES=FALSE)
    }
    setNames(ans, colnames(x))
}
setMethod("colSums", "SVT_SparseMatrix", .colSums_SVT_SparseMatrix)

.rowSums_SVT_SparseMatrix <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "rowSums")
    .colSums_SVT_SparseMatrix(t(x), na.rm=na.rm)
}
setMethod("rowSums", "SVT_SparseMatrix", .rowSums_SVT_SparseMatrix)


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
### colMeans/rowMeans
###
### Like base::colMeans() and base::rowMeans(), the colMeans() and rowMeans()
### generics in BiocGenerics have a 'dims' argument. We do NOT support it.
### base::colMeans() and base::rowMeans() propagate the colnames and rownames.

.colMeans_SVT_SparseMatrix <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "colMeans")
    sums <- .colSums_SVT_SparseMatrix(x, na.rm=na.rm)
    nvals <- rep.int(nrow(x), ncol(x))
    if (na.rm)
        nvals <- nvals - .colCountNAs_SVT_SparseMatrix(x)
    sums / nvals
}
setMethod("colMeans", "SVT_SparseMatrix", .colMeans_SVT_SparseMatrix)

.rowMeans_SVT_SparseMatrix <- function(x, na.rm=FALSE, dims=1)
{
    .check_dims(dims, "rowMeans")
    .colMeans_SVT_SparseMatrix(t(x), na.rm=na.rm)
}
setMethod("rowMeans", "SVT_SparseMatrix", .rowMeans_SVT_SparseMatrix)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colMins/rowMins, colMaxs/rowMaxs, and colRanges/rowRanges
###

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
    if (isTRUE(useNames))
        names(ans) <- colnames(x)
    ans
}

setMethod("colMins", "SVT_SparseMatrix",
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
    {
        .check_rows_cols(rows, cols, "colMins")
        .colMinsMaxs_SVT_SparseMatrix(x, base::min, na.rm=na.rm,
                                      useNames=useNames, ...)
    }
)

setMethod("rowMins", "SVT_SparseMatrix",
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
    {
        .check_rows_cols(rows, cols, "rowMins")
        .colMinsMaxs_SVT_SparseMatrix(t(x), base::min, na.rm=na.rm,
                                      useNames=useNames, ...)
    }
)

setMethod("colMaxs", "SVT_SparseMatrix",
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
    {
        .check_rows_cols(rows, cols, "colMaxs")
        .colMinsMaxs_SVT_SparseMatrix(x, base::max, na.rm=na.rm,
                                      useNames=useNames, ...)
    }
)

setMethod("rowMaxs", "SVT_SparseMatrix",
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
    {
        .check_rows_cols(rows, cols, "rowMaxs")
        .colMinsMaxs_SVT_SparseMatrix(t(x), base::max, na.rm=na.rm,
                                      useNames=useNames, ...)
    }
)

.colRanges_SVT_SparseMatrix <- function(x, na.rm=FALSE, useNames=NA)
{
    mins <- .colMinsMaxs_SVT_SparseMatrix(x, base::min, na.rm=na.rm,
                                          useNames=useNames)
    maxs <- .colMinsMaxs_SVT_SparseMatrix(x, base::max, na.rm=na.rm,
                                          useNames=FALSE)
    cbind(mins, maxs, deparse.level=0L)
}

setMethod("colRanges", "SVT_SparseMatrix",
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
    {
        .check_rows_cols(rows, cols, "colRanges")
        .colRanges_SVT_SparseMatrix(x, na.rm=na.rm, useNames=useNames)
    }
)

setMethod("rowRanges", "SVT_SparseMatrix",
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
    {
        .check_rows_cols(rows, cols, "rowRanges")
        .colRanges_SVT_SparseMatrix(t(x), na.rm=na.rm, useNames=useNames)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colMedians/rowMedians
###

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

setMethod("colMedians", "SVT_SparseMatrix",
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
    {
        .check_rows_cols(rows, cols, "colMedians")
        .colMedians_SVT_SparseMatrix(x, na.rm=na.rm, useNames=useNames, ...)
    }
)

setMethod("rowMedians", "SVT_SparseMatrix",
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
    {
        .check_rows_cols(rows, cols, "rowMedians")
        .colMedians_SVT_SparseMatrix(t(x), na.rm=na.rm, useNames=useNames, ...)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colVars/rowVars
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

setMethod("colVars", "SVT_SparseMatrix",
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL,
             ..., useNames=NA)
    {
        .check_rows_cols(rows, cols, "colVars")
        .colVars_SVT_SparseMatrix(x, na.rm=na.rm, center=center,
                                  useNames=useNames, ...)
    }
)

setMethod("rowVars", "SVT_SparseMatrix",
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL,
             ..., useNames=NA)
    {
        .check_rows_cols(rows, cols, "rowVars")
        .colVars_SVT_SparseMatrix(t(x), na.rm=na.rm, center=center,
                                  useNames=useNames, ...)
    }
)

