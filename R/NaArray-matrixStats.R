### =========================================================================
### matrixStats methods for NaMatrix and NaArray objects
### -------------------------------------------------------------------------
###
### See notes at beginning of NaArray-matrixStats.R for matrixStats usage
### in Bioconductor.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .colStats_NaArray(), .rowStats_NaArray()
###
### Workhorses behind all the matrixStats methods for NaArray objects, with
### the exception of the colMedians()/rowMedians() methods at the moment.
###

### Returns an ordinary array with 'length(dim(x)) - dims' dimensions.
.colStats_NaArray <- function(op, x, na.rm=FALSE, center=NULL, dims=1L,
                              useNames=NA)
{
    stopifnot(isSingleString(op), is(x, "NaArray"))
    check_svt_version(x)

    ## Normalize and check 'dims'.
    dims <- normarg_dims(dims)
    if (dims <= 0L || dims > length(x@dim))
        stop(wmsg("'dims' must be a single integer that is ",
                  "> 0 and <= length(dim(x)) for the col*() functions, and ",
                  ">= 0 and < length(dim(x)) for the row*() functions"))

    ## Check 'na.rm'.
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))

    ## Check and normalize 'center'.
    if (is.null(center)) {
        center <- NA_real_
    } else {
        if (!isSingleNumberOrNA(center))
            stop(wmsg("'center' must be NULL or a single number"))
        if (!is.double(center))
            center <- as.double(center)
    }

    ## Normalize 'useNames'.
    useNames <- normarg_useNames(useNames)

    x_dimnames <- if (useNames) x@dimnames else NULL
    SparseArray.Call("C_colStats_SVT",
                     x@dim, x_dimnames, x@type, x@NaSVT, TRUE,
                     op, na.rm, center, dims)
}

### Returns an ordinary array where the number of dimensions is 'dims'.
.rowStats_NaArray <- function(op, x, na.rm=FALSE, center=NULL, dims=1L,
                              useNames=NA)
{
    stop("not ready yet")

    stopifnot(isSingleString(op), is(x, "NaArray"))
    check_svt_version(x)

    ## Normalize and check 'dims'.
    dims <- normarg_dims(dims)
    if (dims < 0L || dims >= length(x@dim))
        stop(wmsg("'dims' must be a single integer that is ",
                  "> 0 and <= length(dim(x)) for the col*() functions, and ",
                  ">= 0 and < length(dim(x)) for the row*() functions"))

    if (dims == 0L)
        return(.colStats_NaArray(op, x, na.rm=na.rm, center=center,
                                 dims=length(x@dim), useNames=useNames))

    ## Check 'na.rm'.
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))

    ## Check and normalize 'center'.
    if (!is.null(center)) {
        ## Unlike for .colStats_NaArray() where 'center' can only be NULL
        ## or a single number, here it can also be an ordinary numeric array
        ## of the same dimensions as the result of .rowStats_NaArray()
        ## (i.e. of dimensions 'head(dim(x), n=dims)'), or a numeric vector
        ## of the same length as the result of .rowStats_NaArray().
        if (!is.numeric(center))
            stop(wmsg("'center' must be NULL, a single number, ",
                      "or an ordinary array"))
        ans_dim <- head(dim(x), n=dims)
        if (is.array(center)) {
            if (!identical(dim(center), ans_dim))
                stop(wmsg("unexpected 'center' dimensions"))
            if (storage.mode(center) != "double")
                storage.mode(center) <- "double"
        } else if (length(center) %in% c(1L, prod(ans_dim))) {
            center <- array(as.double(center), dim=ans_dim)
        } else {
            stop(wmsg("unexpected 'center' length"))
        }
    }

    ## Normalize 'useNames'.
    useNames <- normarg_useNames(useNames)

    x_dimnames <- if (useNames) x@dimnames else NULL
    SparseArray.Call("C_rowStats_SVT",
                     x@dim, x_dimnames, x@type, x@NaSVT, TRUE,
                     op, na.rm, center, dims)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colCountNAs/rowCountNAs
###
### Not part of the matrixStats API!

.colCountNAs_NaArray <- function(x, dims=1, useNames=NA)
{
    .colStats_NaArray("countNAs", x, dims=dims, useNames=useNames)
}
#setMethod("colCountNAs", "NaArray", .colCountNAs_NaArray)

.rowCountNAs_NaArray <- function(x, dims=1, useNames=NA)
{
    .rowStats_NaArray("countNAs", x, dims=dims, useNames=useNames)
}
#setMethod("rowCountNAs", "NaArray", .rowCountNAs_NaArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .colCountVals_NaArray()/.rowCountVals_NaArray()
###
### Count the number of non-NA vals per column/row.
###
### Both functions return a single value if 'na.rm' is FALSE, or an unnamed
### ordinary vector, matrix, or array if 'na.rm' is TRUE.

.colCountVals_NaArray <- function(x, na.rm=FALSE, dims=1)
{
    stopifnot(is(x, "NaArray"))
    dims <- normarg_dims(dims)
    ans <- prod(head(dim(x), n=dims))
    if (na.rm) {
        count_nas <- .colCountNAs_NaArray(x, dims=dims, useNames=FALSE)
        ans <- ans - count_nas
    }
    ans
}

.rowCountVals_NaArray <- function(x, na.rm=FALSE, dims=1)
{
    stopifnot(is(x, "NaArray"))
    dims <- normarg_dims(dims)
    ans <- prod(tail(dim(x), n=-dims))
    if (na.rm) {
        count_nas <- .rowCountNAs_NaArray(x, dims=dims, useNames=FALSE)
        ans <- ans - count_nas
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colAnyNAs/rowAnyNAs
###

.colAnyNAs_NaArray <-
    function(x, rows=NULL, cols=NULL, dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "colAnyNAs", "NaArray")
    .colStats_NaArray("anyNA", x, dims=dims, useNames=useNames)
}
setMethod("colAnyNAs", "NaArray", .colAnyNAs_NaArray)

.rowAnyNAs_NaArray <-
    function(x, rows=NULL, cols=NULL, dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "rowAnyNAs", "NaArray")
    .rowStats_NaArray("anyNA", x, dims=dims, useNames=useNames)
}
#setMethod("rowAnyNAs", "NaArray", .rowAnyNAs_NaArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colAnys/rowAnys and colAlls/rowAlls
###

.colAnys_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "colAnys", "NaArray")
    .colStats_NaArray("any", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colAnys", "NaArray", .colAnys_NaArray)

.rowAnys_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "rowAnys", "NaArray")
    .rowStats_NaArray("any", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
#setMethod("rowAnys", "NaArray", .rowAnys_NaArray)

.colAlls_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "colAlls", "NaArray")
    .colStats_NaArray("all", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colAlls", "NaArray", .colAlls_NaArray)

.rowAlls_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "rowAlls", "NaArray")
    .rowStats_NaArray("all", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
#setMethod("rowAlls", "NaArray", .rowAlls_NaArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colMins/rowMins, colMaxs/rowMaxs, and colRanges/rowRanges
###

.colMins_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "colMins", "NaArray")
    .colStats_NaArray("min", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colMins", "NaArray", .colMins_NaArray)

.rowMins_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "rowMins", "NaArray")
    .rowStats_NaArray("min", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
#setMethod("rowMins", "NaArray", .rowMins_NaArray)

.colMaxs_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "colMaxs", "NaArray")
    .colStats_NaArray("max", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colMaxs", "NaArray", .colMaxs_NaArray)

.rowMaxs_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "rowMaxs", "NaArray")
    .rowStats_NaArray("max", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
#setMethod("rowMaxs", "NaArray", .rowMaxs_NaArray)

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

.colRanges_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "colRanges", "NaArray")
    ## Using two passes at the moment and binding the two results in R.
    ## TODO: Do all this in a single pass by calling
    ## '.colStats_NaArray("range", ...)' and modifying .Call ENTRY POINT
    ## C_colStats_SVT to perform the binding from the very start at the C level.
    mins <- .colStats_NaArray("min", x, na.rm=na.rm, dims=dims,
                                         useNames=useNames)
    maxs <- .colStats_NaArray("max", x, na.rm=na.rm, dims=dims,
                                         useNames=FALSE)
    .bind_mins_maxs(mins, maxs, dims == length(dim(x)))
}
setMethod("colRanges", "NaArray", .colRanges_NaArray)

.rowRanges_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "rowRanges", "NaArray")
    ## Using two passes at the moment and binding the two results in R.
    ## TODO: Do all this in a single pass by calling
    ## '.rowStats_NaArray("range", ...)' and modifying .Call ENTRY POINT
    ## C_colStats_SVT to perform the binding from the very start at the C level.
    mins <- .rowStats_NaArray("min", x, na.rm=na.rm, dims=dims,
                                         useNames=useNames)
    maxs <- .rowStats_NaArray("max", x, na.rm=na.rm, dims=dims,
                                         useNames=FALSE)
    .bind_mins_maxs(mins, maxs, dims == 0L)
}
#setMethod("rowRanges", "NaArray", .rowRanges_NaArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colSums/rowSums, colProds/rowProds, and colMeans/rowMeans
###
### The colSums/rowSums/colMeans/rowMeans functions in base R propagate the
### dimnames so we do the same.

.colSums_NaArray <- function(x, na.rm=FALSE, dims=1)
{
    .colStats_NaArray("sum", x, na.rm=na.rm, dims=dims)
}
setMethod("colSums", "NaArray", .colSums_NaArray)

.rowSums_NaArray <- function(x, na.rm=FALSE, dims=1)
{
    .rowStats_NaArray("sum", x, na.rm=na.rm, dims=dims)
}
#setMethod("rowSums", "NaArray", .rowSums_NaArray)

.colProds_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "colProds", "NaArray")
    .colStats_NaArray("prod", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colProds", "NaArray", .colProds_NaArray)

.rowProds_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "rowProds", "NaArray")
    .rowStats_NaArray("prod", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
#setMethod("rowProds", "NaArray", .rowProds_NaArray)

.colMeans_NaArray <- function(x, na.rm=FALSE, dims=1)
{
    .colStats_NaArray("mean", x, na.rm=na.rm, dims=dims)
}
setMethod("colMeans", "NaArray", .colMeans_NaArray)

.rowMeans_NaArray <- function(x, na.rm=FALSE, dims=1)
{
    sums <- .rowSums_NaArray(x, na.rm=na.rm, dims=dims)
    nvals <- .rowCountVals_NaArray(x, na.rm=na.rm, dims=dims)
    sums / nvals
}
#setMethod("rowMeans", "NaArray", .rowMeans_NaArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colSums2/rowSums2 and colMeans2/rowMeans2
###

.colSums2_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "colSums2", "NaArray")
    .colStats_NaArray("sum", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colSums2", "NaArray", .colSums2_NaArray)

.rowSums2_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "rowSums2", "NaArray")
    .rowStats_NaArray("sum", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
#setMethod("rowSums2", "NaArray", .rowSums2_NaArray)

.colMeans2_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "colMeans2", "NaArray")
    .colStats_NaArray("mean", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
setMethod("colMeans2", "NaArray", .colMeans2_NaArray)

.rowMeans2_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "rowMeans2", "NaArray")
    .rowStats_NaArray("mean", x, na.rm=na.rm, dims=dims, useNames=useNames)
}
#setMethod("rowMeans2", "NaArray", .rowMeans2_NaArray)


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

.colVars_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL,
                dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "colVars", "NaArray")
    .colStats_NaArray("var1", x, na.rm=na.rm, center=center,
                                  dims=dims, useNames=useNames)
}
setMethod("colVars", "NaArray", .colVars_NaArray)

.rowVars_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL,
                dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "rowVars", "NaArray")
    nvals <- .rowCountVals_NaArray(x, na.rm=na.rm, dims=dims)
    if (is.null(center)) {
        sums <- .rowSums_NaArray(x, na.rm=na.rm, dims=dims)
        center <- sums / nvals
    }
    centered_X2_sums <- .rowStats_NaArray("centered_X2_sum",
                                              x, na.rm=na.rm, center=center,
                                              dims=dims, useNames=useNames)
    centered_X2_sums / (nvals - 1)
}
#setMethod("rowVars", "NaArray", .rowVars_NaArray)

.colSds_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL,
                dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "colSds", "NaArray")
    .colStats_NaArray("sd1", x, na.rm=na.rm, center=center,
                                 dims=dims, useNames=useNames)
}
setMethod("colSds", "NaArray", .colSds_NaArray)

.rowSds_NaArray <-
    function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL,
                dims=1, ..., useNames=NA)
{
    check_unused_arguments(...)
    check_rows_cols(rows, cols, "rowSds", "NaArray")
    row_vars <- .rowVars_NaArray(x, na.rm=na.rm, center=center,
                                     dims=dims, useNames=useNames)
    sqrt(row_vars)
}
#setMethod("rowSds", "NaArray", .rowSds_NaArray)


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

.colMedians_NaMatrix <- function(x, na.rm=FALSE, useNames=NA)
{
    stopifnot(is(x, "NaMatrix"))
    check_svt_version(x)
    if (!isTRUEorFALSE(na.rm))
        stop(wmsg("'na.rm' must be TRUE or FALSE"))
    useNames <- normarg_useNames(useNames)
    x_nrow <- nrow(x)
    x_ncol <- ncol(x)
    if (x_nrow == 0L) {
        ans <- rep.int(NA_real_, x_ncol)
    } else {
        ans <- numeric(x_ncol)
        if (!is.null(x@NaSVT)) {
            ans <- vapply(seq_along(x@NaSVT),
                function(i) {
                    lv <- x@NaSVT[[i]]
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

#setMethod("colMedians", "NaArray",
#    function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
#    {
#        check_unused_arguments(...)
#        stopifnot_2D_object(x, "colMedians", "NaArray", "NaMatrix")
#        check_rows_cols(rows, cols, "colMedians", "NaArray")
#        .colMedians_NaMatrix(x, na.rm=na.rm, useNames=useNames)
#    }
#)

#setMethod("rowMedians", "NaArray",
#    function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
#    {
#        check_unused_arguments(...)
#        stopifnot_2D_object(x, "rowMedians", "NaArray", "NaMatrix")
#        check_rows_cols(rows, cols, "rowMedians", "NaArray")
#        tx <- t(x)
#        .colMedians_NaMatrix(tx, na.rm=na.rm, useNames=useNames, ...)
#    }
#)

