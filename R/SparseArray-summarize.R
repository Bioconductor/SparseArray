### =========================================================================
### Summarization methods for SparseArray objects
### -------------------------------------------------------------------------
###
### Summarization methods:
###   - Summary group generic: min(), max(), range(), sum(), prod(),
###                            any(), all()
###   - mean(), anyNA()
###   - Unary var(), sd()
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Summary group generic
###

.summarize_COO_SparseArray <- function(op, x, na.rm=FALSE)
{
    stopifnot(is(x, "COO_SparseArray"))
    GENERIC <- match.fun(op)
    ## Whether 'x' contains zeros or not doesn't make a difference for
    ## sum() and any().
    if (op %in% c("sum", "any"))
        return(GENERIC(x@nzvals, na.rm=na.rm))
    ## Of course a typical COO_SparseArray object "contains" zeros
    ## (i.e. it would contain zeros if we converted it to a dense
    ## representation with sparse2dense()). However, this is not
    ## guaranteed so we need to make sure to properly handle the case
    ## where it doesn't (admittedly unusual and definitely an inefficient
    ## way to represent dense data!)
    x_has_zeros <- length(x@nzvals) < length(x)
    if (!x_has_zeros)
        return(GENERIC(x@nzvals, na.rm=na.rm))
    x_type <- typeof(x@nzvals)
    if (op == "all") {
        ## Mimic what 'all(sparse2dense(x))' would do.
        if (x_type == "double")
            warning("coercing argument of type 'double' to logical")
        return(FALSE)
    }
    zero <- vector(x_type, length=1L)
    GENERIC(zero, x@nzvals, na.rm=na.rm)
}

setMethod("Summary", "COO_SparseArray",
    function(x, ..., na.rm=FALSE)
    {
        if (length(list(...)) != 0L)
            stop(wmsg(.Generic, "() method for COO_SparseArray objects ",
                      "only accepts a single object"))
        .summarize_COO_SparseArray(.Generic, x, na.rm=na.rm)
    }
)

### 'shift' ignored by all ops except "sum_shifted_X2".
### Returns an integer or numeric vector of length 1 or 2.
### If 'na_rm' is TRUE, then the "na_rm_count" attribute is set on the
### returned vector.
.summarize_SVT_SparseArray <- function(op, x, na.rm=FALSE, shift=0.0)
{
    stopifnot(is(x, "SVT_SparseArray"))
    .Call2("C_summarize_SVT_SparseArray",
           x@dim, x@type, x@SVT, op, na.rm, shift, PACKAGE="SparseArray")
}

setMethod("Summary", "SVT_SparseArray",
    function(x, ..., na.rm=FALSE)
    {
        if (length(list(...)) != 0L)
            stop(wmsg(.Generic, "() method for SVT_SparseArray objects ",
                      "only accepts a single object"))
        ans <- .summarize_SVT_SparseArray(.Generic, x, na.rm=na.rm)
        if (na.rm)
            attr(ans, "na_rm_count") <- NULL
        ans
    }
)

### We override the range() methods defined via the Summary() methods
### above because we want to support the 'finite' argument like S3 method
### base::range.default() does. One might wonder why base::range.default()
### supports the 'finite' argument but min() and max() don't. Or more
### precisely, they seem to take it but they don't do exactly the same thing
### with it:
###
###     > max(c(0, -Inf), finite=TRUE)
###     [1] 1
###
### Another story for another day...

### S3/S4 combo for range.COO_SparseArray
range.COO_SparseArray <- function(..., na.rm=FALSE, finite=FALSE)
{
    objects <- list(...)
    if (length(objects) != 1L)
        stop(wmsg("range() method for COO_SparseArray objects ",
                  "only accepts a single object"))
    x <- objects[[1L]]
    x_has_zeros <- length(x@nzvals) < length(x)
    if (!x_has_zeros)
        return(range(x@nzvals, na.rm=na.rm, finite=finite))
    zero <- vector(typeof(x@nzvals), length=1L)
    range(zero, x@nzvals, na.rm=na.rm, finite=finite)
}
### The signature of all the members in the Summary group generic is
### 'x, ..., na.rm' (see getGeneric("range")) which means that methods
### cannot add arguments after 'na.rm'. So we add the 'finite' argument
### before.
setMethod("range", "COO_SparseArray",
    function(x, ..., finite=FALSE, na.rm=FALSE)
        range.COO_SparseArray(x, ..., na.rm=na.rm, finite=finite)

)

### S3/S4 combo for range.SVT_SparseArray
range.SVT_SparseArray <- function(..., na.rm=FALSE, finite=FALSE)
{
    if (!identical(finite, FALSE))
        stop(wmsg("range() method for SVT_SparseArray objects ",
                  "does not support the 'finite' argument"))
    objects <- list(...)
    if (length(objects) != 1L)
        stop(wmsg("range() method for SVT_SparseArray objects ",
                  "only accepts a single object"))
    x <- objects[[1L]]
    ans <- .summarize_SVT_SparseArray("range", x, na.rm=na.rm)
    if (na.rm)
        attr(ans, "na_rm_count") <- NULL
    ans
}

### The signature of all the members in the Summary group generic is
### 'x, ..., na.rm' (see getGeneric("range")) which means that methods
### cannot add arguments after 'na.rm'. So we add the 'finite' argument
### before.
setMethod("range", "SVT_SparseArray",
    function(x, ..., finite=FALSE, na.rm=FALSE)
        range.SVT_SparseArray(x, ..., na.rm=na.rm, finite=finite)

)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mean()
###

### TODO: Maybe introduce a new generic for this e.g. countNAs()?
.count_SparseArray_NAs <- function(x)
{
    if (is(x, "COO_SparseArray"))
        return(sum(is.na(x@nzvals)))

    if (is(x, "SVT_SparseArray"))
        return(.Call2("C_count_SVT_SparseArray_NAs",
                      x@dim, x@type, x@SVT, PACKAGE="SparseArray"))

    stop(wmsg(class(x)[[1L]], " objects are not supported"))
}

.mean_SparseArray <- function(x, na.rm=FALSE)
{
    nval <- length(x)
    if (is(x, "COO_SparseArray")) {
        sum_X <- sum(x, na.rm=na.rm)
        if (na.rm)
            nval <- nval - .count_SparseArray_NAs(x)
    } else if (is(x, "SVT_SparseArray")) {
        sum_X <- .summarize_SVT_SparseArray("sum", x, na.rm=na.rm)
        if (na.rm)
            nval <- nval - attr(sum_X, "na_rm_count")
    } else {
        stop(wmsg(class(x)[[1L]], " objects are not supported"))
    }
    as.double(sum_X) / nval
}

### S3/S4 combo for mean.SparseArray
mean.SparseArray <- function(x, na.rm=FALSE, ...)
    .mean_SparseArray(x, na.rm=na.rm, ...)
setMethod("mean", "SparseArray", .mean_SparseArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### anyNA()
###

setMethod("anyNA", "COO_SparseArray",
    function(x, recursive=FALSE) anyNA(x@nzvals, recursive=recursive)
)

setMethod("anyNA", "SVT_SparseArray",
    function(x, recursive=FALSE)
    {
        if (!identical(recursive, FALSE))
            stop(wmsg("anyNA() method for SVT_SparseArray objects ",
                      "does not support the 'recursive' argument"))
        .Call2("C_anyNA_SVT_SparseArray",
               x@dim, x@type, x@SVT, PACKAGE="SparseArray")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### var()
###

.var_SparseArray <- function(x, na.rm=FALSE, method=0L)
{
    stopifnot(is(x, "SparseArray"))
    if (method == 0L) {
	## A single pass on 'x'.
        ans <- .summarize_SVT_SparseArray("var", x, na.rm=na.rm)
        if (na.rm)
            attr(ans, "na_rm_count") <- NULL
        return(ans)
    }
    if (method == 1L) {
        ## Uses primary variance formula:
        ##     sum((x - mean(x))^2) / nval
        ## Two passes on 'x'.
        nval <- length(x)
        sum_X <- .summarize_SVT_SparseArray("sum", x, na.rm=na.rm)
        if (na.rm)
            nval <- nval - attr(sum_X, "na_rm_count")
        mu <- sum_X / nval
        sum_shifted_X2 <- .summarize_SVT_SparseArray("sum_shifted_X2",
                                                     x, na.rm=na.rm,
                                                     shift=mu)
        sum_shifted_X2 <- sum_shifted_X2 + mu * mu * (length(x) - nzcount(x))
        attributes(sum_shifted_X2) <- NULL
        return(sum_shifted_X2 / (nval - 1))
    }
    if (method == 2L) {
        ## Uses secondary variance formula:
        ##     (sum(x^2) − (sum(x)^2) / nval) / (nval − 1)
	## A single pass on 'x' so faster than method 1 but numerically
	## instable!
        nval <- length(x)
        sum_X_X2 <- .summarize_SVT_SparseArray("sum_X_X2", x, na.rm=na.rm)
        if (na.rm)
            nval <- nval - attr(sum_X_X2, "na_rm_count")
        if (nval <= 1L)
            return(NA_real_)
        sum_X <- as.double(sum_X_X2[[1L]])
        sum_X2 <- as.double(sum_X_X2[[2L]])
        return((sum_X2 - sum_X * sum_X / nval) / (nval - 1))
    }
    ## Will work out-of-the-box on any object 'x' that supports
    ## .count_SparseArray_NAs(), mean(), `-`, `^`, and sum().
    ## Won't be very efficient though because it performs 5 passes: 3 passes
    ## on 'x' and 2 passes on a modified version of 'x'!
    nval <- length(x)
    if (na.rm)
        nval <- nval - .count_SparseArray_NAs(x)  # 1st pass on 'x'
    if (nval <= 1L)
        return(NA_real_)
    s <- sum((x - mean(x, na.rm=na.rm))^2L, na.rm=na.rm)  # 4 passes on 'x'
    s / (nval - 1L)
}

setMethod("var", "SparseArray",
    function(x, y = NULL, na.rm = FALSE, use)
    {
        if (!is.null(y))
            stop(wmsg("var() method for SparseArray objects ",
                      "does not support the 'y' argument"))
        if (!missing(use))
            stop(wmsg("var() method for SparseArray objects ",
                      "does not support the 'use' argument"))
        .var_SparseArray(x, na.rm=na.rm, method=2L)
    }
)

