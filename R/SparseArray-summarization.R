### =========================================================================
### Summarization methods for SparseArray objects
### -------------------------------------------------------------------------
###
### Summarization methods:
###   - anyNA()
###   - 'Summary' group: any(), all(), min(), max(), range(), sum(), prod()
###   - mean()
###   - Unary var(), sd()
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Workhorse behind all the summarization methods for SVT_SparseArray
### and NaArray objects
###

### 'center' ignored by all ops except "centered_X2_sum".
### Returns an integer or numeric vector of length 1 or 2.
summarize_SVT <- function(op, x, na.rm=FALSE, center=NULL)
{
    stopifnot(isSingleString(op), is(x, "SVT_SparseArray") || is(x, "NaArray"))
    check_svt_version(x)

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

    if (is(x, "NaArray")) {
        SparseArray.Call("C_summarize_SVT",
                         x@dim, x@type, x@NaSVT, TRUE, op, na.rm, center)
    } else {
        SparseArray.Call("C_summarize_SVT",
                         x@dim, x@type, x@SVT, FALSE, op, na.rm, center)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### anyNA(), countNAs()
###

.anyNA_SparseArray <- function(x, recursive=FALSE)
{
    if (!identical(recursive, FALSE))
        stop(wmsg("the anyNA() method for SparseArray objects ",
                  "does not support the 'recursive' argument"))

    if (is(x, "COO_SparseArray"))
        return(anyNA(x@nzdata))

    if (is(x, "SVT_SparseArray"))
        return(summarize_SVT("anyNA", x))

    stop(wmsg(class(x)[[1L]], " objects are not supported"))
}
setMethod("anyNA", "SparseArray", .anyNA_SparseArray)

### NOT USED! There's no countNAs() generic yet!
### TODO: Define the countNAs() in BiocGenerics, and the colCountNAs() and
### rowCountNAs() generics in MatrixGenerics.
.countNAs_SparseArray <- function(x, recursive=FALSE)
{
    if (!identical(recursive, FALSE))
        stop(wmsg("the countNAs() method for SparseArray objects ",
                  "does not support the 'recursive' argument"))

    if (is(x, "COO_SparseArray"))
        return(sum(is.na(x@nzdata)))  # or do 'countNAs(x@nzdata)' when it
                                      # becomes available

    if (is(x, "SVT_SparseArray"))
        return(summarize_SVT("countNAs", x))

    stop(wmsg(class(x)[[1L]], " objects are not supported"))
}
#setMethod("countNAs", "SparseArray", .countNAs_SparseArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Summary' group
###

.summarize_COO <- function(op, x, na.rm=FALSE)
{
    stopifnot(isSingleString(op), is(x, "COO_SparseArray"))
    GENERIC <- match.fun(op)
    ## Whether 'x' contains zeros or not doesn't make a difference for
    ## sum() and any().
    if (op %in% c("sum", "any"))
        return(GENERIC(x@nzdata, na.rm=na.rm))
    ## Of course a typical COO_SparseArray object "contains" zeros
    ## (i.e. it would contain zeros if we converted it to a dense
    ## representation with as.array()). However, this is not guaranteed
    ## so we need to make sure to properly handle the case where it
    ## doesn't (admittedly unusual and definitely an inefficient way
    ## to represent dense data!)
    x_has_zeros <- length(x@nzdata) < length(x)
    if (!x_has_zeros)
        return(GENERIC(x@nzdata, na.rm=na.rm))
    x_type <- typeof(x@nzdata)
    if (op == "all") {
        ## Mimic what 'all(as.array(x))' would do.
        if (x_type == "double")
            warning("coercing argument of type 'double' to logical")
        return(FALSE)
    }
    zero <- vector(x_type, length=1L)
    GENERIC(zero, x@nzdata, na.rm=na.rm)
}

setMethod("Summary", "COO_SparseArray",
    function(x, ..., na.rm=FALSE)
    {
        if (length(list(...)) != 0L)
            stop(wmsg("the ", .Generic, "() method for COO_SparseArray ",
                      "objects only accepts a single object"))
        .summarize_COO(.Generic, x, na.rm=na.rm)
    }
)

setMethod("Summary", "SVT_SparseArray",
    function(x, ..., na.rm=FALSE)
    {
        if (length(list(...)) != 0L)
            stop(wmsg("the ", .Generic, "() method for SVT_SparseArray ",
                      "objects only accepts a single object"))
        summarize_SVT(.Generic, x, na.rm=na.rm)
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
        stop(wmsg("the range() method for COO_SparseArray objects ",
                  "only accepts a single object"))
    x <- objects[[1L]]
    x_has_zeros <- length(x@nzdata) < length(x)
    if (!x_has_zeros)
        return(range(x@nzdata, na.rm=na.rm, finite=finite))
    zero <- vector(typeof(x@nzdata), length=1L)
    range(zero, x@nzdata, na.rm=na.rm, finite=finite)
}
### The signature of all the members in the 'Summary' group generic is
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
        stop(wmsg("the range() method for SVT_SparseArray objects ",
                  "does not support the 'finite' argument"))
    objects <- list(...)
    if (length(objects) != 1L)
        stop(wmsg("the range() method for SVT_SparseArray objects ",
                  "only accepts a single object"))
    x <- objects[[1L]]
    summarize_SVT("range", x, na.rm=na.rm)
}
### The signature of all the members in the 'Summary' group generic is
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

.mean_SparseArray <- function(x, na.rm=FALSE)
{
    summarize_SVT("mean", x, na.rm=na.rm)
}

### S3/S4 combo for mean.SparseArray
mean.SparseArray <- function(x, na.rm=FALSE, ...)
    .mean_SparseArray(x, na.rm=na.rm, ...)
setMethod("mean", "SparseArray", .mean_SparseArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### var(), sd()
###

setMethod("var", c("SparseArray", "ANY"),
    function(x, y=NULL, na.rm=FALSE, use)
    {
        if (!is.null(y))
            stop(wmsg("the var() method for SparseArray objects ",
                      "does not support the 'y' argument"))
        if (!missing(use))
            stop(wmsg("the var() method for SparseArray objects ",
                      "does not support the 'use' argument"))
        summarize_SVT("var1", x, na.rm=na.rm)
    }
)

setMethod("sd", "SparseArray",
    function(x, na.rm=FALSE) summarize_SVT("sd1", x, na.rm=na.rm)
)

