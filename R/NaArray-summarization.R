### =========================================================================
### Summarization methods for NaArray objects
### -------------------------------------------------------------------------
###
### Summarization methods:
###   - anyNA()
###   - 'Summary' group: any(), all(), min(), max(), range(), sum(), prod()
###   - mean()
###   - Unary var(), sd()
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### anyNA(), countNAs()
###

.anyNA_NaArray <- function(x, recursive=FALSE)
{
    if (!identical(recursive, FALSE))
        stop(wmsg("the anyNA() method for NaArray objects ",
                  "does not support the 'recursive' argument"))

    summarize_SVT("anyNA", x)
}
setMethod("anyNA", "NaArray", .anyNA_NaArray)

### NOT USED! There's no countNAs() generic yet!
### TODO: Define the countNAs() in BiocGenerics, and the colCountNAs() and
### rowCountNAs() generics in MatrixGenerics.
.countNAs_NaArray <- function(x, recursive=FALSE)
{
    if (!identical(recursive, FALSE))
        stop(wmsg("the countNAs() method for NaArray objects ",
                  "does not support the 'recursive' argument"))

    summarize_SVT("countNAs", x)
}
#setMethod("countNAs", "NaArray", .countNAs_NaArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Summary' group
###

setMethod("Summary", "NaArray",
    function(x, ..., na.rm=FALSE)
    {
        if (length(list(...)) != 0L)
            stop(wmsg("the ", .Generic, "() method for NaArray ",
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

### S3/S4 combo for range.NaArray
range.NaArray <- function(..., na.rm=FALSE, finite=FALSE)
{
    if (!identical(finite, FALSE))
        stop(wmsg("the range() method for NaArray objects ",
                  "does not support the 'finite' argument"))
    objects <- list(...)
    if (length(objects) != 1L)
        stop(wmsg("the range() method for NaArray objects ",
                  "only accepts a single object"))
    x <- objects[[1L]]
    summarize_SVT("range", x, na.rm=na.rm)
}
### The signature of all the members in the 'Summary' group generic is
### 'x, ..., na.rm' (see getGeneric("range")) which means that methods
### cannot add arguments after 'na.rm'. So we add the 'finite' argument
### before.
setMethod("range", "NaArray",
    function(x, ..., finite=FALSE, na.rm=FALSE)
        range.NaArray(x, ..., na.rm=na.rm, finite=finite)

)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mean()
###

.mean_NaArray <- function(x, na.rm=FALSE)
{
    summarize_SVT("mean", x, na.rm=na.rm)
}

### S3/S4 combo for mean.NaArray
mean.NaArray <- function(x, na.rm=FALSE, ...)
    .mean_NaArray(x, na.rm=na.rm, ...)
setMethod("mean", "NaArray", .mean_NaArray)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### var(), sd()
###

setMethod("var", c("NaArray", "ANY"),
    function(x, y=NULL, na.rm=FALSE, use)
    {
        if (!is.null(y))
            stop(wmsg("the var() method for NaArray objects ",
                      "does not support the 'y' argument"))
        if (!missing(use))
            stop(wmsg("the var() method for NaArray objects ",
                      "does not support the 'use' argument"))
        summarize_SVT("var1", x, na.rm=na.rm)
    }
)

setMethod("sd", "NaArray",
    function(x, na.rm=FALSE) summarize_SVT("sd1", x, na.rm=na.rm)
)

