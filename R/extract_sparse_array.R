### =========================================================================
### extract_sparse_array()
### -------------------------------------------------------------------------
###


### Similar to S4Arrays:::.contact_author_msg1().
.contact_author_msg2 <- function(.Generic, x_class)
{
    msg <- c("Please contact the authors/maintainers of the ",
             x_class, " class")
    class_package <- attr(x_class, "package")
    if (!is.null(class_package))
        msg <- c(msg, " (defined in the ", class_package, " package)")
    c(msg, " about this, and point them to the man page for the ",
           .Generic, "() generic function defined in the SparseArray ",
           "package ('?SparseArray::", .Generic, "').")
}

check_returned_SparseArray <- function(ans, expected_dim, .Generic, x_class)
{
    if (!is(ans, "SparseArray"))
        stop(wmsg("The ", .Generic, "() method for ", x_class, " ",
                  "objects didn't return a SparseArray object. ",
                  .Generic, "() methods should **always** return a ",
                  "SparseArray object. ",
                  .contact_author_msg2(.Generic, x_class)))
    if (!identical(dim(ans), expected_dim))
        stop(wmsg("The ", .Generic, "() method for ", x_class, " objects ",
                  "returned a SparseArray object with incorrect ",
                  "dimensions. ", .contact_author_msg2(.Generic, x_class)))
    ans
}

### extract_sparse_array() is the workhorse behind read_block_as_sparse(),
### and, more generally, behind efficient subsetting of sparse array objects.
### Similar to S4Arrays::extract_array(), except that:
###   (1) The extracted array data must be returned as a SparseArray object.
###       Methods should always operate on the sparse representation of the
###       data and never "expand" it, that is, never turn it into a dense
###       representation, e.g. with as.array(), as this would defeat the
###       purpose of read_block_as_sparse().
###   (2) It should **always** be called on an array-like object 'x' for
###       which 'is_sparse(x)' is TRUE.
###   (3) The subscripts in 'index' should NOT contain duplicates.
### IMPORTANT NOTE: For the sake of efficiency, (2) and (3) are NOT checked
### and are the responsibility of the user. We'll refer to (2) and (3) as
### the "extract_sparse_array() contract".

setGeneric("extract_sparse_array", signature="x",
    function(x, index)
    {
        x_dim <- dim(x)
        if (is.null(x_dim))
            stop(wmsg("first argument to extract_sparse_array() ",
                      "must be an array-like object"))
        ans <- standardGeneric("extract_sparse_array")
        expected_dim <- S4Arrays:::get_Nindex_lengths(index, x_dim)
        check_returned_SparseArray(ans, expected_dim,
                                   "extract_sparse_array", class(x))
    }
)

### S4Arrays:::subset_by_Nindex() uses `[` internally to perform the
### subsetting, so this default extract_sparse_array() method will work
### on any object 'x' that supports `[` and coercion to SparseArray.
setMethod("extract_sparse_array", "ANY",
    function(x, index)
    {
        slice <- S4Arrays:::subset_by_Nindex(x, index)
        as(slice, "SparseArray")
    }
)

