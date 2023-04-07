### =========================================================================
### SparseMatrix crossprod(), tcrossprod(), and %*%
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### crossprod()
###

### Only input types "double" and "integer" are supported at the moment.
### TODO: Add "complex" and "logical" later.
.check_crossprod_input_type <- function(type)
{
    #supported_types <- c("double", "integer", "complex", "logical")
    supported_types <- c("double", "integer")
    if (!(type %in% supported_types))
        stop(wmsg("input objects must be of type() \"double\" or \"integer\""))
                  #"\"integer\", \"complex\", or \"logical\""))
}

.crossprod2_SVT_mat <- function(x, y=NULL)
{
    stopifnot(is(x, "SVT_SparseMatrix"), is.matrix(y))
    if (nrow(x) != nrow(y))
        stop(wmsg("non-conformable arguments"))
    ans_dim <- c(ncol(x), ncol(y))
    if (type(x) == type(y)) {
        .check_crossprod_input_type(type(x))
    } else {
        xy_type <- type(c(vector(type(x)), vector(type(y))))
        .check_crossprod_input_type(xy_type)
        type(x) <- type(y) <- xy_type
    }
    ans_type <- "double"
    ans_dimnames <- list(colnames(x), colnames(y))
    ans_dimnames <- S4Arrays:::simplify_NULL_dimnames(ans_dimnames)
    .Call2("C_crossprod2_SVT_mat", x@dim, x@type, x@SVT, y,
           ans_type, ans_dimnames,
           PACKAGE="SparseArray")
}

.crossprod2_mat_SVT <- function(x, y=NULL)
{
    stopifnot(is.matrix(x), is(y, "SVT_SparseMatrix"))
    if (nrow(x) != nrow(y))
        stop(wmsg("non-conformable arguments"))
    ans_dim <- c(ncol(x), ncol(y))
    if (type(x) == type(y)) {
        .check_crossprod_input_type(type(x))
    } else {
        xy_type <- type(c(vector(type(x)), vector(type(y))))
        .check_crossprod_input_type(xy_type)
        type(x) <- type(y) <- xy_type
    }
    ans_type <- "double"
    ans_dimnames <- list(colnames(x), colnames(y))
    ans_dimnames <- S4Arrays:::simplify_NULL_dimnames(ans_dimnames)
    .Call2("C_crossprod2_mat_SVT", x, y@dim, y@type, y@SVT,
           ans_type, ans_dimnames,
           PACKAGE="SparseArray")
}

.crossprod2_SVT_SVT <- function(x, y=NULL)
{
    stopifnot(is(x, "SVT_SparseMatrix"), is(y, "SVT_SparseMatrix"))
    if (nrow(x) != nrow(y))
        stop(wmsg("non-conformable arguments"))
    ans_dim <- c(ncol(x), ncol(y))
    if (type(x) == type(y)) {
        .check_crossprod_input_type(type(x))
    } else {
        xy_type <- type(c(vector(type(x)), vector(type(y))))
        .check_crossprod_input_type(xy_type)
        type(x) <- type(y) <- xy_type
    }
    ans_type <- "double"
    ans_dimnames <- list(colnames(x), colnames(y))
    ans_dimnames <- S4Arrays:::simplify_NULL_dimnames(ans_dimnames)
    .Call2("C_crossprod2_SVT_SVT", x@dim, x@type, x@SVT, y@dim, y@type, y@SVT,
           ans_type, ans_dimnames,
           PACKAGE="SparseArray")
}

.crossprod1_SVT <- function(x, y=NULL)
{
    stopifnot(is(x, "SVT_SparseMatrix"), is.null(y))
    ans_dim <- c(ncol(x), ncol(x))
    .check_crossprod_input_type(type(x))
    ans_type <- "double"
    ans_dimnames <- list(colnames(x), colnames(x))
    ans_dimnames <- S4Arrays:::simplify_NULL_dimnames(ans_dimnames)
    .Call2("C_crossprod1_SVT", x@dim, x@type, x@SVT, ans_type, ans_dimnames,
           PACKAGE="SparseArray")
}

setMethod("crossprod", c("SVT_SparseMatrix", "matrix"),
    .crossprod2_SVT_mat
)

setMethod("crossprod", c("matrix", "SVT_SparseMatrix"),
    .crossprod2_mat_SVT
)

setMethod("crossprod", c("SVT_SparseMatrix", "SVT_SparseMatrix"),
    .crossprod2_SVT_SVT
)

setMethod("crossprod", c("SVT_SparseMatrix", "ANY"),
    function(x, y=NULL) .crossprod2_SVT_SVT(x, as(y, "SVT_SparseMatrix"))
)

setMethod("crossprod", c("ANY", "SVT_SparseMatrix"),
    function(x, y=NULL) .crossprod2_SVT_SVT(as(x, "SVT_SparseMatrix"), y)
)

setMethod("crossprod", c("SVT_SparseMatrix", "missing"),
    .crossprod1_SVT
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### tcrossprod()
###

setMethod("tcrossprod", c("SVT_SparseMatrix", "matrix"),
    function(x, y=NULL) .crossprod2_SVT_mat(t(x), t(y))
)

setMethod("tcrossprod", c("matrix", "SVT_SparseMatrix"),
    function(x, y=NULL) .crossprod2_mat_SVT(t(x), t(y))
)

setMethod("tcrossprod", c("SVT_SparseMatrix", "SVT_SparseMatrix"),
    function(x, y=NULL) .crossprod2_SVT_SVT(t(x), t(y))
)

setMethod("tcrossprod", c("SVT_SparseMatrix", "ANY"),
    function(x, y=NULL) .crossprod2_SVT_SVT(t(x), t(as(y, "SVT_SparseMatrix")))
)

setMethod("tcrossprod", c("ANY", "SVT_SparseMatrix"),
    function(x, y=NULL) .crossprod2_SVT_SVT(t(as(x, "SVT_SparseMatrix")), t(y))
)

setMethod("tcrossprod", c("SVT_SparseMatrix", "missing"),
    function(x, y=NULL) .crossprod1_SVT(t(x))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Matrix multiplication
###

setMethod("%*%", c("SVT_SparseMatrix", "matrix"),
    function(x, y) .crossprod2_SVT_mat(t(x), y)
)

setMethod("%*%", c("matrix", "SVT_SparseMatrix"),
    function(x, y) .crossprod2_mat_SVT(t(x), y)
)

setMethod("%*%", c("SVT_SparseMatrix", "SVT_SparseMatrix"),
    function(x, y) .crossprod2_SVT_SVT(t(x), y)
)

setMethod("%*%", c("SVT_SparseMatrix", "ANY"),
    function(x, y) .crossprod2_SVT_SVT(t(x), as(y, "SVT_SparseMatrix"))
)

setMethod("%*%", c("ANY", "SVT_SparseMatrix"),
    function(x, y) .crossprod2_SVT_SVT(t(as(x, "SVT_SparseMatrix")), y)
)

