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

.crossprod2_SparseMatrix_matrix <- function(x, y, transpose.y=FALSE)
{
    if (is(x, "SVT_SparseMatrix")) {
        check_svt_version(x)
    } else {
        x <- as(x, "SVT_SparseMatrix")
    }
    stopifnot(is.matrix(y), isTRUEorFALSE(transpose.y))
    if (transpose.y) {
        if (nrow(x) != ncol(y))
            stop(wmsg("non-conformable arguments"))
        ans_dim <- c(ncol(x), nrow(y))
        ans_dimnames <- list(colnames(x), rownames(y))
    } else {
        if (nrow(x) != nrow(y))
            stop(wmsg("non-conformable arguments"))
        ans_dim <- c(ncol(x), ncol(y))
        ans_dimnames <- list(colnames(x), colnames(y))
    }
    if (type(x) == type(y)) {
        .check_crossprod_input_type(type(x))
    } else {
        xy_type <- type(c(vector(type(x)), vector(type(y))))
        .check_crossprod_input_type(xy_type)
        type(x) <- type(y) <- xy_type
    }
    ans_type <- "double"
    ans_dimnames <- S4Arrays:::simplify_NULL_dimnames(ans_dimnames)
    SparseArray.Call("C_crossprod2_SVT_mat",
                     x@dim, x@type, x@SVT, y, transpose.y,
                     ans_type, ans_dimnames)
}

.crossprod2_matrix_SparseMatrix <- function(x, y, transpose.x=FALSE)
{
    if (is(y, "SVT_SparseMatrix")) {
        check_svt_version(y)
    } else {
        y <- as(y, "SVT_SparseMatrix")
    }
    stopifnot(is.matrix(x), isTRUEorFALSE(transpose.x))
    if (transpose.x) {
        if (ncol(x) != nrow(y))
            stop(wmsg("non-conformable arguments"))
        ans_dim <- c(nrow(x), ncol(y))
        ans_dimnames <- list(rownames(x), colnames(y))
    } else {
        if (nrow(x) != nrow(y))
            stop(wmsg("non-conformable arguments"))
        ans_dim <- c(ncol(x), ncol(y))
        ans_dimnames <- list(colnames(x), colnames(y))
    }
    if (type(x) == type(y)) {
        .check_crossprod_input_type(type(x))
    } else {
        xy_type <- type(c(vector(type(x)), vector(type(y))))
        .check_crossprod_input_type(xy_type)
        type(x) <- type(y) <- xy_type
    }
    ans_type <- "double"
    ans_dimnames <- S4Arrays:::simplify_NULL_dimnames(ans_dimnames)
    SparseArray.Call("C_crossprod2_mat_SVT",
                     x, y@dim, y@type, y@SVT, transpose.x,
                     ans_type, ans_dimnames)
}

.crossprod2_SparseMatrix_SparseMatrix <- function(x, y=NULL)
{
    if (is(x, "SVT_SparseMatrix")) {
        check_svt_version(x)
    } else {
        x <- as(x, "SVT_SparseMatrix")
    }
    if (is(y, "SVT_SparseMatrix")) {
        check_svt_version(y)
    } else {
        y <- as(y, "SVT_SparseMatrix")
    }
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
    SparseArray.Call("C_crossprod2_SVT_SVT",
                     x@dim, x@type, x@SVT, y@dim, y@type, y@SVT,
                     ans_type, ans_dimnames)
}

.crossprod1_SparseMatrix <- function(x, y=NULL)
{
    if (is(x, "SVT_SparseMatrix")) {
        check_svt_version(x)
    } else {
        x <- as(x, "SVT_SparseMatrix")
    }
    stopifnot(is.null(y))
    ans_dim <- c(ncol(x), ncol(x))
    .check_crossprod_input_type(type(x))
    ans_type <- "double"
    ans_dimnames <- list(colnames(x), colnames(x))
    ans_dimnames <- S4Arrays:::simplify_NULL_dimnames(ans_dimnames)
    SparseArray.Call("C_crossprod1_SVT",
                     x@dim, x@type, x@SVT,
                     ans_type, ans_dimnames)
}

setMethod("crossprod", c("SparseMatrix", "matrix"),
    .crossprod2_SparseMatrix_matrix
)

setMethod("crossprod", c("matrix", "SparseMatrix"),
    .crossprod2_matrix_SparseMatrix
)

setMethod("crossprod", c("SparseMatrix", "SparseMatrix"),
    .crossprod2_SparseMatrix_SparseMatrix
)

setMethod("crossprod", c("SparseMatrix", "ANY"),
    function(x, y=NULL) .crossprod2_SparseMatrix_SparseMatrix(x, y)
)

setMethod("crossprod", c("ANY", "SparseMatrix"),
    function(x, y=NULL) .crossprod2_SparseMatrix_SparseMatrix(x, y)
)

setMethod("crossprod", c("SparseMatrix", "missing"),
    .crossprod1_SparseMatrix
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### tcrossprod()
###

setMethod("tcrossprod", c("SparseMatrix", "matrix"),
    function(x, y=NULL) .crossprod2_SparseMatrix_matrix(t(x), y,
                                                        transpose.y=TRUE)
)

setMethod("tcrossprod", c("matrix", "SparseMatrix"),
    function(x, y=NULL) .crossprod2_matrix_SparseMatrix(x, t(y),
                                                        transpose.x=TRUE)
)

setMethod("tcrossprod", c("SparseMatrix", "SparseMatrix"),
    function(x, y=NULL) .crossprod2_SparseMatrix_SparseMatrix(t(x), t(y))
)

setMethod("tcrossprod", c("SparseMatrix", "ANY"),
    function(x, y=NULL) .crossprod2_SparseMatrix_SparseMatrix(t(x), t(y))
)

setMethod("tcrossprod", c("ANY", "SparseMatrix"),
    function(x, y=NULL) .crossprod2_SparseMatrix_SparseMatrix(t(x), t(y))
)

setMethod("tcrossprod", c("SparseMatrix", "missing"),
    function(x, y=NULL) .crossprod1_SparseMatrix(t(x))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Matrix multiplication
###

setMethod("%*%", c("SparseMatrix", "matrix"),
    function(x, y) .crossprod2_SparseMatrix_matrix(t(x), y)
)

setMethod("%*%", c("matrix", "SparseMatrix"),
    function(x, y) .crossprod2_matrix_SparseMatrix(x, y, transpose.x=TRUE)
)

setMethod("%*%", c("SparseMatrix", "SparseMatrix"),
    function(x, y) .crossprod2_SparseMatrix_SparseMatrix(t(x), y)
)

setMethod("%*%", c("SparseMatrix", "ANY"),
    function(x, y) .crossprod2_SparseMatrix_SparseMatrix(t(x), y)
)

setMethod("%*%", c("ANY", "SparseMatrix"),
    function(x, y) .crossprod2_SparseMatrix_SparseMatrix(t(x), y)
)

