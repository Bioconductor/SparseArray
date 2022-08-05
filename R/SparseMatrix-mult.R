### =========================================================================
### SparseMatrix crossprod(), tcrossprod(), and %*%
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### crossprod()
###

.check_crossprod_input_type <- function(type)
{
    supported_types <- c("double", "integer", "complex", "logical")
    if (!(type %in% supported_types))
        stop(wmsg("operands must be of type() \"double\", ",
                  "\"integer\", \"complex\", or \"logical\""))
}

.unary_SVT_SparseMatrix_crossprod <- function(x, y=NULL)
{
    stopifnot(is(x, "SVT_SparseMatrix"), is.null(y))
    ans_dim <- c(ncol(x), ncol(x))
    ans_dimnames <- list(colnames(x), colnames(x))
    ans_type <- type(x)
    .check_crossprod_input_type(ans_type)
    ans_SVT <- .Call2("C_SVT_SparseMatrix_crossprod",
                      x@dim, x@SVT, x@dim, x@SVT, ans_type,
                      PACKAGE="SparseArray")
    new_SVT_SparseArray(ans_dim, ans_dimnames, ans_type, ans_SVT, check=FALSE)
}

.binary_SVT_SparseMatrix_crossprod <- function(x, y=NULL)
{
    stopifnot(is(x, "SVT_SparseMatrix"), is(y, "SVT_SparseMatrix"))
    if (nrow(y) != nrow(x))
        stop(wmsg("non-conformable arguments"))
    ans_dim <- c(ncol(x), ncol(y))
    ans_dimnames <- list(colnames(x), colnames(y))
    if (type(x) == type(y)) {
        ans_type <- type(x)
        .check_crossprod_input_type(ans_type)
    } else {
        ans_type <- type(c(vector(type(x)), vector(type(y))))
        .check_crossprod_input_type(ans_type)
        type(x) <- type(y) <- ans_type
    }
    ans_SVT <- .Call2("C_SVT_SparseMatrix_crossprod",
                      x@dim, x@SVT, y@dim, y@SVT, ans_type,
                      PACKAGE="SparseArray")
    new_SVT_SparseArray(ans_dim, ans_dimnames, ans_type, ans_SVT, check=FALSE)
}

setMethod("crossprod", c("SVT_SparseMatrix", "missing"),
    .unary_SVT_SparseMatrix_crossprod
)
setMethod("crossprod", c("SVT_SparseMatrix", "SVT_SparseMatrix"),
    .binary_SVT_SparseMatrix_crossprod
)

