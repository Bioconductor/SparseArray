### =========================================================================
### randomSparseArray() and poissonSparseArray()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### randomSparseArray()
###

### Returns an SVT_SparseArray object of type "double".
randomSparseArray <- function(dim, density=0.05)
{
    if (!is.numeric(dim))
        stop(wmsg("'dim' must be an integer vector"))
    if (!is.integer(dim))
        dim <- as.integer(dim)
    if (!isSingleNumber(density) || density < 0 || density > 1)
        stop(wmsg("'density' must be a number >= 0 and <= 1"))

    ## Start with an empty sparse array.
    ans <- new_SVT_SparseArray(dim, type="double")

    ## Add the nonzero values to it.
    ans_len <- length(ans)
    nzcount <- as.integer(ans_len * density)
    Lindex <- sample.int(ans_len, nzcount)
    nzvals <- signif(rnorm(nzcount), 2)
    if (nzcount <= .Machine$integer.max) {
        ## Using an Mindex seems to be slightly faster (4%-5%) than using an
        ## Lindex but we can only do this when the resulting Mindex matrix
        ## has < 2^31 rows.
        ans[Lindex2Mindex(Lindex, dim(ans))] <- nzvals
    } else {
        ans[Lindex] <- nzvals
    }

    ans
}

randomSparseMatrix <- function(nrow=1L, ncol=1L, density=0.05)
{
    if (!isSingleNumber(nrow) || !isSingleNumber(ncol))
        stop(wmsg("'nrow' and 'ncol' must be single integers"))
    randomSparseArray(c(nrow, ncol), density=density)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### poissonSparseArray()
###

### Like stats::rpois() but slightly faster and implementation is much
### simpler. Only for 0 <= 'lambda' <= 4.
simple_rpois <- function(n, lambda)
    .Call2("C_simple_rpois", n, lambda, PACKAGE="SparseArray")

### Returns an SVT_SparseArray object of type "integer".
### Density of the returned object is expected to be about '1 - exp(-lambda)'.
### Default for 'lambda' is set to -log(0.95) which should produce an object
### with an expected density of 0.05.
poissonSparseArray <- function(dim, lambda=-log(0.95), density=NA)
{
    dim <- S4Arrays:::normarg_dim(dim)

    if (!missing(lambda) && !identical(density, NA))
        stop(wmsg("only one of 'lambda' and 'density' can be specified"))
    if (!missing(lambda)) {
        if (!isSingleNumber(lambda) || lambda < 0)
            stop(wmsg("'lambda' must be a non-negative number"))
    } else if (!identical(density, NA)) {
        if (!isSingleNumber(density) || density < 0 || density >= 1)
            stop(wmsg("'density' must be a number >= 0 and < 1"))
        lambda <- -log(1 - density)
    }

    ans_SVT <- .Call2("C_poissonSparseArray", dim, lambda,
                      PACKAGE="SparseArray")
    new_SVT_SparseArray(dim, type="integer", SVT=ans_SVT, check=FALSE)
}

### Replacement for rpois() when 'n' is big and 'lambda' is small.
### For example:
###     .sparse_rpois(3e9, 0.005)  # takes about 1 min. and uses < 1G of RAM
###     rpois(3e9, 0.005)          # takes about 55 sec. and uses 12G of RAM
.sparse_rpois <- function(n, lambda, chunksize=5e6L)
{
    if (n == 0L)
        return(list(integer(0), integer(0)))
    nzidx_env <- new.env(parent=emptyenv())
    nzvals_env <- new.env(parent=emptyenv())
    offset <- 0  # double to avoid integer overflow when n >= 2^31
    k <- 1L
    while (offset < n) {
        nn <- n - offset
        if (nn > chunksize)
            nn <- chunksize
        vals <- rpois(nn, lambda)
        nzidx <- which(vals != 0L)
        key <- sprintf("%04d", k)
        assign(key, offset + nzidx, envir=nzidx_env)
        assign(key, vals[nzidx], envir=nzvals_env)
        offset <- offset + nn
        k <- k + 1L
    }
    nzidx <- as.list(nzidx_env, all.names=TRUE, sorted=TRUE)
    nzidx <- unlist(nzidx, recursive=FALSE, use.names=FALSE)
    nzvals <- as.list(nzvals_env, all.names=TRUE, sorted=TRUE)
    nzvals <- unlist(nzvals, recursive=FALSE, use.names=FALSE)
    list(nzidx, nzvals)
}

### NOT exported.
### Solution based on .sparse_rpois(). Equivalent to poissonSparseArray()
### but slower and uses more memory e.g.
###
###     poissonSparseArray2(c(1e5, 2e4), density=0.02)
###
### is about 3x slower uses about 2.5x more memory than
###
###     poissonSparseArray(c(1e5, 2e4), density=0.02)
###
poissonSparseArray2 <- function(dim, lambda=-log(0.95), density=NA)
{
    dim <- S4Arrays:::normarg_dim(dim)

    if (!missing(lambda) && !identical(density, NA))
        stop(wmsg("only one of 'lambda' and 'density' can be specified"))
    if (!missing(lambda)) {
        if (!isSingleNumber(lambda) || lambda < 0)
            stop(wmsg("'lambda' must be a non-negative number"))
    } else {
        if (!isSingleNumber(density) || density < 0 || density >= 1)
            stop(wmsg("'density' must be a number >= 0 and < 1"))
        lambda <- -log(1 - density)
    }

    ## Start with an empty sparse array.
    ans <- new_SVT_SparseArray(dim, type="integer")

    ## Add the nonzero values to it.
    ans_len <- length(ans)
    srp <- .sparse_rpois(ans_len, lambda)
    ans[srp[[1L]]] <- srp[[2L]]

    ans
}

poissonSparseMatrix <- function(nrow=1L, ncol=1L, lambda=-log(0.95), density=NA)
{
    if (!isSingleNumber(nrow) || !isSingleNumber(ncol))
        stop(wmsg("'nrow' and 'ncol' must be single integers"))
    if (missing(lambda)) {
        poissonSparseArray(c(nrow, ncol), density=density)
    } else {
        poissonSparseArray(c(nrow, ncol), lambda=lambda, density=density)
    }
}

