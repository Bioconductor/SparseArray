
### A silly trick used only to trigger an error when the function is called
### with arguments passed to it.
check_unused_arguments <- function() NULL

coercion_can_introduce_zeros <- function(from_type, to_type)
{
    if (!isSingleString(from_type))
        stop(wmsg("'from_type' must be a single string"))
    if (!isSingleString(to_type))
        stop(wmsg("'to_type' must be a single string"))
    if (!(to_type %in% c("double", "logical")))
        stop(wmsg("'to_type' must be \"double\" or \"logical\""))
    .Call2("C_coercion_can_introduce_zeros", from_type, to_type,
                                             PACKAGE="SparseArray")
}

coercion_can_introduce_NAs <- function(from_type, to_type)
{
    if (!isSingleString(from_type))
        stop(wmsg("'from_type' must be a single string"))
    if (!isSingleString(to_type))
        stop(wmsg("'to_type' must be a single string"))
    .Call2("C_coercion_can_introduce_NAs", from_type, to_type,
                                           PACKAGE="SparseArray")
}

Nindex2Noffs <- function(Nindex)
{
    stopifnot(is.list(Nindex))
    lapply(Nindex,
        function(subscript)
            if (is.null(subscript)) NULL else subscript - 1L
    )
}

vector_of_zeros <- function(mode="logical", length=0L)
{
    vector(mode=mode, length=length)
}

vector_of_ones <- function(mode="logical", length=0L)
{
    as.fun <- base::get(paste0("as.", mode), envir=asNamespace("base"),
                        mode="function")
    rep.int(as.fun(1L), length)
}

### Can be used on the leaves of a SparseArray or NaArray object.
get_leaf_nzvals <- function(leaf, type)
{
    stopifnot(is.list(leaf), length(leaf) == 2L)
    nzvals <- leaf[[1L]]
    if (!is.null(nzvals))
        return(nzvals)
    vector_of_ones(type, length(leaf[[2L]]))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### padded_median()
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
        return(sort.int(x, partial=partial)[partial])
    }
    i1 <- n %/% 2L - padding
    i2 <- i1 + 1L
    mean(sort.int(x, partial=i1:i2)[i1:i2])
}

### Equivalent to 'median(c(x, integer(padding)), ...)' but doesn't actually
### realize the padding with zeros.
padded_median <- function(x, padding=0L, na.rm=FALSE)
{
    stopifnot(is.numeric(x), isSingleInteger(padding), isTRUEorFALSE(na.rm))
    if (na.rm) {
        x <- x[!is.na(x)]
    } else {
        if (anyNA(x))
            return(NA_real_)
    }
    x_len <- length(x)
    n <- x_len + padding
    if (n == 0L)
        return(NA_real_)
    if (padding > x_len)
        return(0)

    ## Handle case where we have more positive values than non-positive values.
    is_pos <- x > 0L
    pos_count <- sum(is_pos)
    nonpos_count <- n - pos_count
    if (pos_count > nonpos_count) {
        ans <- .positive_padded_median(x[is_pos], padding=nonpos_count)
        return(ans)
    }

    ## Handle case where we have more negative values than non-negative values.
    neg_count <- x_len - pos_count
    nonneg_count <- n - neg_count
    if (neg_count > nonneg_count) {
        ans <- - .positive_padded_median(-x[!is_pos], padding=nonneg_count)
        return(ans)
    }

    if (n %% 2L == 1L)
        return(0)

    half <- n %/% 2L
    if (pos_count == half) {
        right <- min(x[is_pos])
    } else {
        right <- 0
    }
    if (neg_count == half) {
        left <- max(x[!is_pos])
    } else {
        left <- 0
    }
    (left + right) * 0.5
}

