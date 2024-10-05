IS_INTEL_MAC <- Sys.info()[["sysname"]] == "Darwin" &&
                Sys.info()[["machine"]] == "x86_64"

IS_ARM64_MAC <- Sys.info()[["sysname"]] == "Darwin" &&
                Sys.info()[["machine"]] == "arm64"

make_3D_logical_array <- function(background=FALSE)
{
    a <- array(background, 6:4,
               dimnames=list(letters[1:6], NULL, LETTERS[1:4]))
    a[c(2:3, 6), 2, 1] <- TRUE
    a[ , 3, 1] <- c(TRUE, FALSE, NA, FALSE, FALSE, TRUE)
    a[2:3, 4:5, 1] <- NA
    a[c(1, 6), 1 , 2] <- FALSE
    a[ 2, , 2] <- FALSE
    a[ , 3, 2] <- TRUE
    a[3:5, 4:5, 2] <- NA
    a[c(1, 6), -5, 3] <- NA
    a[3:5, 1:3, 3] <- FALSE
    a[-3, 5, 3] <- TRUE
    a
}

make_3D_integer_array <- function(background=0L)
{
    a <- array(background, 6:4,
               dimnames=list(letters[1:6], NULL, LETTERS[1:4]))
    a[c(2:3, 6), 2, 1] <- 101:103
    a[c(1, 6), 1 , 2] <- 201:202
    a[c(3:4, 6), c(2, 5), 2] <- a[ , 3, 2] <- a[6, 4, 2] <- 1L
    a[4, 5, 2] <- NA
    a[1:5, 5, 3] <- -(301:305)
    a[6, 5, 4] <- NA
    a
}

make_3D_double_array <- function(background=0.0)
{
    a <- array(background, c(7, 10, 3),
               dimnames=list(NULL, letters[1:10], LETTERS[1:3]))
    set.seed(123)
    a[5*(1:26)] <- runif(26, min=-5, max=10)
    a[2, c(1:4, 7:9), 1] <- c(NA, NaN, Inf, 3e9, 256, -0.999, -1)
    a
}

make_3D_complex_array <- function(background=0L)
{
    a <- make_3D_integer_array(background)
    set.seed(123)
    a[6*(1:20)] <- a[6*(1:20)] + runif(20, min=-5, max=10) * 1i
    a[2, , 4] <- c(NA, NaN,  Inf,    3e9, 256)
    a[3, , 4] <- c(NA, NaN, -Inf, -0.999,  -1) * 1i
    a
}

make_lacunar_leaf <- function(mode, nzoffs)
{
    stopifnot(isSingleString(mode))
    stopifnot(is.integer(nzoffs))
    if (SparseArray:::lacunar_mode_is_on()) {
        nzvals <- NULL
    } else {
        as.fun <- base::get(paste0("as.", mode), envir=asNamespace("base"),
                            mode="function")
        nzvals <- rep.int(as.fun(1L), length(nzoffs))
    }
    list(nzvals, nzoffs)
}

check_array_like_object <- function(object, expected_class, a0, strict=TRUE)
{
    expect_true(class(object) == expected_class)
    expect_true(validObject(object))
    expect_identical(dim(object), dim(a0))
    expect_identical(dimnames(object), dimnames(a0))
    expect_identical(type(object), type(a0))
    EXPECT_FUN <- if (strict) expect_identical else expect_equal
    EXPECT_FUN(as.array(object), a0)
}

check_SVT_SparseArray_object <- function(svt, a0, strict=TRUE)
{
    Class0 <- S4Vectors:::capitalize(class(a0)[[1]])
    expected_class <- paste0("SVT_Sparse", Class0)
    check_array_like_object(svt, expected_class, a0, strict=strict)
    if (strict)
        expect_identical(svt, as(a0, expected_class))
}

check_NaArray_object <- function(naa, a0, strict=TRUE)
{
    Class0 <- S4Vectors:::capitalize(class(a0)[[1]])
    expected_class <- paste0("Na", Class0)
    check_array_like_object(naa, expected_class, a0, strict=strict)
    if (strict)
        expect_identical(naa, as(a0, expected_class))
}

### Subsetting a 1D ordinary array always preserves the "dim" and "dimnames"
### attributes i.e. it always returns another 1D array. This is inconsistent
### with the multidimensional case where, for example, linear subsetting (i.e.
### subsetting by a numeric vector or matrix) returns an ordinary vector
### (atomic or list). At the root of the problem is the behavior of
### base::drop() on a 1D array (see drop_even_if_1D() function in S4Arrays,
### in file R/dim-tuning-utils.R). So we "fix" subsetting of a 1D ordinary
### array by passing the result of the subsetting operation thru
### S4Arrays:::drop_even_if_1D().
test_linear_subsetting <- function(a, object, Mindex0)
{
    stopifnot(is.array(a))
    Lindex0 <- Mindex2Lindex(Mindex0, dim(a))

    expected <- S4Arrays:::drop_even_if_1D(a[Mindex0])
    expect_identical(object[Lindex0], expected)
    expect_identical(object[Mindex0], expected)
    expect_identical(object[Lindex0 + 0.99], expected)
    expect_identical(object[Mindex0 + 0.99], expected)

    revLindex <- rev(Lindex0)
    revMindex <- Lindex2Mindex(revLindex, dim(a))
    expected <- S4Arrays:::drop_even_if_1D(a[revLindex])
    expect_identical(object[revLindex], expected)
    expect_identical(object[revMindex], expected)
    type(revLindex) <- "double"
    type(revMindex) <- "double"
    expect_identical(object[revLindex], expected)
    expect_identical(object[revMindex], expected)
    expect_identical(object[revLindex + 0.99], expected)
    expect_identical(object[revMindex + 0.99], expected)

    Lindex <- c(NA, Lindex0, NA, NA, Lindex0, NA, rev(Lindex0))
    ## Lindex2Mindex() and 'object[Mindex]' don't accept NAs at the moment.
    #Mindex <- Lindex2Mindex(Lindex, dim(a))
    expected <- S4Arrays:::drop_even_if_1D(a[Lindex])
    expect_identical(object[Lindex], expected)
    #expect_identical(object[Mindex], expected)
    Lindex[1L] <- NaN      # this coerces 'Lindex' to "numeric"
    #Mindex[1L, 1L] <- NaN  # this coerces 'Mindex' to "numeric"
    expect_identical(object[Lindex], expected)
    #expect_identical(object[Mindex], expected)
}

test_summarize_op1 <- function(a, object, op)
{
    FUN <- match.fun(op)
    expected <- FUN(as.vector(a))
    current <- FUN(object)
    expect_identical(current, expected)
}

test_summarize_op2 <- function(a, object, op)
{
    FUN <- match.fun(op)
    if (op %in% c("var", "sd") ||
        is.double(a) && op %in% c("sum", "prod", "mean"))
    {
        EXPECT_FUN <- expect_equal
    } else {
        EXPECT_FUN <- expect_identical
    }
    expected <- FUN(as.vector(a))
    current <- FUN(object)
    if (op == "prod" && is.integer(current))
        expected <- as.integer(expected)
    EXPECT_FUN(current, expected)
    expected <- FUN(as.vector(a), na.rm=TRUE)
    current <- FUN(object, na.rm=TRUE)
    if (op == "prod" && is.integer(current))
        expected <- as.integer(expected)
    EXPECT_FUN(current, expected)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### test_3D_colrowMinsMaxs()
###

.fix_simple_colMinsMaxs_result <- function(ans, x, dims)
{
    if (type(x) == "integer" && type(ans) == "double")
        type(ans) <- "integer"
    ans <- S4Arrays:::set_dim(ans, tail(dim(x), n=-dims))
    x_dimnames <- dimnames(x)
    if (!is.null(x_dimnames))
        dimnames(ans) <- tail(x_dimnames, n=-dims)
    S4Arrays:::drop_even_if_1D(ans)
}

.fix_simple_rowMinsMaxs_result <- function(ans, x, dims)
{
    if (type(x) == "integer" && type(ans) == "double")
        type(ans) <- "integer"
    ans <- S4Arrays:::set_dim(ans, head(dim(x), n=dims))
    x_dimnames <- dimnames(x)
    if (!is.null(x_dimnames))
        dimnames(ans) <- head(x_dimnames, n=dims)
    S4Arrays:::drop_even_if_1D(ans)
}

.simple_3D_colMins <- function(x, na.rm=FALSE, dims=1)
{
    stopifnot(length(dim(x)) == 3L, isSingleNumber(dims))
    if (dims == 1) {
        ans <- apply(x, MARGIN=3, colMins, na.rm=na.rm)
    } else if (dims == 2) {
        ans <- suppressWarnings(apply(x, MARGIN=3, min, na.rm=na.rm))
    } else {
        stop("unsupported 'dims'")
    }
    .fix_simple_colMinsMaxs_result(ans, x, dims)
}

.simple_3D_colMaxs <- function(x, na.rm=FALSE, dims=1)
{
    stopifnot(length(dim(x)) == 3L, isSingleNumber(dims))
    if (dims == 1) {
        ans <- apply(x, MARGIN=3, colMaxs, na.rm=na.rm)
    } else if (dims == 2) {
        ans <- suppressWarnings(apply(x, MARGIN=3, max, na.rm=na.rm))
    } else {
        stop("unsupported 'dims'")
    }
    .fix_simple_colMinsMaxs_result(ans, x, dims)
}

.simple_3D_rowMins <- function(x, na.rm=FALSE, dims=1)
{
    stopifnot(length(dim(x)) == 3L, isSingleNumber(dims))
    if (dims == 1) {
        ans <- suppressWarnings(apply(x, MARGIN=1, min, na.rm=na.rm))
    } else if (dims == 2) {
        ans <- apply(x, MARGIN=2, rowMins, na.rm=na.rm)
    } else {
        stop("unsupported 'dims'")
    }
    .fix_simple_rowMinsMaxs_result(ans, x, dims)
}

.simple_3D_rowMaxs <- function(x, na.rm=FALSE, dims=1)
{
    stopifnot(length(dim(x)) == 3L, isSingleNumber(dims))
    if (dims == 1) {
        ans <- suppressWarnings(apply(x, MARGIN=1, max, na.rm=na.rm))
    } else if (dims == 2) {
        ans <- apply(x, MARGIN=2, rowMaxs, na.rm=na.rm)
    } else {
        stop("unsupported 'dims'")
    }
    .fix_simple_rowMinsMaxs_result(ans, x, dims)
}

### Tests *Mins() and *Maxs() methods on a 3D array-like object.
test_3D_colrowMinsMaxs <- function(object)
{
    a <- as.array(object)

    ## Base R does NOT allow an ordinay array with dimensions of extent
    ## zero to carry a character(0) in its dimnames, only a NULL. This means
    ## that, if 'object' is an Array derivative, then 'dimnames(object)'
    ## and 'dimnames(as.array(object))' are not guaranteed to be identical
    ## because a character(0) in the former will be replaced with a NULL in
    ## the latter.
    ## As a consequence, calling colMins/Maxs() or rowMins/Maxs()
    ## on 'object' won't necessarily produce the exact same result as
    ## calling .simple_3D_colMins/Maxs3D() or .simple_3D_rowMins/Maxs3D()
    ## on 'as.array(object)' when the result is a vector of length 0.
    ## More precisely, one can be named while the other is not.

    ## Does not look at the names if 'a1' and 'a2' are vectors of length 0.
    expect_almost_identical <- function(a1, a2) {
        if (is.vector(a1) && length(a1) == 0L) {
            a1 <- unname(a1)
            a2 <- unname(a2)
        }
        expect_identical(a1, a2)
    }

    ## dims == 1 (default)

    expected <- .simple_3D_colMins(a)
    expect_almost_identical(colMins(object), expected)
    expected <- .simple_3D_colMaxs(a)
    expect_almost_identical(colMaxs(object), expected)

    expected <- .simple_3D_rowMins(a)
    expect_almost_identical(rowMins(object), expected)
    expected <- .simple_3D_rowMaxs(a)
    expect_almost_identical(rowMaxs(object), expected)

    expected <- .simple_3D_colMins(a, na.rm=TRUE)
    expect_almost_identical(colMins(object, na.rm=TRUE), expected)
    expected <- .simple_3D_colMaxs(a, na.rm=TRUE)
    expect_almost_identical(colMaxs(object, na.rm=TRUE), expected)

    expected <- .simple_3D_rowMins(a, na.rm=TRUE)
    expect_almost_identical(rowMins(object, na.rm=TRUE), expected)
    expected <- .simple_3D_rowMaxs(a, na.rm=TRUE)
    expect_almost_identical(rowMaxs(object, na.rm=TRUE), expected)

    ## dims == 2

    expected <- .simple_3D_colMins(a, dims=2)
    expect_almost_identical(colMins(object, dims=2), expected)
    expected <- .simple_3D_colMaxs(a, dims=2)
    expect_almost_identical(colMaxs(object, dims=2), expected)

    expected <- .simple_3D_rowMins(a, dims=2)
    expect_almost_identical(rowMins(object, dims=2), expected)
    expected <- .simple_3D_rowMaxs(a, dims=2)
    expect_almost_identical(rowMaxs(object, dims=2), expected)

    expected <- .simple_3D_colMins(a, na.rm=TRUE, dims=2)
    expect_almost_identical(colMins(object, na.rm=TRUE, dims=2), expected)
    expected <- .simple_3D_colMaxs(a, na.rm=TRUE, dims=2)
    expect_almost_identical(colMaxs(object, na.rm=TRUE, dims=2), expected)

    expected <- .simple_3D_rowMins(a, na.rm=TRUE, dims=2)
    expect_almost_identical(rowMins(object, na.rm=TRUE, dims=2), expected)
    expected <- .simple_3D_rowMaxs(a, na.rm=TRUE, dims=2)
    expect_almost_identical(rowMaxs(object, na.rm=TRUE, dims=2), expected)
}

