make_3D_integer_array <- function()
{
    a <- array(0L, 6:4,
               dimnames=list(letters[1:6], NULL, LETTERS[1:4]))
    a[c(2:3, 6), 2, 1] <- 101:103
    a[c(1, 6), 1 , 2] <- 201:202
    a[c(3:4, 6), c(2, 5), 2] <- a[ , 3, 2] <- a[6, 4, 2] <- 1L
    a[4, 5, 2] <- NA
    a[1:5, 5, 3] <- -(301:305)
    a[6, 5, 4] <- NA
    a
}

make_3D_double_array <- function()
{
    a <- array(0.0, c(7, 10, 3),
               dimnames=list(NULL, letters[1:10], LETTERS[1:3]))
    set.seed(123)
    a[5*(1:26)] <- runif(26, min=-5, max=10)
    a[2, c(1:4, 7:9), 1] <- c(NA, NaN, Inf, 3e9, 256, -0.999, -1)
    a
}

make_3D_complex_array <- function()
{
    a <- make_3D_integer_array()
    set.seed(123)
    a[6*(1:20)] <- a[6*(1:20)] + runif(20, min=-5, max=10) * 1i
    a[2, , 4] <- c(NA, NaN,  Inf,    3e9, 256)
    a[3, , 4] <- c(NA, NaN, -Inf, -0.999,  -1) * 1i
    a
}

check_SparseArray_object <- function(object, expected_class, a0)
{
    expect_true(class(object) == expected_class)
    expect_true(validObject(object))
    expect_identical(dim(object), dim(a0))
    expect_identical(dimnames(object), dimnames(a0))
    expect_identical(type(object), type(a0))
    expect_identical(as.array(object), a0)
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

test_summarize_op1 <- function(a, svt, op)
{
    FUN <- match.fun(op)
    expected <- FUN(as.vector(a))
    current <- FUN(svt)
    expect_identical(current, expected)
}

test_summarize_op2 <- function(a, svt, op)
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
    current <- FUN(svt)
    if (op == "prod" && is.integer(current))
        expected <- as.integer(expected)
    EXPECT_FUN(current, expected)
    expected <- FUN(as.vector(a), na.rm=TRUE)
    current <- FUN(svt, na.rm=TRUE)
    if (op == "prod" && is.integer(current))
        expected <- as.integer(expected)
    EXPECT_FUN(current, expected)
}

