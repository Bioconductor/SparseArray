
### When called on an ordinary matrix 'm' that contains a mix of NAs and NaNs,
### 'crossprod(m)' can return a square matrix 'cp' where the NA/NaN pattern
### is not symetric with respect to the diagonal. More precisely, there might
### be some row/col indices for which 'cp[i, j]' is NA_real_ and 'cp[j, i]'
### is NaN. This utility fixes that by setting the latter to NA_real_.
### This ensures that 'cp' is strictly symetric.
.fix_sym_mat_NA_NaN_pattern <- function(cp)
{
    is_NaN <- is.nan(cp)
    is_not_sym <- is_NaN != t(is_NaN)
    stopifnot(all(is.na(cp[is_not_sym])))
    cp[is_not_sym] <- NA_real_
    stopifnot(identical(cp, t(cp)))
    cp
}

.test_sym_crossprod_SparseMatrix <- function(m, svt)
{
    stopifnot(is.matrix(m),
              is(svt, "SparseMatrix"),
              identical(dim(m), dim(svt)))

    ## crossprod()
    expected <- crossprod(m)
    if (type(m) == "double") {
        expected <- .fix_sym_mat_NA_NaN_pattern(expected)
        EXPECT_FUN <- expect_equal
    } else {
        EXPECT_FUN <- expect_identical
    }
    cp <- crossprod(svt)
    EXPECT_FUN(cp, expected)
    EXPECT_FUN(cp, t(cp))
    EXPECT_FUN(cp, crossprod(svt, svt))
    EXPECT_FUN(cp, crossprod(svt, m))
    EXPECT_FUN(cp, crossprod(m, svt))

    ## tcrossprod()
    expected <- tcrossprod(m)
    if (type(m) == "double")
        expected <- .fix_sym_mat_NA_NaN_pattern(expected)
    tcp <- tcrossprod(svt)
    EXPECT_FUN(tcp, expected)
    EXPECT_FUN(tcp, t(tcp))
    EXPECT_FUN(tcp, tcrossprod(svt, svt))
    EXPECT_FUN(tcp, tcrossprod(svt, m))
    EXPECT_FUN(tcp, tcrossprod(m, svt))
}

test_that("crossprod()/tcrossprod() on input objects of type \"double\"", {
    m0 <- matrix(0, nrow=5, ncol=3)
    m0[3, 1] <- Inf
    m0[2, 3] <- -11.99
    svt0 <- as(m0, "SVT_SparseMatrix")
    .test_sym_crossprod_SparseMatrix(m0, svt0)

    m1 <- matrix(c(0, -4.5, 7, NA, 0, NaN, Inf, -Inf), nrow=1)
    svt1 <- as(m1, "SVT_SparseMatrix")
    .test_sym_crossprod_SparseMatrix(m1, svt1)

    m2 <- matrix(0, nrow=6, ncol=4, dimnames=list(letters[1:6], LETTERS[1:4]))
    m2[c(24, 1:2, 8, 10, 15:17)] <- 1:8 - 3.5
    svt2 <- as(m2, "SVT_SparseMatrix")
    .test_sym_crossprod_SparseMatrix(m2, svt2)

    m3 <- matrix(0, nrow=6, ncol=7,
                 dimnames=list(letters[21:26], LETTERS[20:26]))
    m3[3 + 4*(0:9)] <- 2.4^(1:10)
    m3[4 + 4*(0:9)] <- -(101:110)
    m3[1, 5] <- NaN
    m3[5, 3] <- Inf
    svt3 <- as(m3, "SVT_SparseMatrix")
    .test_sym_crossprod_SparseMatrix(m3, svt3)

    expected <- crossprod(m2, m3)
    cp <- crossprod(svt2, svt3)
    expect_identical(cp, expected)
    expect_identical(cp, crossprod(svt2, m3))
    expect_identical(cp, crossprod(m2, svt3))
    expect_identical(crossprod(svt3, svt2), t(expected))
    tm2 <- t(m2)
    tm3 <- t(m3)
    tsvt2 <- t(svt2)
    tsvt3 <- t(svt3)
    expected <- tcrossprod(tm2, tm3)
    tcp <- tcrossprod(tsvt2, tsvt3)
    expect_identical(tcp, expected)
    expect_identical(tcp, tcrossprod(tsvt2, tm3))
    expect_identical(tcp, tcrossprod(tm2, tsvt3))
    expect_identical(tcrossprod(tsvt3, tsvt2), t(expected))

    ## Zero rows.
    m4 <- matrix(double(0), ncol=3, dimnames=list(NULL, letters[1:3]))
    svt4 <- as(m4, "SVT_SparseMatrix")
    .test_sym_crossprod_SparseMatrix(m4, svt4)

    ## Zero cols.
    m5 <- matrix(double(0), nrow=6, dimnames=list(LETTERS[1:6], NULL))
    svt5 <- as(m5, "SVT_SparseMatrix")
    .test_sym_crossprod_SparseMatrix(m5, svt5)

    expected <- crossprod(m3, m5)
    cp <- crossprod(svt3, svt5)
    expect_identical(cp, expected)
    expect_identical(cp, crossprod(svt3, m5))
    expect_identical(cp, crossprod(m3, svt5))
    expect_identical(crossprod(svt5, svt3), t(expected))
    tm5 <- t(m5)
    tsvt5 <- t(svt5)
    expected <- tcrossprod(tm3, tm5)
    tcp <- tcrossprod(tsvt3, tsvt5)
    expect_identical(tcp, expected)
    expect_identical(tcp, tcrossprod(tsvt3, tm5))
    expect_identical(tcp, tcrossprod(tm3, tsvt5))
    expect_identical(tcrossprod(tsvt5, tsvt3), t(expected))
})

.test_crossprod_int_SVT_SparseMatrix <- function(m1, m2=NULL)
{
    stopifnot(is.matrix(m1))
    svt1 <- as(m1, "SVT_SparseMatrix")
    stopifnot(type(svt1) == "integer")
    expected <- crossprod(m1)
    cp <- crossprod(svt1)
    expect_identical(cp, expected)
    expect_identical(cp, t(cp))
    expect_identical(cp, crossprod(svt1, svt1))
    svt <- `type<-`(svt1, "double")
    expect_identical(cp, crossprod(svt))
    expect_identical(cp, crossprod(svt, svt))

    if (is.null(m2))
        return()

    stopifnot(is.matrix(m2))
    svt2 <- as(m2, "SVT_SparseMatrix")
    stopifnot(type(svt2) == "integer")
    expected <- crossprod(m2)
    cp <- crossprod(svt2)
    expect_identical(cp, expected)
    expect_identical(cp, t(cp))
    expect_identical(cp, crossprod(svt2, svt2))
    svt <- `type<-`(svt2, "double")
    expect_identical(cp, crossprod(svt))
    expect_identical(cp, crossprod(svt, svt))

    expected <- crossprod(m1, m2)
    cp <- crossprod(svt1, svt2)
    expect_identical(cp, expected)
    expect_identical(crossprod(svt2, svt1), t(expected))
}

test_that("crossprod() on input objects of type \"integer\"", {
    m1 <- matrix(c(0L, -4L, 7L, NA_integer_, 0L, NA_integer_), nrow=1)
    .test_crossprod_int_SVT_SparseMatrix(m1)

    m2 <- matrix(0L, nrow=6, ncol=4, dimnames=list(letters[1:6], LETTERS[1:4]))
    m2[c(24, 1:2, 8, 10, 15:17)] <- (1:8) * 10L - 35L

    m3 <- matrix(0L, nrow=6, ncol=7,
                 dimnames=list(letters[21:26], LETTERS[20:26]))
    m3[3 + 4*(0:9)] <- 1:10
    m3[4 + 4*(0:9)] <- -(101:110)
    .test_crossprod_int_SVT_SparseMatrix(m2, m3)

    m2[2, 4] <- NA_integer_
    m3[1, 5] <- NA_integer_
    .test_crossprod_int_SVT_SparseMatrix(m2, m3)

    ## Zero rows.
    m4 <- matrix(integer(0), ncol=3, dimnames=list(NULL, letters[1:3]))
    .test_crossprod_int_SVT_SparseMatrix(m4)

    ## Zero cols.
    m5 <- matrix(integer(0), nrow=6, dimnames=list(LETTERS[1:6], NULL))
    .test_crossprod_int_SVT_SparseMatrix(m5)
    .test_crossprod_int_SVT_SparseMatrix(m3, m5)
})

test_that("matrix multiplication", {
    m1 <- matrix(0L, nrow=15, ncol=6)
    m1[c(2, 6, 12:17, 22:33, 55, 59:62, 90)] <- 101:126
    svt1 <- as(m1, "SVT_SparseMatrix")

    set.seed(333)
    m2 <- matrix(runif(12), nrow=6)
    svt2 <- as(m2, "SVT_SparseMatrix")

    expected <- m1 %*% m2
    expect_identical(svt1 %*% svt2, expected)
    expect_identical(svt1 %*% m2, expected)
    expect_identical(m1 %*% svt2, expected)
    expect_identical(t(m1) %*% m1, crossprod(m1))
    expect_identical(m1 %*% t(m1), tcrossprod(m1))
})

