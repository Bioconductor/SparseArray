
test_that("crossprod() on input objects of type \"double\"", {
    m1 <- matrix(c(0, -4.5, 7, NA, 0, NaN, Inf, -Inf), nrow=1)
    svt1 <- as(m1, "SVT_SparseMatrix")

    expected <- crossprod(m1)
    cp <- crossprod(svt1)
    expect_equal(cp, expected)
    expect_identical(t(cp), cp)
    expect_identical(crossprod(svt1, svt1), cp)

    m2 <- matrix(0, nrow=6, ncol=4, dimnames=list(letters[1:6], LETTERS[1:4]))
    m2[c(24, 1:2, 8, 10, 15:17)] <- 1:8 - 3.5
    svt2 <- as(m2, "SVT_SparseMatrix")

    expected <- crossprod(m2)
    cp <- crossprod(svt2)
    expect_identical(cp, expected)
    expect_identical(t(cp), cp)
    expect_identical(crossprod(svt2, svt2), expected)

    m3 <- matrix(0, nrow=6, ncol=7,
		 dimnames=list(letters[21:26], LETTERS[20:26]))
    m3[3 + 4*(0:9)] <- 1:10
    m3[4 + 4*(0:9)] <- -(101:110)
    svt3 <- as(m3, "SVT_SparseMatrix")

    expected <- crossprod(m3)
    cp <- crossprod(svt3)
    expect_identical(cp, expected)
    expect_identical(t(cp), cp)
    expect_identical(crossprod(svt3, svt3), expected)

    expected <- crossprod(m2, m3)
    cp <- crossprod(svt2, svt3)
    expect_identical(cp, expected)
    expect_identical(crossprod(svt3, svt2), t(expected))

    ## Zero rows.
    m4 <- matrix(double(0), ncol=3, dimnames=list(NULL, letters[1:3]))
    svt4 <- as(m4, "SVT_SparseMatrix")

    expected <- crossprod(m4)
    cp <- crossprod(svt4)
    expect_identical(cp, expected)
    expect_identical(t(cp), cp)
    expect_identical(crossprod(svt4, svt4), expected)

    ## Zero cols.
    m5 <- matrix(double(0), nrow=6, dimnames=list(LETTERS[1:6], NULL))
    svt5 <- as(m5, "SVT_SparseMatrix")

    expected <- crossprod(m5)
    cp <- crossprod(svt5)
    expect_identical(cp, expected)
    expect_identical(t(cp), cp)
    expect_identical(crossprod(svt5, svt5), expected)

    expected <- crossprod(m3, m5)
    cp <- crossprod(svt3, svt5)
    expect_identical(cp, expected)
    expect_identical(crossprod(svt5, svt3), t(expected))
})


