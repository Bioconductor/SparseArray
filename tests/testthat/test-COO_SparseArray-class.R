
test_that("nzwhich(), nzvals(), `nzvals<-`() on COO_SparseArray objects", {

    dim <- c(5L, 8L, 2L)
    nzcoo  <- Lindex2Mindex(c(5:12, 15, 3:1), dim)
    nzdata <- c(101:107, 0, 109:110, 0, 112)
    coo0 <- COO_SparseArray(dim, nzcoo, nzdata)
    a0 <- as.array(coo0)

    check_SparseArray_object(coo0, "COO_SparseArray", a0)
    expect_identical(nzwhich(coo0), nzwhich(a0))
    expect_identical(nzvals(coo0), nzvals(a0))

    coo <- coo0
    a <- a0

    nzvals(coo) <- nzvals(a) <- nzvals(a) + 0.99
    check_SparseArray_object(coo, "COO_SparseArray", a)
    expect_identical(nzwhich(coo0), nzwhich(a0))
    expect_identical(nzvals(coo0), nzvals(a0))

    nzvals(coo) <- nzvals(a) <- -0.1 * (11:15)
    check_SparseArray_object(coo, "COO_SparseArray", a)
    expect_identical(nzwhich(coo0), nzwhich(a0))
    expect_identical(nzvals(coo0), nzvals(a0))
})

