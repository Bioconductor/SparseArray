
test_that("nzcount(), nzwhich() and nzvals() on CsparseMatrix objects", {
    ## With a "long object" (i.e. length(x) >= 2^31).
    set.seed(2009)
    x <- Matrix::rsparsematrix(1e5, 3e4, 0.0001)
    expect_equal(nzcount(x), 300000)
    nzidx <- nzwhich(x)
    expect_equal(length(nzidx), 300000)
    expect_false(is.unsorted(nzidx, strict=TRUE))
    expect_identical(nzidx, SparseArray:::default_nzwhich(x))
    nzvals0 <- x[nzidx]
    expect_true(all(nzvals0 != 0))
    expect_identical(nzvals(x), nzvals0)
})

