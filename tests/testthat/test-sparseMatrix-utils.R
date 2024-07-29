
test_that("CsparseMatrix() and RsparseMatrix()", {
    CsparseMatrix <- SparseArray:::CsparseMatrix
    RsparseMatrix <- SparseArray:::RsparseMatrix

    i <- c(1:3, 1L, 10:7)
    j <- c(1:3, 1L, 7:4)
    nzdata <- c(11:14, NA, 0:2)
    m <- matrix(0, nrow=10, ncol=7)
    m[cbind(i, j)] <- nzdata

    x <- CsparseMatrix(dim(m), i, j, nzdata)
    expect_true(is(x, "dgCMatrix"))
    expect_true(validObject(x))
    expect_identical(as.matrix(x), m)

    x <- RsparseMatrix(dim(m), i, j, nzdata)
    expect_true(is(x, "dgRMatrix"))
    expect_true(validObject(x))
    expect_identical(as.matrix(x), m)
})

test_that("colStats_dgCMatrix", {
    colMins_dgCMatrix   <- SparseArray:::colMins_dgCMatrix
    colMaxs_dgCMatrix   <- SparseArray:::colMaxs_dgCMatrix
    colRanges_dgCMatrix <- SparseArray:::colRanges_dgCMatrix
    colVars_dgCMatrix   <- SparseArray:::colVars_dgCMatrix

    do_other_tests <- function(m, m0) {
        ## colMins_dgCMatrix()
        expected <- apply(m, 2, min)
        expect_identical(colMins(m), expected)
        expect_identical(colMins_dgCMatrix(m0), expected)
        expected <- suppressWarnings(apply(m, 2, min, na.rm=TRUE))
        expect_identical(colMins(m, na.rm=TRUE), expected)
        expect_identical(colMins_dgCMatrix(m0, na.rm=TRUE), expected)

        ## colMaxs_dgCMatrix()
        expected <- apply(m, 2, max)
        expect_identical(colMaxs(m), expected)
        expect_identical(colMaxs_dgCMatrix(m0), expected)
        expected <- suppressWarnings(apply(m, 2, max, na.rm=TRUE))
        expect_identical(colMaxs(m, na.rm=TRUE), expected)
        expect_identical(colMaxs_dgCMatrix(m0, na.rm=TRUE), expected)

        ## colRanges_dgCMatrix()
        expected <- t(apply(m, 2, range))
        expect_identical(colRanges(m), expected)
        expect_identical(colRanges_dgCMatrix(m0), expected)
        expected <- suppressWarnings(t(apply(m, 2, range, na.rm=TRUE)))
        expect_identical(colRanges(m, na.rm=TRUE), expected)
        expect_identical(colRanges_dgCMatrix(m0, na.rm=TRUE), expected)

        ## colVars_dgCMatrix()
        expected <- apply(m, 2, var)
        expect_equal(colVars(m), expected)
        expect_equal(colVars_dgCMatrix(m0), expected)
        expected <- apply(m, 2, var, na.rm=TRUE)
        expect_equal(colVars(m, na.rm=TRUE), expected)
        expect_equal(colVars_dgCMatrix(m0, na.rm=TRUE), expected)
    }

    set.seed(123)
    m0 <- rsparsematrix(22, 10, density=0.25)
    m0[3L, 1L] <- NA
    m0[3L, 2L] <- NaN
    m0[4L, 4L] <- NA
    m0[8L, 4L] <- NaN
    m0[4L, 8L] <- NaN
    m0[22L, 8L] <- NA
    m0[ , 9L] <- NA
    m0[ , 10L] <- NaN
    m <- as.matrix(m0)

    ## rowsum()
    rgroup <- sample(9, nrow(m0), replace=TRUE)  # define groups of rows
    expect_identical(rowsum(m0, rgroup), rowsum(m, rgroup))
    expect_identical(rowsum(m0, rgroup, reorder=FALSE),
                     rowsum(m, rgroup, reorder=FALSE))
    expect_identical(rowsum(m0, rgroup, na.rm=TRUE),
                     rowsum(m, rgroup, na.rm=TRUE))

    ## colsum()
    cgroup <- sample(5, ncol(m0), replace=TRUE)  # define groups of cols
    expect_equal(colsum(m0, cgroup), colsum(m, cgroup))
    expect_equal(colsum(m0, cgroup, reorder=FALSE),
                 colsum(m, cgroup, reorder=FALSE))
    expect_equal(colsum(m0, cgroup, na.rm=TRUE),
                 colsum(m, cgroup, na.rm=TRUE))

    do_other_tests(m, m0)
})

