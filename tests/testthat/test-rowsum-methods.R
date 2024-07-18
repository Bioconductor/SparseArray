
### Set the rownames of a zero-row matrix to NULL.
### When the input is a zero-row matrix, base::rowsum() returns a zero-row
### matrix where the rownames are set to character(0). Note that this is very
### unusual to see on an ordinary matrix where explicitly trying to set the
### rownames to character(0) (using the rownames() or dimnames() setter) will
### always set them to NULL instead.
.fix_rownames <- function(x)
{
    if (nrow(x) == 0L)
        rownames(x) <- NULL
    x
}

.test_rowsum_methods <- function(m, group, FUN=rowsum)
{
    stopifnot(is.matrix(m))
    FUN <- match.fun(FUN)
    svt <- as(m, "SVT_SparseMatrix")
    dgcm <- as(m, "dgCMatrix")
    coo <- as(svt, "COO_SparseMatrix")

    check_rs1_rs2 <- function(expected, rs1, rs2) {
        expected <- .fix_rownames(expected)
        expect_true(is.matrix(rs1))
        expect_identical(typeof(rs1), typeof(expected))
        if (typeof(expected) == "double") {
            expect_equal(rs1, expected)
	} else {
            expect_identical(rs1, expected)
        }
        expect_true(is.matrix(rs2))
        expect_equal(rs2, expected)  # 'rs2' always of type "double"
    }

    expected <- FUN(m, group)
    rs1 <- FUN(svt, group)
    rs2 <- FUN(dgcm, group)
    check_rs1_rs2(expected, rs1, rs2)
    expect_identical(FUN(coo, group), rs1)

    expected <- FUN(m, group, na.rm=TRUE)
    rs1 <- FUN(svt, group, na.rm=TRUE)
    rs2 <- FUN(dgcm, group, na.rm=TRUE)
    check_rs1_rs2(expected, rs1, rs2)
    expect_identical(FUN(coo, group, na.rm=TRUE), rs1)

    expected <- FUN(m, group, reorder=FALSE)
    rs1 <- FUN(svt, group, reorder=FALSE)
    rs2 <- FUN(dgcm, group, reorder=FALSE)
    check_rs1_rs2(expected, rs1, rs2)
    expect_identical(FUN(coo, group, reorder=FALSE), rs1)

    expected <- FUN(m, group, reorder=FALSE, na.rm=TRUE)
    rs1 <- FUN(svt, group, reorder=FALSE, na.rm=TRUE)
    rs2 <- FUN(dgcm, group, reorder=FALSE, na.rm=TRUE)
    check_rs1_rs2(expected, rs1, rs2)
    expect_identical(FUN(coo, group, reorder=FALSE, na.rm=TRUE), rs1)
}

test_that("rowsum()/colsum() on a SVT_SparseMatrix or dgCMatrix object", {
    ## type "double"
    m1 <- matrix(0, nrow=6, ncol=4)
    group <- c("B", "A", "B", "B", "B", "A")
    colnames(m1) <- letters[1:4]
    .test_rowsum_methods(m1, group)
    .test_rowsum_methods(t(m1), group, FUN=colsum)

    m1[ , 1] <- c(8.55, Inf, NA_real_, 0, NaN, -Inf)
    m1[ , 3] <- c(0.6, -11.99, 0, 4.44, 0, 0)
    m1[ , 4] <- 1:6
    .test_rowsum_methods(m1, group)
    .test_rowsum_methods(t(m1), group, FUN=colsum)
    .test_rowsum_methods(m1[0L,   ], integer(0))
    .test_rowsum_methods(t(m1[0L,   ]), integer(0), FUN=colsum)
    .test_rowsum_methods(m1[  , 0L], group)
    .test_rowsum_methods(t(m1[  , 0L]), group, FUN=colsum)
    .test_rowsum_methods(m1[0L, 0L], integer(0))
    .test_rowsum_methods(t(m1[0L, 0L]), integer(0), FUN=colsum)

    ## type "integer"
    m2 <- matrix(0L, nrow=6, ncol=4)
    dimnames(m2) <- list(letters[21:26], letters[1:4])
    m2[1, 2] <- NA_integer_
    m2[3, 2] <- 99L
    m2[ , 4] <- 1:6
    .test_rowsum_methods(m2, group)
    .test_rowsum_methods(t(m2), group, FUN=colsum)
})

