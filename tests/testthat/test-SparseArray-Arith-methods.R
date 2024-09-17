
.test_Arith_SVT1_v2 <- function(a1, svt1, v2, relax.MOD.and.IDIV=FALSE)
{
    if (is.infinite(v2)) {
        expect_error(svt1 * v2, "not supported")
        expect_error(v2 * svt1, "not supported")
    } else {
        a <- a1 * v2
        svt <- svt1 * v2
        check_SVT_SparseArray_object(svt, a)
        svt <- v2 * svt1
        check_SVT_SparseArray_object(svt, a)
    }

    if (v2 == 0)
        return()

    a <- a1 / v2
    svt <- svt1 / v2
    check_SVT_SparseArray_object(svt, a)

    if (v2 > 0) {
        a <- a1 ^ v2
        svt <- svt1 ^ v2
        expect_true(is(svt, "SVT_SparseArray"))
        expect_true(validObject(svt))
        expect_equal(as.array(svt), a)
    }

    if (relax.MOD.and.IDIV) {
        reconstructed <- as.array(svt1 %% v2) + v2 * as.array(svt1 %/% v2)
        expect_equal(reconstructed, a1)
    } else {
        a <- a1 %% v2
        svt <- svt1 %% v2
        expect_true(is(svt, "SVT_SparseArray"))
        expect_true(validObject(svt))
        expect_equal(as.array(svt), a)

        a <- a1 %/% v2
        svt <- svt1 %/% v2
        check_SVT_SparseArray_object(svt, a)
    }
}

.test_Arith_SVT1_SVT2 <- function(a1, a2, svt1, svt2)
{
    a <- a1 + a2
    svt <- svt1 + svt2
    check_SVT_SparseArray_object(svt, a)
    expect_identical(svt1 + a2, svt)
    expect_identical(a1 + svt2, svt)

    a <- a1 - a2
    svt <- svt1 - svt2
    check_SVT_SparseArray_object(svt, a)
    expect_identical(svt1 - a2, svt)
    expect_identical(a1 - svt2, svt)

    a <- a1 * a2
    svt <- svt1 * svt2
    check_SVT_SparseArray_object(svt, a)
    expect_identical(svt1 * a2, svt)
    expect_identical(a1 * svt2, svt)
}

test_that("'Arith' ops between SVT_SparseArray object and single value", {

    ## --- 3D ---

    a1 <- make_3D_integer_array()
    svt1 <- as(a1, "SVT_SparseArray")

    .test_Arith_SVT1_v2(a1, svt1, 0L)
    .test_Arith_SVT1_v2(a1, svt1, 5L)
    .test_Arith_SVT1_v2(a1, svt1, -5L)
    .test_Arith_SVT1_v2(a1, svt1, 1L)
    .test_Arith_SVT1_v2(a1, svt1, -1L)
    .test_Arith_SVT1_v2(a1, svt1, 5.1)
    .test_Arith_SVT1_v2(a1, svt1, -5.1, relax.MOD.and.IDIV=TRUE)
    .test_Arith_SVT1_v2(a1, svt1, 0.001, relax.MOD.and.IDIV=TRUE)
    .test_Arith_SVT1_v2(a1, svt1, -0.001, relax.MOD.and.IDIV=TRUE)
    .test_Arith_SVT1_v2(a1, svt1, 134)
    expect_warning(svt1 * 10650000L, "integer overflow")
    expect_warning(10650000L * svt1, "integer overflow")

    a0 <- a1[ , 0, ]
    svt0 <- as(a0, "SVT_SparseArray")
    .test_Arith_SVT1_v2(a0, svt0, 5L)

    ## --- 2D ---

    m1 <- a1[ , , 2]
    svt1 <- as(m1, "SVT_SparseMatrix")

    .test_Arith_SVT1_v2(m1, svt1, 0L)
    .test_Arith_SVT1_v2(m1, svt1, 5L)
    .test_Arith_SVT1_v2(m1, svt1, -5L)
    .test_Arith_SVT1_v2(m1, svt1, 1L)
    .test_Arith_SVT1_v2(m1, svt1, -1L)

    m2 <- matrix(1:20, nrow=5)
    m2[2, ] <- 0L
    m2[3:5, 1] <- c(18L, -18L, -15L)
    m2[3, 2:3] <- 0L
    svt2 <- as(m2, "SVT_SparseMatrix")

    .test_Arith_SVT1_v2(m2, svt2, 0L)
    .test_Arith_SVT1_v2(m2, svt2, 5L)
    .test_Arith_SVT1_v2(m2, svt2, -5L)
    .test_Arith_SVT1_v2(m2, svt2, 1L)
    .test_Arith_SVT1_v2(m2, svt2, -1L)

    m2[1, ] <- m2[1, ] + 0.5
    m2[2, ] <- c(NA, NaN, Inf, -Inf)
    m2[3, ] <- c(-1, -0.1, 0.1, 1)
    svt2 <- as(m2, "SVT_SparseMatrix")

    .test_Arith_SVT1_v2(m2, svt2, 0L)
    .test_Arith_SVT1_v2(m2, svt2, 5L)
    .test_Arith_SVT1_v2(m2, svt2, -5L)
    .test_Arith_SVT1_v2(m2, svt2, 1L)
    .test_Arith_SVT1_v2(m2, svt2, -1L)
    .test_Arith_SVT1_v2(m2, svt2, Inf)
    .test_Arith_SVT1_v2(m2, svt2, -Inf)

    ## --- Not expected to work ---

    expect_error(5 / svt2, "not supported")
    expect_error(5 ^ svt2, "not supported")
    expect_error(svt2 + "A", "not supported")
    expect_error(svt2 + 1, "not supported")
    expect_error(svt2 - 1, "not supported")
    expect_error(svt2 * 1:2, "not supported")
    expect_error(svt2 * NA_integer_, "not supported")
    expect_error(svt2 * NA_real_, "not supported")
    expect_error(svt2 * NaN, "not supported")
    expect_error(svt2 * Inf, "not supported")
    expect_error(svt2 * -Inf, "not supported")
    expect_error(svt2 ^ 0, "not supported")
    expect_error(svt2 ^ -2, "not supported")
    expect_error(svt2 / 0, "not supported")
    expect_error(svt2 %% 0, "not supported")
    expect_error(svt2 %/% 0, "not supported")
})

test_that("'Arith' ops between 2 SVT_SparseArray objects", {
    a1 <- a2 <- array(0L, 6:4)
    dimnames(a1) <- list(letters[1:6], NULL, LETTERS[1:4])
    dimnames(a2) <- list(NULL, letters[22:26], LETTERS[23:26])

    svt1 <- as(a1, "SVT_SparseArray")
    svt2 <- as(a2, "SVT_SparseArray")

    .test_Arith_SVT1_SVT2(a1, a2, svt1, svt2)
    .test_Arith_SVT1_SVT2(a2, a1, svt2, svt1)

    a1[c(2:3, 6), 2, 1] <- 101:103
    a2[2:4, 2, 1] <- 1001:1003
    a1[c(1, 6), 1 , 2] <- 201:202
    a2[c(3, 5), 2 , 2] <- 2001:2002
    a1[1:5, 5, 3] <- 301:305
    a2[c(2:4, 6), 5, 3] <- 3001:3004
    a2[6, 5, 4] <- NA
    svt1 <- as(a1, "SVT_SparseArray")
    svt2 <- as(a2, "SVT_SparseArray")

    .test_Arith_SVT1_SVT2(a1, a2, svt1, svt2)
    .test_Arith_SVT1_SVT2(a2, a1, svt2, svt1)
    .test_Arith_SVT1_SVT2(a1, a1, svt1, svt1)
    .test_Arith_SVT1_SVT2(a2, a2, svt2, svt2)
    expect_identical((svt1 + svt2) - svt1, `dimnames<-`(svt2, dimnames(a1)))
    expect_identical(svt1 + (svt2 - svt1), `dimnames<-`(svt2, dimnames(a1)))
    expect_identical(svt1 - (svt1 - svt2), `dimnames<-`(svt2, dimnames(a1)))
    .test_Arith_SVT1_SVT2(a1, a1 + a2, svt1, svt1 + svt2)
    .test_Arith_SVT1_SVT2(a1, a2 - a1, svt1, svt2 - svt1)
    .test_Arith_SVT1_SVT2(a2, a1 - a2, svt2, svt1 - svt2)
    .test_Arith_SVT1_SVT2(a1, a1 * a2, svt1, svt1 * svt2)

    dimnames(a1) <- dimnames(svt1) <- NULL
    .test_Arith_SVT1_SVT2(a1, a2, svt1, svt2)
    .test_Arith_SVT1_SVT2(a2, a1, svt2, svt1)

    a1[ , , ] <- 0L
    svt1 <- as(a1, "SVT_SparseArray")

    .test_Arith_SVT1_SVT2(a1, a2, svt1, svt2)
    .test_Arith_SVT1_SVT2(a2, a1, svt2, svt1)
    .test_Arith_SVT1_SVT2(a1, a1, svt1, svt1)
    .test_Arith_SVT1_SVT2(a2, a2, svt2, svt2)

    a1[ , , ] <- 0.0
    svt1 <- as(a1, "SVT_SparseArray")

    .test_Arith_SVT1_SVT2(a1, a2, svt1, svt2)
    .test_Arith_SVT1_SVT2(a2, a1, svt2, svt1)
    .test_Arith_SVT1_SVT2(a1, a1, svt1, svt1)
    .test_Arith_SVT1_SVT2(a2, a2, svt2, svt2)

    a0 <- a2[ , 0, ]
    svt0 <- as(a0, "SVT_SparseArray")
    .test_Arith_SVT1_SVT2(a0, a0, svt0, svt0)
    .test_Arith_SVT1_SVT2(a0, unname(a0), svt0, unname(svt0))
    .test_Arith_SVT1_SVT2(unname(a0), a0, unname(svt0), svt0)

    ## --- Not expected to work ---

    expect_error(svt1 + svt2[ , , -1], "non-conformable")
    expect_error(svt1 - svt2[ , , -1], "non-conformable")
    expect_error(svt1 * svt2[ , , -1], "non-conformable")
    expect_error(svt1 / svt2, "not supported")
    expect_error(svt1 ^ svt2, "not supported")
    expect_error(svt1 %% svt2, "not supported")
    expect_error(svt1 %/% svt2, "not supported")
})

test_that("unary minus on a SVT_SparseArray object", {
    a1 <- make_3D_integer_array()
    svt1 <- as(a1, "SVT_SparseArray")
    check_SVT_SparseArray_object(- svt1, - a1)
    expect_identical(- (- svt1), svt1)

    a2 <- make_3D_double_array()
    svt2 <- as(a2, "SVT_SparseArray")
    check_SVT_SparseArray_object(- svt2, - a2)
    expect_identical(- (- svt2), svt2)

    a3 <- make_3D_complex_array()
    svt3 <- as(a3, "SVT_SparseArray")
    check_SVT_SparseArray_object(- svt3, - a3)
    expect_identical(- (- svt3), svt3)
})

