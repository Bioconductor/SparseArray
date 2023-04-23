.test_Arith_SVT1_v2 <- function(a1, svt1, v2, relax.MOD.and.IDIV=FALSE)
{
    expected <- a1 * v2
    svt <- svt1 * v2
    expect_true(is(svt, "SVT_SparseArray"))
    expect_true(validObject(svt))
    expect_identical(as.array(svt), expected)
    svt <- v2 * svt1
    expect_true(is(svt, "SVT_SparseArray"))
    expect_true(validObject(svt))
    expect_identical(as.array(svt), expected)

    if (v2 == 0)
        return()

    expected <- a1 / v2
    svt <- svt1 / v2
    expect_true(is(svt, "SVT_SparseArray"))
    expect_true(validObject(svt))
    expect_identical(as.array(svt), expected)

    if (v2 > 0) {
        expected <- a1 ^ v2
        svt <- svt1 ^ v2
        expect_true(is(svt, "SVT_SparseArray"))
        expect_true(validObject(svt))
        expect_identical(as.array(svt), expected)
    }

    if (relax.MOD.and.IDIV) {
        reconstructed <- as.array(svt1 %% v2) + v2 * as.array(svt1 %/% v2)
        expect_equal(reconstructed, a1)
    } else {
        expected <- a1 %% v2
        svt <- svt1 %% v2
        expect_true(is(svt, "SVT_SparseArray"))
        expect_true(validObject(svt))
        expect_equal(as.array(svt), expected)

        expected <- a1 %/% v2
        svt <- svt1 %/% v2
        expect_true(is(svt, "SVT_SparseArray"))
        expect_true(validObject(svt))
        expect_identical(as.array(svt), expected)
    }
}

.test_Arith_SVT1_SVT2 <- function(a1, a2, svt1, svt2)
{
    expected <- a1 + a2
    svt <- svt1 + svt2
    expect_true(is(svt, "SVT_SparseArray"))
    expect_true(validObject(svt))
    expect_identical(as.array(svt), expected)
    svt <- svt1 + a2
    expect_true(is(svt, "SVT_SparseArray"))
    expect_true(validObject(svt))
    expect_identical(as.array(svt), expected)
    svt <- a1 + svt2
    expect_true(is(svt, "SVT_SparseArray"))
    expect_true(validObject(svt))
    expect_identical(as.array(svt), expected)

    expected <- a1 - a2
    svt <- svt1 - svt2
    expect_true(is(svt, "SVT_SparseArray"))
    expect_true(validObject(svt))
    expect_identical(as.array(svt), expected)
    svt <- svt1 - a2
    expect_true(is(svt, "SVT_SparseArray"))
    expect_true(validObject(svt))
    expect_identical(as.array(svt), expected)
    svt <- a1 - svt2
    expect_true(is(svt, "SVT_SparseArray"))
    expect_true(validObject(svt))
    expect_identical(as.array(svt), expected)

    expected <- a1 * a2
    svt <- svt1 * svt2
    expect_true(is(svt, "SVT_SparseArray"))
    expect_true(validObject(svt))
    expect_identical(as.array(svt), expected)
    svt <- svt1 * a2
    expect_true(is(svt, "SVT_SparseArray"))
    expect_true(validObject(svt))
    expect_identical(as.array(svt), expected)
    svt <- a1 * svt2
    expect_true(is(svt, "SVT_SparseArray"))
    expect_true(validObject(svt))
    expect_identical(as.array(svt), expected)
}

test_that("Arith ops between an SVT_SparseArray object and a single value", {
    a1 <- array(0L, 6:4)
    dimnames(a1) <- list(letters[1:6], NULL, LETTERS[1:4])
    a1[c(2:3, 6), 2, 1] <- 101:103
    a1[c(1, 6), 1 , 2] <- 201:202
    a1[1:5, 5, 3] <- -(301:305)
    a1[6, 5, 4] <- NA
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

    a0 <- a1[ , 0, ]
    svt0 <- as(a0, "SVT_SparseArray")
    .test_Arith_SVT1_v2(a0, svt0, 5L)

    m1 <- matrix(1:20, nrow=5)
    m1[2, ] <- 0L
    m1[3:5, 1] <- c(18L, -18L, -15L)
    m1[3, 2:3] <- 0L
    svt1 <- as(m1, "SVT_SparseArray")

    .test_Arith_SVT1_v2(m1, svt1, 0L)
    .test_Arith_SVT1_v2(m1, svt1, 5L)
    .test_Arith_SVT1_v2(m1, svt1, -5L)
    .test_Arith_SVT1_v2(m1, svt1, 1L)
    .test_Arith_SVT1_v2(m1, svt1, -1L)

    m1[1, ] <- m1[1, ] + 0.5
    m1[2, ] <- c(NA, NaN, Inf, -Inf)
    svt1 <- as(m1, "SVT_SparseArray")

    .test_Arith_SVT1_v2(m1, svt1, 0L)
    .test_Arith_SVT1_v2(m1, svt1, 5L)
    .test_Arith_SVT1_v2(m1, svt1, -5L)
    .test_Arith_SVT1_v2(m1, svt1, 1L)
    .test_Arith_SVT1_v2(m1, svt1, -1L)

    ## Not expected to work.
    expect_error(5 / svt1, "not supported")
    expect_error(5 ^ svt1, "not supported")
    expect_error(svt1 + "A", "not supported")
    expect_error(svt1 + 1, "not supported")
    expect_error(svt1 - 1, "not supported")
    expect_error(svt1 * 1:2, "not supported")
    expect_error(svt1 * NA_integer_, "not supported")
    expect_error(svt1 * NA_real_, "not supported")
    expect_error(svt1 * NaN, "not supported")
    expect_error(svt1 * Inf, "not supported")
    expect_error(svt1 * -Inf, "not supported")
    expect_error(svt1 ^ 0, "not supported")
    expect_error(svt1 ^ -2, "not supported")
    expect_error(svt1 / 0, "not supported")
    expect_error(svt1 %% 0, "not supported")
    expect_error(svt1 %/% 0, "not supported")
})

test_that("Arith ops between 2 SVT_SparseArray objects", {
    a1 <- a2 <- array(0L, 6:4)
    dimnames(a1) <- list(letters[1:6], NULL, LETTERS[1:4])
    dimnames(a2) <- list(NULL, letters[22:26], LETTERS[23:26])
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

    a0 <- a2[ , 0, ]
    svt0 <- as(a0, "SVT_SparseArray")
    .test_Arith_SVT1_SVT2(a0, a0, svt0, svt0)
    .test_Arith_SVT1_SVT2(a0, unname(a0), svt0, unname(svt0))
    .test_Arith_SVT1_SVT2(unname(a0), a0, unname(svt0), svt0)

    ## Not expected to work.
    expect_error(svt1 + svt2[ , , -1], "non-conformable")
    expect_error(svt1 - svt2[ , , -1], "non-conformable")
    expect_error(svt1 * svt2[ , , -1], "non-conformable")
    expect_error(svt1 / svt2, "not supported")
    expect_error(svt1 ^ svt2, "not supported")
    expect_error(svt1 %% svt2, "not supported")
    expect_error(svt1 %/% svt2, "not supported")
})

