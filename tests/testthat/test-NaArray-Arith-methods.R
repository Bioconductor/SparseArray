
.test_Arith_NaSVT1_v2 <- function(a1, naa1, v2, strict=TRUE,
                                  relax.MOD.and.IDIV=FALSE)
{
    a <- a1 + v2
    naa <- naa1 + v2
    check_NaArray_object(naa, a, strict=strict)

    a <- a1 - v2
    naa <- naa1 - v2
    check_NaArray_object(naa, a, strict=strict)

    a <- a1 * v2
    naa <- naa1 * v2
    check_NaArray_object(naa, a, strict=strict)

    a <- a1 / v2
    naa <- naa1 / v2
    check_NaArray_object(naa, a, strict=strict)

    if (!(v2 %in% c(0, NaN))) {
        a <- a1 ^ v2
        naa <- naa1 ^ v2
        check_NaArray_object(naa, a, strict=FALSE)
    }

    if (relax.MOD.and.IDIV) {
        reconstructed <- as.array((naa1 %% v2) + v2 * (naa1 %/% v2))
        expect_equal(reconstructed, a1)
    } else {
        if (is.na(v2) || v2 != 0) {
            a <- a1 %% v2
            naa <- naa1 %% v2
            check_NaArray_object(naa, a, strict=FALSE)
        }
        a <- a1 %/% v2
        naa <- naa1 %/% v2
        check_NaArray_object(naa, a, strict=strict)
    }
}

.test_Arith_v1_NaSVT2 <- function(v1, a2, naa2, strict=TRUE,
                                  relax.MOD.and.IDIV=FALSE)
{
    a <- v1 + a2
    naa <- v1 + naa2
    check_NaArray_object(naa, a, strict=strict)

    a <- v1 - a2
    naa <- v1 - naa2
    check_NaArray_object(naa, a, strict=strict)

    a <- v1 * a2
    naa <- v1 * naa2
    check_NaArray_object(naa, a, strict=strict)

    a <- v1 / a2
    naa <- v1 / naa2
    check_NaArray_object(naa, a, strict=strict)

    if (is.na(v1) || v1 != 1) {
        a <- v1 ^ a2
        naa <- v1 ^ naa2
        check_NaArray_object(naa, a, strict=FALSE)
    }

    if (relax.MOD.and.IDIV) {
        reconstructed <- as.array((v1 %% naa2) + naa2 * (v1 %/% naa2))
        expect_true(all(reconstructed == v1 | !is.finite(reconstructed)))
    } else {
        a <- v1 %% a2
        naa <- v1 %% naa2
        check_NaArray_object(naa, a, strict=FALSE)

        a <- v1 %/% a2
        naa <- v1 %/% naa2
        check_NaArray_object(naa, a, strict=strict)
    }
}

### We also test Arith ops between:
### - an NaArray and an ordinary array;
### - an NaArray and SVT_SparseArray object.
.test_Arith_NaSVT1_NaSVT2 <- function(a1, a2, naa1, naa2)
{
    svt1 <- as(a1, "SVT_SparseArray")
    svt2 <- as(a2, "SVT_SparseArray")

    a <- a1 + a2
    naa <- naa1 + naa2
    check_NaArray_object(naa, a)
    expect_identical(naa, naa1 + a2)
    expect_identical(naa, a1 + naa2)
    expect_identical(naa, naa1 + svt2)
    expect_identical(naa, svt1 + naa2)

    a <- a1 - a2
    naa <- naa1 - naa2
    check_NaArray_object(naa, a)
    expect_identical(naa, naa1 - a2)
    expect_identical(naa, a1 - naa2)
    expect_identical(naa, naa1 - svt2)
    expect_identical(naa, svt1 - naa2)

    a <- a1 * a2
    naa <- naa1 * naa2
    check_NaArray_object(naa, a)
    expect_identical(naa, naa1 * a2)
    expect_identical(naa, a1 * naa2)
    expect_identical(naa, naa1 * svt2)
    expect_identical(naa, svt1 * naa2)

    a <- a1 / a2
    naa <- naa1 / naa2
    check_NaArray_object(naa, a)
    expect_identical(naa, naa1 / a2)
    expect_identical(naa, a1 / naa2)
    expect_identical(naa, naa1 / svt2)
    expect_identical(naa, svt1 / naa2)

    a <- a1 ^ a2
    naa <- naa1 ^ naa2
    check_NaArray_object(naa, a)
    expect_identical(naa, naa1 ^ a2)
    expect_identical(naa, a1 ^ naa2)
    expect_error(naa1 ^ svt2, "not supported")
    expect_identical(naa, svt1 ^ naa2)

    a <- a1 %% a2
    naa <- naa1 %% naa2
    check_NaArray_object(naa, a)
    expect_identical(naa, naa1 %% a2)
    expect_identical(naa, a1 %% naa2)
    if (type(naa1) == "double" || type(svt2) == "double") {
        expect_error(naa1 %% svt2, "not supported")
    } else {
        expect_identical(naa, naa1 %% svt2)
    }
    expect_identical(naa, svt1 %% naa2)

    a <- a1 %/% a2
    naa <- naa1 %/% naa2
    check_NaArray_object(naa, a)
    expect_identical(naa, naa1 %/% a2)
    expect_identical(naa, a1 %/% naa2)
    expect_identical(naa, naa1 %/% svt2)
    expect_identical(naa, svt1 %/% naa2)
}

test_that("'Arith' ops between NaArray object and single value", {

    ## --- 3D ---

    a1 <- make_3D_integer_array(NA_integer_)
    a1[2:5, 4, 2] <- a1[ , c(1:2, 4), 3] <- 0L
    a1[c(4, 6), 1, 3] <- 1L
    naa1 <- as(a1, "NaArray")

    .test_Arith_NaSVT1_v2(a1, naa1, 0L)
    .test_Arith_v1_NaSVT2(0L, a1, naa1)
    .test_Arith_NaSVT1_v2(a1, naa1, 5L)
    .test_Arith_v1_NaSVT2(5L, a1, naa1)
    .test_Arith_NaSVT1_v2(a1, naa1, -5L)
    .test_Arith_v1_NaSVT2(-5L, a1, naa1)
    .test_Arith_NaSVT1_v2(a1, naa1, 1L)
    .test_Arith_v1_NaSVT2(1L, a1, naa1)
    .test_Arith_NaSVT1_v2(a1, naa1, -1L)
    .test_Arith_v1_NaSVT2(-1L, a1, naa1)
    .test_Arith_NaSVT1_v2(a1, naa1, 5.1)
    .test_Arith_v1_NaSVT2(5.1, a1, naa1)
    .test_Arith_NaSVT1_v2(a1, naa1, -5.1, relax.MOD.and.IDIV=TRUE)
    .test_Arith_v1_NaSVT2(-5.1, a1, naa1)
    .test_Arith_NaSVT1_v2(a1, naa1, 0.001, relax.MOD.and.IDIV=TRUE)
    .test_Arith_v1_NaSVT2(0.001, a1, naa1)
    .test_Arith_NaSVT1_v2(a1, naa1, -0.001, relax.MOD.and.IDIV=TRUE)
    .test_Arith_v1_NaSVT2(-0.001, a1, naa1)
    .test_Arith_NaSVT1_v2(a1, naa1, 134)
    .test_Arith_v1_NaSVT2(134, a1, naa1)
    .test_Arith_NaSVT1_v2(a1, naa1, Inf)
    .test_Arith_v1_NaSVT2(Inf, a1, naa1)
    .test_Arith_NaSVT1_v2(a1, naa1, -Inf)
    .test_Arith_v1_NaSVT2(-Inf, a1, naa1)
    .test_Arith_NaSVT1_v2(a1, naa1, NA_integer_)
    .test_Arith_v1_NaSVT2(NA_integer_, a1, naa1)
    .test_Arith_NaSVT1_v2(a1, naa1, NA_real_)
    .test_Arith_v1_NaSVT2(NA_real_, a1, naa1)
    .test_Arith_NaSVT1_v2(a1, naa1, NaN, strict=FALSE)
    .test_Arith_v1_NaSVT2(NaN, a1, naa1, strict=FALSE)
    expect_warning(naa1 * 10650000L, "integer overflow")
    expect_warning(10650000L * naa1, "integer overflow")

    a0 <- a1[ , 0, ]
    naa0 <- as(a0, "NaArray")
    .test_Arith_NaSVT1_v2(a0, naa0, 5L)
    .test_Arith_v1_NaSVT2(5L, a0, naa0)

    ## --- 2D ---

    m1 <- a1[ , , 2]
    naa1 <- as(m1, "NaMatrix")

    .test_Arith_NaSVT1_v2(m1, naa1, 0L)
    .test_Arith_v1_NaSVT2(0L, m1, naa1)
    .test_Arith_NaSVT1_v2(m1, naa1, 5L)
    .test_Arith_v1_NaSVT2(5L, m1, naa1)
    .test_Arith_NaSVT1_v2(m1, naa1, -5L)
    .test_Arith_v1_NaSVT2(-5L, m1, naa1)
    .test_Arith_NaSVT1_v2(m1, naa1, 1L)
    .test_Arith_v1_NaSVT2(1L, m1, naa1)
    .test_Arith_NaSVT1_v2(m1, naa1, -1L)
    .test_Arith_v1_NaSVT2(-1L, m1, naa1)
    .test_Arith_NaSVT1_v2(m1, naa1, Inf)
    .test_Arith_v1_NaSVT2(Inf, m1, naa1)
    .test_Arith_NaSVT1_v2(m1, naa1, -Inf)
    .test_Arith_v1_NaSVT2(-Inf, m1, naa1)
    .test_Arith_NaSVT1_v2(m1, naa1, NA_integer_)
    .test_Arith_v1_NaSVT2(NA_integer_, m1, naa1)
    .test_Arith_NaSVT1_v2(m1, naa1, NA_real_)
    .test_Arith_v1_NaSVT2(NA_real_, m1, naa1)
    .test_Arith_NaSVT1_v2(m1, naa1, NaN, strict=FALSE)
    .test_Arith_v1_NaSVT2(NaN, m1, naa1, strict=FALSE)

    m2 <- matrix(1:20, nrow=5)
    m2[2, ] <- m2[-5, 4] <- NA_integer_
    m2[3:5, 1] <- c(18L, -18L, -15L)
    m2[3, 2:3] <- 0L
    naa2 <- as(m2, "NaMatrix")

    .test_Arith_NaSVT1_v2(m2, naa2, 0L)
    .test_Arith_v1_NaSVT2(0L, m2, naa2)
    .test_Arith_NaSVT1_v2(m2, naa2, 5L)
    .test_Arith_v1_NaSVT2(5L, m2, naa2)
    .test_Arith_NaSVT1_v2(m2, naa2, -5L)
    .test_Arith_v1_NaSVT2(-5L, m2, naa2)
    .test_Arith_NaSVT1_v2(m2, naa2, 1L)
    .test_Arith_v1_NaSVT2(1L, m2, naa2)
    .test_Arith_NaSVT1_v2(m2, naa2, -1L)
    .test_Arith_v1_NaSVT2(-1L, m2, naa2)
    .test_Arith_NaSVT1_v2(m2, naa2, Inf)
    .test_Arith_v1_NaSVT2(Inf, m2, naa2)
    .test_Arith_NaSVT1_v2(m2, naa2, -Inf)
    .test_Arith_v1_NaSVT2(-Inf, m2, naa2)
    .test_Arith_NaSVT1_v2(m2, naa2, NA_integer_)
    .test_Arith_v1_NaSVT2(NA_integer_, m2, naa2)
    .test_Arith_NaSVT1_v2(m2, naa2, NA_real_)
    .test_Arith_v1_NaSVT2(NA_real_, m2, naa2)
    .test_Arith_NaSVT1_v2(m2, naa2, NaN, strict=FALSE)
    .test_Arith_v1_NaSVT2(NaN, m2, naa2, strict=FALSE)

    m2[1, ] <- m2[1, ] + 0.5
    m2[2, ] <- c(NA, NaN, Inf, -Inf)
    m2[3, ] <- c(-1, -0.1, 0.1, 1)
    m2[5, 2:3] <- 0
    naa2 <- as(m2, "NaMatrix")

    .test_Arith_NaSVT1_v2(m2, naa2, 0L)
    .test_Arith_v1_NaSVT2(0L, m2, naa2)
    .test_Arith_NaSVT1_v2(m2, naa2, 5L)
    .test_Arith_v1_NaSVT2(5L, m2, naa2, relax.MOD.and.IDIV=TRUE)
    .test_Arith_NaSVT1_v2(m2, naa2, -5L)
    .test_Arith_v1_NaSVT2(-5L, m2, naa2, relax.MOD.and.IDIV=TRUE)
    .test_Arith_NaSVT1_v2(m2, naa2, 1L)
    .test_Arith_v1_NaSVT2(1L, m2, naa2, relax.MOD.and.IDIV=TRUE)
    .test_Arith_NaSVT1_v2(m2, naa2, -1L)
    .test_Arith_v1_NaSVT2(-1L, m2, naa2, relax.MOD.and.IDIV=TRUE)
    .test_Arith_NaSVT1_v2(m2, naa2, Inf)
    .test_Arith_v1_NaSVT2(Inf, m2, naa2)
    .test_Arith_NaSVT1_v2(m2, naa2, -Inf)
    .test_Arith_v1_NaSVT2(-Inf, m2, naa2)
    .test_Arith_NaSVT1_v2(m2, naa2, NA_integer_, strict=FALSE)
    .test_Arith_v1_NaSVT2(NA_integer_, m2, naa2, strict=FALSE)
    .test_Arith_NaSVT1_v2(m2, naa2, NA_real_, strict=FALSE)
    .test_Arith_v1_NaSVT2(NA_real_, m2, naa2, strict=FALSE)
    .test_Arith_NaSVT1_v2(m2, naa2, NaN, strict=FALSE)
    .test_Arith_v1_NaSVT2(NaN, m2, naa2, strict=FALSE)

    ## --- Not expected to work ---

    expect_error(naa2 + "A", "not supported")
    expect_error(naa2 * 1:2, "not supported")
    expect_error(naa2 ^ 0, "not supported")
    expect_error(naa2 ^ NaN, "not supported")
    expect_error(naa2 %% 0, "not supported")
    expect_error(1 ^ naa2, "not supported")
})

test_that("'Arith' ops between 2 NaArray objects", {
    a1 <- a2 <- array(NA_integer_, 6:4)
    dimnames(a1) <- list(letters[1:6], NULL, LETTERS[1:4])
    dimnames(a2) <- list(NULL, letters[22:26], LETTERS[23:26])

    naa1 <- as(a1, "NaArray")
    naa2 <- as(a2, "NaArray")

    .test_Arith_NaSVT1_NaSVT2(a1, a2, naa1, naa2)
    .test_Arith_NaSVT1_NaSVT2(a2, a1, naa2, naa1)

    a1[c(2:3, 6), 2, 1] <- 101:103
    a2[2:4, 2, 1] <- 1001:1003
    a1[c(1, 6), 1 , 2] <- 201:202
    a2[c(3, 5), 2 , 2] <- 2001:2002
    a1[1:5, 5, 3] <- 301:305
    a2[c(2:4, 6), 5, 3] <- 3001:3004
    a2[6, 5, 4] <- 0L
    naa1 <- as(a1, "NaArray")
    naa2 <- as(a2, "NaArray")

    .test_Arith_NaSVT1_NaSVT2(a1, a2, naa1, naa2)
    .test_Arith_NaSVT1_NaSVT2(a2, a1, naa2, naa1)
    .test_Arith_NaSVT1_NaSVT2(a1, a1, naa1, naa1)
    .test_Arith_NaSVT1_NaSVT2(a2, a2, naa2, naa2)
    .test_Arith_NaSVT1_NaSVT2(a1, a1 + a2, naa1, naa1 + naa2)
    .test_Arith_NaSVT1_NaSVT2(a1, a2 - a1, naa1, naa2 - naa1)
    .test_Arith_NaSVT1_NaSVT2(a2, a1 - a2, naa2, naa1 - naa2)
    .test_Arith_NaSVT1_NaSVT2(a1, a1 * a2, naa1, naa1 * naa2)

    dimnames(a1) <- dimnames(naa1) <- NULL
    .test_Arith_NaSVT1_NaSVT2(a1, a2, naa1, naa2)
    .test_Arith_NaSVT1_NaSVT2(a2, a1, naa2, naa1)

    a1[ , , ] <- NA_integer_
    naa1 <- as(a1, "NaArray")

    .test_Arith_NaSVT1_NaSVT2(a1, a2, naa1, naa2)
    .test_Arith_NaSVT1_NaSVT2(a2, a1, naa2, naa1)
    .test_Arith_NaSVT1_NaSVT2(a1, a1, naa1, naa1)
    .test_Arith_NaSVT1_NaSVT2(a2, a2, naa2, naa2)

    a1[ , , ] <- NA_real_
    naa1 <- as(a1, "NaArray")

    .test_Arith_NaSVT1_NaSVT2(a1, a2, naa1, naa2)
    .test_Arith_NaSVT1_NaSVT2(a2, a1, naa2, naa1)
    .test_Arith_NaSVT1_NaSVT2(a1, a1, naa1, naa1)
    .test_Arith_NaSVT1_NaSVT2(a2, a2, naa2, naa2)

    a0 <- a2[ , 0, ]
    naa0 <- as(a0, "NaArray")
    .test_Arith_NaSVT1_NaSVT2(a0, a0, naa0, naa0)
    .test_Arith_NaSVT1_NaSVT2(a0, unname(a0), naa0, unname(naa0))
    .test_Arith_NaSVT1_NaSVT2(unname(a0), a0, unname(naa0), naa0)

    ## --- Not expected to work ---

    expect_error(naa1 + naa2[ , , -1], "non-conformable")
    expect_error(naa1 - naa2[ , , -1], "non-conformable")
    expect_error(naa1 * naa2[ , , -1], "non-conformable")
    expect_error(naa1 / naa2[ , , -1], "non-conformable")
    expect_error(naa1 ^ naa2[ , , -1], "non-conformable")
    expect_error(naa1 %% naa2[ , , -1], "non-conformable")
    expect_error(naa1 %/% naa2[ , , -1], "non-conformable")
})

test_that("unary minus on an NaArray object", {
    a1 <- make_3D_integer_array(NA_integer_)
    naa1 <- as(a1, "NaArray")
    check_NaArray_object(- naa1, - a1)
    expect_identical(- (- naa1), naa1)

    a2 <- make_3D_double_array(NA_real_)
    naa2 <- as(a2, "NaArray")
    check_NaArray_object(- naa2, - a2)
    expect_identical(- (- naa2), naa2)

    a3 <- make_3D_complex_array(NA_complex_)
    naa3 <- as(a3, "NaArray")
    check_NaArray_object(- naa3, - a3)
    expect_identical(- (- naa3), naa3)
})

