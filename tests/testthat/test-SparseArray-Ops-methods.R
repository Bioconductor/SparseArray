### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Arith' ops

.test_Arith_SVT1_v2 <- function(a1, svt1, v2, relax.MOD.and.IDIV=FALSE)
{
    if (is.infinite(v2)) {
        expect_error(svt1 * v2, "not supported")
        expect_error(v2 * svt1, "not supported")
    } else {
        expected <- a1 * v2
        svt <- svt1 * v2
        expect_true(is(svt, "SVT_SparseArray"))
        expect_true(validObject(svt))
        expect_identical(as.array(svt), expected)
        svt <- v2 * svt1
        expect_true(is(svt, "SVT_SparseArray"))
        expect_true(validObject(svt))
        expect_identical(as.array(svt), expected)
    }

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
        expect_equal(as.array(svt), expected)
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
    check_SparseArray_object(- svt1, "SVT_SparseArray", - a1)
    expect_identical(- svt1, as(- a1, "SVT_SparseArray"))
    expect_identical(- (- svt1), svt1)

    a2 <- make_3D_double_array()
    svt2 <- as(a2, "SVT_SparseArray")
    check_SparseArray_object(- svt2, "SVT_SparseArray", - a2)
    expect_identical(- svt2, as(- a2, "SVT_SparseArray"))
    expect_identical(- (- svt2), svt2)

    a3 <- make_3D_complex_array()
    svt3 <- as(a3, "SVT_SparseArray")
    check_SparseArray_object(- svt3, "SVT_SparseArray", - a3)
    expect_identical(- svt3, as(- a3, "SVT_SparseArray"))
    expect_identical(- (- svt3), svt3)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 'Compare' ops

.test_Compare_SVT_val <- function(a, svt, y)
{
    ## svt == y
    if (y == 0) {
	expect_error(svt == y, "not supported")
	expect_error(y == svt, "not supported")
    } else {
        expected <- a == y
        current <- svt == y
        expect_true(is(current, "SVT_SparseArray"))
        expect_true(validObject(current))
        expect_identical(type(current), "logical")
        expect_identical(as.array(current), expected)
        current <- y == svt
        expect_true(is(current, "SVT_SparseArray"))
        expect_true(validObject(current))
        expect_identical(type(current), "logical")
        expect_identical(as.array(current), expected)
    }

    ## svt != y
    if (y != 0) {
	expect_error(svt != y, "not supported")
	expect_error(y != svt, "not supported")
    } else {
        expected <- a != y
        current <- svt != y
        expect_true(is(current, "SVT_SparseArray"))
        expect_true(validObject(current))
        expect_identical(type(current), "logical")
        expect_identical(as.array(current), expected)
        current <- y != svt
        expect_true(is(current, "SVT_SparseArray"))
        expect_true(validObject(current))
        expect_identical(type(current), "logical")
        expect_identical(as.array(current), expected)
    }

    ## svt <= y
    if (type(svt) == "complex" || type(y) == "complex") {
	expect_error(svt <= y, "invalid comparison with complex values")
	expect_error(y >= svt, "invalid comparison with complex values")
    } else if (y >= 0) {
	expect_error(svt <= y, "not supported")
	expect_error(y >= svt, "not supported")
    } else {
        expected <- a <= y
        current <- svt <= y
        expect_true(is(current, "SVT_SparseArray"))
        expect_true(validObject(current))
        expect_identical(type(current), "logical")
        expect_identical(as.array(current), expected)
        current <- y >= svt
        expect_true(is(current, "SVT_SparseArray"))
        expect_true(validObject(current))
        expect_identical(type(current), "logical")
        expect_identical(as.array(current), expected)
    }

    ## svt >= y
    if (type(svt) == "complex" || type(y) == "complex") {
	expect_error(svt >= y, "invalid comparison with complex values")
	expect_error(y <= svt, "invalid comparison with complex values")
    } else if (y <= 0) {
	expect_error(svt >= y, "not supported")
	expect_error(y <= svt, "not supported")
    } else {
        expected <- a >= y
        current <- svt >= y
        expect_true(is(current, "SVT_SparseArray"))
        expect_true(validObject(current))
        expect_identical(type(current), "logical")
        expect_identical(as.array(current), expected)
        current <- y <= svt
        expect_true(is(current, "SVT_SparseArray"))
        expect_true(validObject(current))
        expect_identical(type(current), "logical")
        expect_identical(as.array(current), expected)
    }

    ## svt < y
    if (type(svt) == "complex" || type(y) == "complex") {
        expect_error(svt < y, "invalid comparison with complex values")
        expect_error(y > svt, "invalid comparison with complex values")
    } else if (type(y) %in% c("raw", "logical") || y > 0) {
	expect_error(svt < y, "not supported")
        expect_error(y > svt, "not supported")
    } else {
        expected <- a < y
        current <- svt < y
        expect_true(is(current, "SVT_SparseArray"))
        expect_true(validObject(current))
        expect_identical(type(current), "logical")
        expect_identical(as.array(current), expected)
        current <- y > svt
        expect_true(is(current, "SVT_SparseArray"))
        expect_true(validObject(current))
        expect_identical(type(current), "logical")
        expect_identical(as.array(current), expected)
    }

    ## svt > y
    if (type(svt) == "complex" || type(y) == "complex") {
        expect_error(svt > y, "invalid comparison with complex values")
        expect_error(y < svt, "invalid comparison with complex values")
    } else if (y < 0) {
        expect_error(svt > y, "not supported")
        expect_error(y < svt, "not supported")
    } else {
        expected <- a > y
        current <- svt > y
        expect_true(is(current, "SVT_SparseArray"))
        expect_true(validObject(current))
        expect_identical(type(current), "logical")
        expect_identical(as.array(current), expected)
        current <- y < svt
        expect_true(is(current, "SVT_SparseArray"))
        expect_true(validObject(current))
        expect_identical(type(current), "logical")
        expect_identical(as.array(current), expected)
    }
}

.test_Compare_SVT1_SVT2 <- function(a1, a2, svt1, svt2)
{
    ## svt1 != svt2
    expected <- a1 != a2
    current <- svt1 != svt2
    expect_true(is(current, "SVT_SparseArray"))
    expect_true(validObject(current))
    expect_identical(type(current), "logical")
    expect_identical(as.array(current), expected)

    ## svt1 < svt2
    if (type(svt1) == "complex" || type(svt2) == "complex") {
        expect_error(svt1 < svt2, "invalid comparison with complex values")
    } else {
        expected <- a1 < a2
        current <- svt1 < svt2
        expect_true(is(current, "SVT_SparseArray"))
        expect_true(validObject(current))
        expect_identical(type(current), "logical")
        expect_identical(as.array(current), expected)
    }

    ## svt1 > svt2
    if (type(svt1) == "complex" || type(svt2) == "complex") {
        expect_error(svt1 > svt2, "invalid comparison with complex values")
    } else {
        expected <- a1 > a2
        current <- svt1 > svt2
        expect_true(is(current, "SVT_SparseArray"))
        expect_true(validObject(current))
        expect_identical(type(current), "logical")
        expect_identical(as.array(current), expected)
    }
}

test_that("'Compare' ops between SVT_SparseArray object and single value", {
    ## "integer" array.
    a1 <- array(0L, 6:4)
    dimnames(a1) <- list(letters[1:6], NULL, LETTERS[1:4])
    a1[c(2:3, 6), 2, 1] <- 101:103
    a1[c(1, 6), 1 , 2] <- 201:202
    a1[1:5, 5, 3] <- -(301:305)
    a1[1, 3, 2] <- NA

    ## "double" array.
    a2 <- a1
    a2[3:6, 3, 2] <- c(NA, NaN, Inf, -Inf)

    ## "raw" array.
    a3 <- a1
    suppressWarnings(type(a3) <- "raw")

    ## "logical" array.
    a4 <- a1
    type(a4) <- "logical"

    ## "complex" array.
    a5 <- a2
    a5[4, , 1] <- (0.2 - 3i)^(1:5)
    a5[4:6, 4, 2] <- c(NaN, Inf, -Inf) + 3i
    a5[4:6, 4, 2] <- 0.2 + 1i * c(NaN, Inf, -Inf)

    ## Empty array.
    a0 <- a1[ , 0, ]

    for (a in list(a1, a2, a3, a4, a5, a0)) {
        svt <- as(a, "SVT_SparseArray")
        .test_Compare_SVT_val(a, svt, integer(1))
        .test_Compare_SVT_val(a, svt, 102L)
        .test_Compare_SVT_val(a, svt, -302L)
        .test_Compare_SVT_val(a, svt, raw(1))
        .test_Compare_SVT_val(a, svt, as.raw(102L))
        .test_Compare_SVT_val(a, svt, FALSE)
        .test_Compare_SVT_val(a, svt, TRUE)
        .test_Compare_SVT_val(a, svt, double(1))
        .test_Compare_SVT_val(a, svt, 102.5)
        .test_Compare_SVT_val(a, svt, -302.0)
        .test_Compare_SVT_val(a, svt, Inf)
        .test_Compare_SVT_val(a, svt, -Inf)
        .test_Compare_SVT_val(a, svt, complex(1))
        .test_Compare_SVT_val(a, svt, 0+9.5i)
        .test_Compare_SVT_val(a, svt, -302+0i)
        .test_Compare_SVT_val(a, svt, -302+9.5i)
        .test_Compare_SVT_val(a, svt, Inf+9.5i)
        .test_Compare_SVT_val(a, svt, -Inf+9.5i)
    }

    ## Not expected to work.
    svt1 <- as(a1, "SVT_SparseArray")
    for (y in list(NA, NaN, 11:15, numeric(0), list(-0.22))) {
        expect_error(svt1 == y, "not supported")
        expect_error(y == svt1, "not supported")
        expect_error(svt1 != y, "not supported")
        expect_error(y != svt1, "not supported")
        expect_error(svt1 <= y, "not supported")
        expect_error(y >= svt1, "not supported")
        expect_error(svt1 >= y, "not supported")
        expect_error(y <= svt1, "not supported")
        expect_error(svt1 < y, "not supported")
        expect_error(y > svt1, "not supported")
        expect_error(svt1 > y, "not supported")
        expect_error(y < svt1, "not supported")
    }
})

test_that("'Compare' ops between 2 SVT_SparseArray objects", {
    ## "integer" arrays.
    a1 <- array(0L, 6:4)
    dimnames(a1) <- list(letters[1:6], NULL, LETTERS[1:4])
    a1[c(2:3, 6), 2, 1] <- 101:103
    a1[c(1, 6), 1 , 2] <- 201:202
    a1[1:5, 5, 3] <- -(301:305)
    a1[ , 3, 2] <- NA

    a2 <- a1
    dimnames(a2) <- list(NULL, letters[22:26], LETTERS[23:26])
    a2[c(1, 2, 6), 2, 1] <- c(-3L, 100L, 104L)
    a2[2:6, 3, 1] <- -2:2
    a2[ , 3, 2] <- c(-1:3, NA)

    svt1 <- as(a1, "SVT_SparseArray")
    svt2 <- as(a2, "SVT_SparseArray")
    .test_Compare_SVT1_SVT2(a1, a2, svt1, svt2)
    .test_Compare_SVT1_SVT2(a2, a1, svt2, svt1)
    .test_Compare_SVT1_SVT2(a1, a1, svt1, svt1)
    .test_Compare_SVT1_SVT2(a2, a2, svt2, svt2)

    ## "double" arrays.
    a3 <- a1
    a3[ , 1, 3] <- NaN
    a3[ , 2, 3] <- Inf
    a3[ , 3, 3] <- -Inf
    a4 <- a2
    a4[ , 1:3, 3] <- c(-1, NA, 2, NaN, Inf, -Inf)

    svt3 <- as(a3, "SVT_SparseArray")
    svt4 <- as(a4, "SVT_SparseArray")
    .test_Compare_SVT1_SVT2(a3, a4, svt3, svt4)
    .test_Compare_SVT1_SVT2(a4, a3, svt4, svt3)
    .test_Compare_SVT1_SVT2(a3, a3, svt3, svt3)
    .test_Compare_SVT1_SVT2(a4, a4, svt4, svt4)

    ## "raw" arrays.
    a5 <- a1
    a6 <- a2
    suppressWarnings(type(a5) <- type(a6) <- "raw")

    svt5 <- as(a5, "SVT_SparseArray")
    svt6 <- as(a6, "SVT_SparseArray")
    .test_Compare_SVT1_SVT2(a5, a6, svt5, svt6)
    .test_Compare_SVT1_SVT2(a6, a5, svt6, svt5)
    .test_Compare_SVT1_SVT2(a5, a5, svt5, svt5)
    .test_Compare_SVT1_SVT2(a6, a6, svt6, svt6)

    ## "logical" arrays.
    a7 <- a1
    a8 <- a2
    type(a7) <- type(a8) <- "logical"

    svt7 <- as(a7, "SVT_SparseArray")
    svt8 <- as(a8, "SVT_SparseArray")
    .test_Compare_SVT1_SVT2(a7, a8, svt7, svt8)
    .test_Compare_SVT1_SVT2(a8, a7, svt8, svt7)
    .test_Compare_SVT1_SVT2(a7, a7, svt7, svt7)
    .test_Compare_SVT1_SVT2(a8, a8, svt8, svt8)

    ## "complex" arrays.
    a9 <- a3
    a10 <- a4
    a9[4, , 1] <- (0.2 - 3i)^(1:5)
    a9[4:6, 4, 2] <- 0.2 + 1i * c(NaN, Inf, -Inf)
    a10[4:6, 3, 2] <- c(NaN, Inf, -Inf) + 3i
    a10[ , 4, 2] <- NaN

    svt9 <- as(a9, "SVT_SparseArray")
    svt10 <- as(a10, "SVT_SparseArray")
    .test_Compare_SVT1_SVT2(a9, a10, svt9, svt10)
    .test_Compare_SVT1_SVT2(a10, a9, svt10, svt9)
    .test_Compare_SVT1_SVT2(a9, a9, svt9, svt9)
    .test_Compare_SVT1_SVT2(a10, a10, svt10, svt10)

    ## "integer" <op> "double".
    .test_Compare_SVT1_SVT2(a1, a4, svt1, svt4)
    .test_Compare_SVT1_SVT2(a4, a1, svt4, svt1)
    .test_Compare_SVT1_SVT2(a2, a3, svt2, svt3)
    .test_Compare_SVT1_SVT2(a3, a2, svt3, svt2)

    ## "integer" <op> "raw".
    .test_Compare_SVT1_SVT2(a1, a6, svt1, svt6)
    .test_Compare_SVT1_SVT2(a6, a1, svt6, svt1)
    .test_Compare_SVT1_SVT2(a2, a5, svt2, svt5)
    .test_Compare_SVT1_SVT2(a5, a2, svt5, svt2)

    ## "integer" <op> "logical".
    .test_Compare_SVT1_SVT2(a1, a8, svt1, svt8)
    .test_Compare_SVT1_SVT2(a8, a1, svt8, svt1)
    .test_Compare_SVT1_SVT2(a2, a7, svt2, svt7)
    .test_Compare_SVT1_SVT2(a7, a2, svt7, svt2)

    ## "integer" <op> "complex".
    .test_Compare_SVT1_SVT2(a1, a10, svt1, svt10)
    .test_Compare_SVT1_SVT2(a10, a1, svt10, svt1)
    .test_Compare_SVT1_SVT2(a2, a9, svt2, svt9)
    .test_Compare_SVT1_SVT2(a9, a2, svt9, svt2)

    ## "double" <op> "raw".
    .test_Compare_SVT1_SVT2(a3, a6, svt3, svt6)
    .test_Compare_SVT1_SVT2(a6, a3, svt6, svt3)
    .test_Compare_SVT1_SVT2(a4, a5, svt4, svt5)
    .test_Compare_SVT1_SVT2(a5, a4, svt5, svt4)

    ## "double" <op> "logical".
    .test_Compare_SVT1_SVT2(a3, a8, svt3, svt8)
    .test_Compare_SVT1_SVT2(a8, a3, svt8, svt3)
    .test_Compare_SVT1_SVT2(a4, a7, svt4, svt7)
    .test_Compare_SVT1_SVT2(a7, a4, svt7, svt4)

    ## "double" <op> "complex".
    .test_Compare_SVT1_SVT2(a3, a10, svt3, svt10)
    .test_Compare_SVT1_SVT2(a10, a3, svt10, svt3)
    .test_Compare_SVT1_SVT2(a4, a9, svt4, svt9)
    .test_Compare_SVT1_SVT2(a9, a4, svt9, svt4)

    ## "raw" <op> "logical".
    .test_Compare_SVT1_SVT2(a5, a8, svt5, svt8)
    .test_Compare_SVT1_SVT2(a8, a5, svt8, svt5)
    .test_Compare_SVT1_SVT2(a6, a7, svt6, svt7)
    .test_Compare_SVT1_SVT2(a7, a6, svt7, svt6)

    ## "raw" <op> "complex".
    .test_Compare_SVT1_SVT2(a5, a10, svt5, svt10)
    .test_Compare_SVT1_SVT2(a10, a5, svt10, svt5)
    .test_Compare_SVT1_SVT2(a6, a9, svt6, svt9)
    .test_Compare_SVT1_SVT2(a9, a6, svt9, svt6)

    ## "logical" <op> "complex".
    .test_Compare_SVT1_SVT2(a7, a10, svt7, svt10)
    .test_Compare_SVT1_SVT2(a10, a7, svt10, svt7)
    .test_Compare_SVT1_SVT2(a8, a9, svt8, svt9)
    .test_Compare_SVT1_SVT2(a9, a8, svt9, svt8)

    ## Not expected to work.
    expect_error(svt1 != svt2[ , , -1], "non-conformable")
    expect_error(svt1 < svt2[ , , -1], "non-conformable")
    expect_error(svt1 > svt2[ , , -1], "non-conformable")
    expect_error(svt1 == svt2, "not supported")
    expect_error(svt1 <= svt2, "not supported")
    expect_error(svt1 >= svt2, "not supported")
})

