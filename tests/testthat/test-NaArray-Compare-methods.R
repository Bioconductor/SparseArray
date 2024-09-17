
.test_Compare_NaSVT1_v2 <- function(a1, naa1, y)
{
    ## naa1 == y
    a <- a1 == y
    naa <- naa1 == y
    check_NaArray_object(naa, a)
    naa <- y == naa1
    check_NaArray_object(naa, a)

    ## naa1 != y
    a <- a1 != y
    naa <- naa1 != y
    check_NaArray_object(naa, a)
    naa <- y != naa1
    check_NaArray_object(naa, a)

    ## naa1 <= y
    if (type(naa1) == "complex" || type(y) == "complex") {
	expect_error(naa1 <= y, "invalid comparison with complex values")
	expect_error(y >= naa1, "invalid comparison with complex values")
    } else {
        a <- a1 <= y
        naa <- naa1 <= y
        check_NaArray_object(naa, a)
        naa <- y >= naa1
        check_NaArray_object(naa, a)
    }

    ## naa1 >= y
    if (type(naa1) == "complex" || type(y) == "complex") {
	expect_error(naa1 >= y, "invalid comparison with complex values")
	expect_error(y <= naa1, "invalid comparison with complex values")
    } else {
        a <- a1 >= y
        naa <- naa1 >= y
        check_NaArray_object(naa, a)
        naa <- y <= naa1
        check_NaArray_object(naa, a)
    }

    ## naa1 < y
    if (type(naa1) == "complex" || type(y) == "complex") {
        expect_error(naa1 < y, "invalid comparison with complex values")
        expect_error(y > naa1, "invalid comparison with complex values")
    } else {
        a <- a1 < y
        naa <- naa1 < y
        check_NaArray_object(naa, a)
        naa <- y > naa1
        check_NaArray_object(naa, a)
    }

    ## naa1 > y
    if (type(naa1) == "complex" || type(y) == "complex") {
        expect_error(naa1 > y, "invalid comparison with complex values")
        expect_error(y < naa1, "invalid comparison with complex values")
    } else {
        a <- a1 > y
        naa <- naa1 > y
        check_NaArray_object(naa, a)
        naa <- y < naa1
        check_NaArray_object(naa, a)
    }
}

### We also test Compare ops between an NaArray and SVT_SparseArray object!
.test_Compare_NaSVT1_NaSVT2 <- function(a1, a2, naa1, naa2)
{
    svt1 <- as(a1, "SVT_SparseArray")
    svt2 <- as(a2, "SVT_SparseArray")

    ## naa1 == naa2
    a <- a1 == a2
    naa <- naa1 == naa2
    check_NaArray_object(naa, a)
    expect_identical(naa, naa1 == svt2)
    expect_identical(naa, svt1 == naa2)

    ## naa1 != naa2
    a <- a1 != a2
    naa <- naa1 != naa2
    check_NaArray_object(naa, a)
    expect_identical(naa, naa1 != svt2)
    expect_identical(naa, svt1 != naa2)

    ## naa1 <= naa2
    if (type(naa1) == "complex" || type(naa2) == "complex") {
        expect_error(naa1 <= naa2, "invalid comparison with complex values")
    } else {
        a <- a1 <= a2
        naa <- naa1 <= naa2
        check_NaArray_object(naa, a)
        expect_identical(naa, naa1 <= svt2)
        expect_identical(naa, svt1 <= naa2)
    }

    ## naa1 >= naa2
    if (type(naa1) == "complex" || type(naa2) == "complex") {
        expect_error(naa1 >= naa2, "invalid comparison with complex values")
    } else {
        a <- a1 >= a2
        naa <- naa1 >= naa2
        check_NaArray_object(naa, a)
        expect_identical(naa, naa1 >= svt2)
        expect_identical(naa, svt1 >= naa2)
    }

    ## naa1 < naa2
    if (type(naa1) == "complex" || type(naa2) == "complex") {
        expect_error(naa1 < naa2, "invalid comparison with complex values")
    } else {
        a <- a1 < a2
        naa <- naa1 < naa2
        check_NaArray_object(naa, a)
        expect_identical(naa, naa1 < svt2)
        expect_identical(naa, svt1 < naa2)
    }

    ## naa1 > naa2
    if (type(naa1) == "complex" || type(naa2) == "complex") {
        expect_error(naa1 > naa2, "invalid comparison with complex values")
    } else {
        a <- a1 > a2
        naa <- naa1 > naa2
        check_NaArray_object(naa, a)
        expect_identical(naa, naa1 > svt2)
        expect_identical(naa, svt1 > naa2)
    }
}

test_that("'Compare' ops between NaArray object and single value", {
    ## "integer" array.
    a1 <- array(NA_integer_, 6:4)
    dimnames(a1) <- list(letters[1:6], NULL, LETTERS[1:4])
    a1[c(2:3, 6), 2, 1] <- 101:103
    a1[c(1, 6), 1 , 2] <- 201:202
    a1[1:5, 5, 3] <- -(301:305)
    a1[1, 3, 2] <- NA

    ## "double" array.
    a2 <- a1
    a2[3:6, 3, 2] <- c(NA, NaN, Inf, -Inf)

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

    for (a in list(a1, a2, a4, a5, a0)) {
        naa <- as(a, "NaArray")
        .test_Compare_NaSVT1_v2(a, naa, NA)
        .test_Compare_NaSVT1_v2(a, naa, NA_integer_)
        .test_Compare_NaSVT1_v2(a, naa, NA_real_)
        .test_Compare_NaSVT1_v2(a, naa, NaN)
        .test_Compare_NaSVT1_v2(a, naa, NA_complex_)
        .test_Compare_NaSVT1_v2(a, naa, integer(1))
        .test_Compare_NaSVT1_v2(a, naa, 102L)
        .test_Compare_NaSVT1_v2(a, naa, -302L)
        .test_Compare_NaSVT1_v2(a, naa, raw(1))
        .test_Compare_NaSVT1_v2(a, naa, as.raw(102L))
        .test_Compare_NaSVT1_v2(a, naa, FALSE)
        .test_Compare_NaSVT1_v2(a, naa, TRUE)
        .test_Compare_NaSVT1_v2(a, naa, double(1))
        .test_Compare_NaSVT1_v2(a, naa, 102.5)
        .test_Compare_NaSVT1_v2(a, naa, -302.0)
        .test_Compare_NaSVT1_v2(a, naa, Inf)
        .test_Compare_NaSVT1_v2(a, naa, -Inf)
        .test_Compare_NaSVT1_v2(a, naa, complex(1))
        .test_Compare_NaSVT1_v2(a, naa, 0+9.5i)
        .test_Compare_NaSVT1_v2(a, naa, -302+0i)
        .test_Compare_NaSVT1_v2(a, naa, -302+9.5i)
        .test_Compare_NaSVT1_v2(a, naa, Inf+9.5i)
        .test_Compare_NaSVT1_v2(a, naa, -Inf+9.5i)
    }

    ## Not expected to work.
    naa1 <- as(a1, "NaArray")
    for (y in list(11:15, numeric(0), list(-0.22))) {
        expect_error(naa1 == y, "not[\\s]+supported", perl=TRUE)
        expect_error(y == naa1, "not[\\s]+supported", perl=TRUE)
        expect_error(naa1 != y, "not[\\s]+supported", perl=TRUE)
        expect_error(y != naa1, "not[\\s]+supported", perl=TRUE)
        expect_error(naa1 <= y, "not[\\s]+supported", perl=TRUE)
        expect_error(y >= naa1, "not[\\s]+supported", perl=TRUE)
        expect_error(naa1 >= y, "not[\\s]+supported", perl=TRUE)
        expect_error(y <= naa1, "not[\\s]+supported", perl=TRUE)
        expect_error(naa1 < y, "not[\\s]+supported", perl=TRUE)
        expect_error(y > naa1, "not[\\s]+supported", perl=TRUE)
        expect_error(naa1 > y, "not[\\s]+supported", perl=TRUE)
        expect_error(y < naa1, "not[\\s]+supported", perl=TRUE)
    }
})

test_that("'Compare' ops between 2 NaArray objects", {
    ## "integer" arrays.
    a1 <- array(NA_integer_, 6:4)
    dimnames(a1) <- list(letters[1:6], NULL, LETTERS[1:4])
    a1[c(2:3, 6), 2, 1] <- 101:103
    a1[c(1, 6), 1 , 2] <- 201:202
    a1[1:5, 5, 3] <- -(301:305)

    a2 <- a1
    dimnames(a2) <- list(NULL, letters[22:26], LETTERS[23:26])
    a2[c(1, 2, 6), 2, 1] <- c(-3L, 100L, 104L)
    a2[2:6, 3, 1] <- -2:2
    a2[-6, 3, 2] <- -1:3
    a2[5, , 1] <- a2[ , 1:2, 3] <- a2[ , , 4] <- 0L

    naa1 <- as(a1, "NaArray")
    naa2 <- as(a2, "NaArray")
    .test_Compare_NaSVT1_NaSVT2(a1, a2, naa1, naa2)
    .test_Compare_NaSVT1_NaSVT2(a2, a1, naa2, naa1)
    .test_Compare_NaSVT1_NaSVT2(a1, a1, naa1, naa1)
    .test_Compare_NaSVT1_NaSVT2(a2, a2, naa2, naa2)

    ## "double" arrays.
    a3 <- a1
    a3[ , 1, 3] <- NaN
    a3[ , 2, 3] <- Inf
    a3[ , 3, 3] <- -Inf
    a4 <- a2
    a4[ , 1:3, 3] <- c(-1, NA, 2, NaN, Inf, -Inf)

    naa3 <- as(a3, "NaArray")
    naa4 <- as(a4, "NaArray")
    .test_Compare_NaSVT1_NaSVT2(a3, a4, naa3, naa4)
    .test_Compare_NaSVT1_NaSVT2(a4, a3, naa4, naa3)
    .test_Compare_NaSVT1_NaSVT2(a3, a3, naa3, naa3)
    .test_Compare_NaSVT1_NaSVT2(a4, a4, naa4, naa4)

    ## "raw" arrays.
    #a5 <- a1
    #a6 <- a2
    #suppressWarnings(type(a5) <- type(a6) <- "raw")

    #naa5 <- as(a5, "NaArray")
    #naa6 <- as(a6, "NaArray")
    #.test_Compare_NaSVT1_NaSVT2(a5, a6, naa5, naa6)
    #.test_Compare_NaSVT1_NaSVT2(a6, a5, naa6, naa5)
    #.test_Compare_NaSVT1_NaSVT2(a5, a5, naa5, naa5)
    #.test_Compare_NaSVT1_NaSVT2(a6, a6, naa6, naa6)

    ## "logical" arrays.
    a7 <- a1
    a8 <- a2
    type(a7) <- type(a8) <- "logical"

    naa7 <- as(a7, "NaArray")
    naa8 <- as(a8, "NaArray")
    .test_Compare_NaSVT1_NaSVT2(a7, a8, naa7, naa8)
    .test_Compare_NaSVT1_NaSVT2(a8, a7, naa8, naa7)
    .test_Compare_NaSVT1_NaSVT2(a7, a7, naa7, naa7)
    .test_Compare_NaSVT1_NaSVT2(a8, a8, naa8, naa8)

    ## "complex" arrays.
    a9 <- a3
    a10 <- a4
    a9[4, , 1] <- (0.2 - 3i)^(1:5)
    a9[4:6, 4, 2] <- 0.2 + 1i * c(NaN, Inf, -Inf)
    a10[4:6, 3, 2] <- c(NaN, Inf, -Inf) + 3i
    a10[ , 4, 2] <- NaN

    naa9 <- as(a9, "NaArray")
    naa10 <- as(a10, "NaArray")
    .test_Compare_NaSVT1_NaSVT2(a9, a10, naa9, naa10)
    .test_Compare_NaSVT1_NaSVT2(a10, a9, naa10, naa9)
    .test_Compare_NaSVT1_NaSVT2(a9, a9, naa9, naa9)
    .test_Compare_NaSVT1_NaSVT2(a10, a10, naa10, naa10)

    ## "integer" <op> "double".
    .test_Compare_NaSVT1_NaSVT2(a1, a4, naa1, naa4)
    .test_Compare_NaSVT1_NaSVT2(a4, a1, naa4, naa1)
    .test_Compare_NaSVT1_NaSVT2(a2, a3, naa2, naa3)
    .test_Compare_NaSVT1_NaSVT2(a3, a2, naa3, naa2)

    ## "integer" <op> "raw".
    #.test_Compare_NaSVT1_NaSVT2(a1, a6, naa1, naa6)
    #.test_Compare_NaSVT1_NaSVT2(a6, a1, naa6, naa1)
    #.test_Compare_NaSVT1_NaSVT2(a2, a5, naa2, naa5)
    #.test_Compare_NaSVT1_NaSVT2(a5, a2, naa5, naa2)

    ## "integer" <op> "logical".
    .test_Compare_NaSVT1_NaSVT2(a1, a8, naa1, naa8)
    .test_Compare_NaSVT1_NaSVT2(a8, a1, naa8, naa1)
    .test_Compare_NaSVT1_NaSVT2(a2, a7, naa2, naa7)
    .test_Compare_NaSVT1_NaSVT2(a7, a2, naa7, naa2)

    ## "integer" <op> "complex".
    .test_Compare_NaSVT1_NaSVT2(a1, a10, naa1, naa10)
    .test_Compare_NaSVT1_NaSVT2(a10, a1, naa10, naa1)
    .test_Compare_NaSVT1_NaSVT2(a2, a9, naa2, naa9)
    .test_Compare_NaSVT1_NaSVT2(a9, a2, naa9, naa2)

    ## "double" <op> "raw".
    #.test_Compare_NaSVT1_NaSVT2(a3, a6, naa3, naa6)
    #.test_Compare_NaSVT1_NaSVT2(a6, a3, naa6, naa3)
    #.test_Compare_NaSVT1_NaSVT2(a4, a5, naa4, naa5)
    #.test_Compare_NaSVT1_NaSVT2(a5, a4, naa5, naa4)

    ## "double" <op> "logical".
    .test_Compare_NaSVT1_NaSVT2(a3, a8, naa3, naa8)
    .test_Compare_NaSVT1_NaSVT2(a8, a3, naa8, naa3)
    .test_Compare_NaSVT1_NaSVT2(a4, a7, naa4, naa7)
    .test_Compare_NaSVT1_NaSVT2(a7, a4, naa7, naa4)

    ## "double" <op> "complex".
    .test_Compare_NaSVT1_NaSVT2(a3, a10, naa3, naa10)
    .test_Compare_NaSVT1_NaSVT2(a10, a3, naa10, naa3)
    .test_Compare_NaSVT1_NaSVT2(a4, a9, naa4, naa9)
    .test_Compare_NaSVT1_NaSVT2(a9, a4, naa9, naa4)

    ## "raw" <op> "logical".
    #.test_Compare_NaSVT1_NaSVT2(a5, a8, naa5, naa8)
    #.test_Compare_NaSVT1_NaSVT2(a8, a5, naa8, naa5)
    #.test_Compare_NaSVT1_NaSVT2(a6, a7, naa6, naa7)
    #.test_Compare_NaSVT1_NaSVT2(a7, a6, naa7, naa6)

    ## "raw" <op> "complex".
    #.test_Compare_NaSVT1_NaSVT2(a5, a10, naa5, naa10)
    #.test_Compare_NaSVT1_NaSVT2(a10, a5, naa10, naa5)
    #.test_Compare_NaSVT1_NaSVT2(a6, a9, naa6, naa9)
    #.test_Compare_NaSVT1_NaSVT2(a9, a6, naa9, naa6)

    ## "logical" <op> "complex".
    .test_Compare_NaSVT1_NaSVT2(a7, a10, naa7, naa10)
    .test_Compare_NaSVT1_NaSVT2(a10, a7, naa10, naa7)
    .test_Compare_NaSVT1_NaSVT2(a8, a9, naa8, naa9)
    .test_Compare_NaSVT1_NaSVT2(a9, a8, naa9, naa8)

    ## Not expected to work.
    expect_error(naa1 != naa2[ , , -1], "non-conformable")
    expect_error(naa1 < naa2[ , , -1], "non-conformable")
    expect_error(naa1 > naa2[ , , -1], "non-conformable")
})

