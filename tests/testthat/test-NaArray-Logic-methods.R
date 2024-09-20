
.test_Logic_NaSVT1_v2 <- function(a1, naa1, y)
{
    ## naa1 & y
    if (isFALSE(y)) {
        expect_error(naa1 & y, "not supported")
        expect_error(y & naa1, "not supported")
    } else {
        a <- a1 & y
        naa <- naa1 & y
        check_NaArray_object(naa, a)
        naa <- y & naa1
        check_NaArray_object(naa, a)
    }

    ## naa1 | y
    if (isTRUE(y)) {
	expect_error(naa1 | y, "not supported")
	expect_error(y | naa1, "not supported")
    } else {
        a <- a1 | y
        naa <- naa1 | y
        check_NaArray_object(naa, a)
        naa <- y | naa1
        check_NaArray_object(naa, a)
    }
}

### We also test Logic ops between an NaArray and SVT_SparseArray object!
.test_Logic_NaSVT1_NaSVT2 <- function(a1, a2, naa1, naa2)
{
    svt1 <- as(a1, "SVT_SparseArray")
    svt2 <- as(a2, "SVT_SparseArray")

    ## naa1 & naa2
    a <- a1 & a2
    naa <- naa1 & naa2
    check_NaArray_object(naa, a)
    svt <- naa1 & svt2
    check_SVT_SparseArray_object(svt, a)
    expect_identical(svt, svt1 & naa2)

    ## naa1 | naa2
    a <- a1 | a2
    naa <- naa1 | naa2
    check_NaArray_object(naa, a)
    expect_identical(naa, naa1 | svt2)
    expect_identical(naa, svt1 | naa2)
}

test_that("'Logic' ops between NaArray object and single value", {
    a0 <- matrix(NA, nrow=4, ncol=7)
    naa0 <- as(a0, "NaArray")
    .test_Logic_NaSVT1_v2(a0, naa0, TRUE)
    .test_Logic_NaSVT1_v2(a0, naa0, FALSE)
    .test_Logic_NaSVT1_v2(a0, naa0, NA)

    a1 <- make_3D_logical_array(NA)
    naa1 <- as(a1, "NaArray")
    .test_Logic_NaSVT1_v2(a1, naa1, TRUE)
    .test_Logic_NaSVT1_v2(a1, naa1, FALSE)
    .test_Logic_NaSVT1_v2(a1, naa1, NA)
})

test_that("'Logic' ops between 2 NaArray objects", {
    a1 <- array(NA, dim=c(1, 4:2),
                dimnames=list(NULL, letters[1:4], LETTERS[1:3], NULL))
    a2 <- array(FALSE, dim=dim(a1),
                dimnames=list("A", NULL, letters[24:26], c("Y", "Z")))
    naa1 <- as(a1, "NaArray")
    naa2 <- as(a2, "NaArray")
    .test_Logic_NaSVT1_NaSVT2(a1, a2, naa1, naa2)
    .test_Logic_NaSVT1_NaSVT2(a2, a1, naa2, naa1)
    .test_Logic_NaSVT1_NaSVT2(a1, a1, naa1, naa1)
    .test_Logic_NaSVT1_NaSVT2(a2, a2, naa2, naa2)

    a1 <- a1[0 , , , 0]
    a2 <- a2[0 , , , 0]
    naa1 <- as(a1, "NaArray")
    naa2 <- as(a2, "NaArray")
    .test_Logic_NaSVT1_NaSVT2(a1, a2, naa1, naa2)
    .test_Logic_NaSVT1_NaSVT2(a2, a1, naa2, naa1)
    .test_Logic_NaSVT1_NaSVT2(a1, a1, naa1, naa1)
    .test_Logic_NaSVT1_NaSVT2(a2, a2, naa2, naa2)

    a1 <- make_3D_logical_array(NA)
    a2 <- array(FALSE, dim(a1),
                dimnames=list(NULL, letters[22:26], LETTERS[23:26]))
    a2[c(1, 2, 6), 2, 1] <- a2[-3, 5, 1] <-TRUE
    a2[2:6, 3, 1] <- NA
    a2[c(1:3, 6), , 2] <- c(TRUE, NA, NA, TRUE, NA)
    a2[ , , 3] <- c(TRUE, NA)
    naa1 <- as(a1, "NaArray")
    naa2 <- as(a2, "NaArray")
    .test_Logic_NaSVT1_NaSVT2(a1, a2, naa1, naa2)
    .test_Logic_NaSVT1_NaSVT2(a2, a1, naa2, naa1)
    .test_Logic_NaSVT1_NaSVT2(a1, a1, naa1, naa1)
    .test_Logic_NaSVT1_NaSVT2(a2, a2, naa2, naa2)

    ## Not expected to work.
    expect_error(naa1 & naa2[ , , -1], "non-conformable")
    expect_error(naa1 | naa2[ , , -1], "non-conformable")
})

