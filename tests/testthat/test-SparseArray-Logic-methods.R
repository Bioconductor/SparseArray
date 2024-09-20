
.test_Logic_SVT1_v2 <- function(a1, svt1, y)
{
    ## svt1 & y
    a <- a1 & y
    svt <- svt1 & y
    check_SVT_SparseArray_object(svt, a)
    svt <- y & svt1
    check_SVT_SparseArray_object(svt, a)

    ## svt1 | y
    if (isTRUE(y)) {
	expect_error(svt1 | y, "not supported")
	expect_error(y | svt1, "not supported")
    } else {
        a <- a1 | y
        svt <- svt1 | y
        check_SVT_SparseArray_object(svt, a)
        svt <- y | svt1
        check_SVT_SparseArray_object(svt, a)
    }
}

.test_Logic_SVT1_SVT2 <- function(a1, a2, svt1, svt2)
{
    ## svt1 & svt2
    a <- a1 & a2
    svt <- svt1 & svt2
    check_SVT_SparseArray_object(svt, a)

    ## svt1 | svt2
    a <- a1 | a2
    svt <- svt1 | svt2
    check_SVT_SparseArray_object(svt, a)
}

test_that("'Logic' ops between SVT_SparseArray object and single value", {
    a1 <- array(FALSE, 6:4)
    dimnames(a1) <- list(letters[1:6], NULL, LETTERS[1:4])
    a1[c(2:3, 6), 2, 1] <- TRUE
    a1[c(1, 6), 1 , 2] <- NA
    a1[1:5, 5, 3] <- TRUE
    a1[1, 3, 2] <- NA

    svt1 <- as(a1, "SVT_SparseArray")
    .test_Logic_SVT1_v2(a1, svt1, TRUE)
    .test_Logic_SVT1_v2(a1, svt1, FALSE)

    ## Not expected to work.
    expect_error(svt1 & NA, "not supported")
    expect_error(NA & svt1, "not supported")
    expect_error(svt1 | NA, "not supported")
    expect_error(NA | svt1, "not supported")
})

test_that("'Logic' ops between 2 SVT_SparseArray objects", {
    a1 <- array(FALSE, 6:4)
    dimnames(a1) <- list(letters[1:6], NULL, LETTERS[1:4])
    a1[c(2:3, 6), 2, 1] <- TRUE
    a1[c(1, 6), 1 , 2] <- NA
    a1[1:5, 5, 3] <- TRUE
    a1[1, 3, 2] <- NA

    a2 <- a1
    dimnames(a2) <- list(NULL, letters[22:26], LETTERS[23:26])
    a2[c(1, 2, 6), 2, 1] <- TRUE
    a2[2:6, 3, 1] <- NA
    a2[ , 3, 2] <- c(TRUE, NA)

    svt1 <- as(a1, "SVT_SparseArray")
    svt2 <- as(a2, "SVT_SparseArray")
    .test_Logic_SVT1_SVT2(a1, a2, svt1, svt2)
    .test_Logic_SVT1_SVT2(a2, a1, svt2, svt1)
    .test_Logic_SVT1_SVT2(a1, a1, svt1, svt1)
    .test_Logic_SVT1_SVT2(a2, a2, svt2, svt2)

    ## Not expected to work.
    expect_error(svt1 & svt2[ , , -1], "non-conformable")
    expect_error(svt1 | svt2[ , , -1], "non-conformable")
})

