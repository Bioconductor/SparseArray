
.test_Arith_SparseArray <- function(a1, a2, svt1, svt2)
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

test_that("\"+\", \"-\", \"*\" between SVT_SparseArray objects", {
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

    .test_Arith_SparseArray(a1, a2, svt1, svt2)
    .test_Arith_SparseArray(a2, a1, svt2, svt1)
    .test_Arith_SparseArray(a1, a1, svt1, svt1)
    .test_Arith_SparseArray(a2, a2, svt2, svt2)
    expect_identical((svt1 + svt2) - svt1, `dimnames<-`(svt2, dimnames(a1)))
    expect_identical(svt1 + (svt2 - svt1), `dimnames<-`(svt2, dimnames(a1)))
    expect_identical(svt1 - (svt1 - svt2), `dimnames<-`(svt2, dimnames(a1)))
    .test_Arith_SparseArray(a1, a1 + a2, svt1, svt1 + svt2)
    .test_Arith_SparseArray(a1, a2 - a1, svt1, svt2 - svt1)
    .test_Arith_SparseArray(a2, a1 - a2, svt2, svt1 - svt2)
    .test_Arith_SparseArray(a1, a1 * a2, svt1, svt1 * svt2)

    dimnames(a1) <- dimnames(svt1) <- NULL
    .test_Arith_SparseArray(a1, a2, svt1, svt2)
    .test_Arith_SparseArray(a2, a1, svt2, svt1)

    a0 <- a2[ , 0, ]
    svt0 <- as(a0, "SVT_SparseArray")
    .test_Arith_SparseArray(a0, a0, svt0, svt0)
    .test_Arith_SparseArray(a0, unname(a0), svt0, unname(svt0))
    .test_Arith_SparseArray(unname(a0), a0, unname(svt0), svt0)

    ## Not expected to work.
    expect_error(svt1 + svt2[ , , -1], "non-conformable")
    expect_error(svt1 - svt2[ , , -1], "non-conformable")
    expect_error(svt1 * svt2[ , , -1], "non-conformable")
    expect_error(svt1 / svt2, regexp="not implemented")
    expect_error(svt1 ^ svt2, regexp="not implemented")
    expect_error(svt1 %% svt2, regexp="not implemented")
    expect_error(svt1 %/% svt2, regexp="not implemented")
})

