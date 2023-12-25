
test_that(".aperm_SVT()", {
    ## --- with 2 dimensions ---

    m0 <- array(1:120, c(8, 15))
    svt0 <- SparseArray(m0)

    expect_identical(aperm(svt0, 1:2), svt0)

    expected <- aperm(m0)
    current <- aperm(svt0)
    check_SparseArray_object(current, "SVT_SparseMatrix", expected)
    expect_identical(aperm(current, 2:1), svt0)

    ## --- with 3 dimensions ---

    a0 <- array(1:360, c(8, 3, 15))
    a0[ , , c(11,15)] <- 0L
    a0[ , 1:2, 14] <- 0L
    a0[c(1:4, 6), 3, 13] <- 0L
    svt0 <- SparseArray(a0)

    expect_identical(aperm(svt0, 1:3), svt0)

    expected <- aperm(a0)
    current <- aperm(svt0)
    check_SparseArray_object(current, "SVT_SparseArray", expected)
    expect_identical(aperm(current, 3:1), svt0)

    perm <- c(1, 3, 2)
    expected <- aperm(a0, perm)
    current <- aperm(svt0, perm)
    check_SparseArray_object(current, "SVT_SparseArray", expected)
    expect_identical(aperm(current, perm), svt0)

    perm <- c(2, 1, 3)
    expected <- aperm(a0, perm)
    current <- aperm(svt0, perm)
    check_SparseArray_object(current, "SVT_SparseArray", expected)
    expect_identical(aperm(current, perm), svt0)

    perm <- c(2, 3, 1)
    expected <- aperm(a0, perm)
    current <- aperm(svt0, perm)
    check_SparseArray_object(current, "SVT_SparseArray", expected)
    expect_identical(aperm(current, perm[perm]), svt0)

    ## --- with 5 dimensions and dimnames ---

    a0 <- array(1:184800, c(55, 8, 10, 3, 14),
                dimnames=list(NULL, NULL, letters[1:10], NULL, LETTERS[1:14]))
    svt0 <- SparseArray(a0)

    expect_identical(aperm(svt0, 1:5), svt0)

    expected <- aperm(a0)
    current <- aperm(svt0)
    check_SparseArray_object(current, "SVT_SparseArray", expected)
    expect_identical(aperm(current, 5:1), svt0)

    perm <- c(1, 2, 3, 5, 4)
    expected <- aperm(a0, perm)
    current <- aperm(svt0, perm)
    check_SparseArray_object(current, "SVT_SparseArray", expected)
    expect_identical(aperm(current, perm), svt0)

    perm <- c(1, 2, 4, 3, 5)
    expected <- aperm(a0, perm)
    current <- aperm(svt0, perm)
    check_SparseArray_object(current, "SVT_SparseArray", expected)
    expect_identical(aperm(current, perm), svt0)

    perm <- c(3, 1, 2, 4, 5)
    expected <- aperm(a0, perm)
    current <- aperm(svt0, perm)
    check_SparseArray_object(current, "SVT_SparseArray", expected)
    rperm <- c(2, 3, 1, 4, 5)  # reverse permutation
    expect_identical(aperm(current, rperm), svt0)

    perm <- c(1, 2, 5, 3, 4)
    expected <- aperm(a0, perm)
    current <- aperm(svt0, perm)
    check_SparseArray_object(current, "SVT_SparseArray", expected)
    rperm <- c(1, 2, 4, 5, 3)  # reverse permutation
    expect_identical(aperm(current, rperm), svt0)

    perm <- c(4, 3, 5, 1, 2)
    expected <- aperm(a0, perm)
    current <- aperm(svt0, perm)
    check_SparseArray_object(current, "SVT_SparseArray", expected)
    rperm <- c(4, 5, 2, 1, 3)  # reverse permutation
    expect_identical(aperm(current, rperm), svt0)

    expect_identical(aperm(svt0, perm[perm]), aperm(current, perm))

    a <- a0 * 0.1
    svt <- SparseArray(a)

    expect_identical(aperm(svt, 1:5), svt)

    expected <- aperm(a)
    current <- aperm(svt)
    check_SparseArray_object(current, "SVT_SparseArray", expected)
    expect_identical(aperm(current, 5:1), svt)

    perm <- c(1, 2, 5, 3, 4)
    expected <- aperm(a, perm)
    current <- aperm(svt, perm)
    check_SparseArray_object(current, "SVT_SparseArray", expected)
    rperm <- c(1, 2, 4, 5, 3)  # reverse permutation
    expect_identical(aperm(current, rperm), svt)

    perm <- c(4, 3, 5, 1, 2)
    expected <- aperm(a, perm)
    current <- aperm(svt, perm)
    check_SparseArray_object(current, "SVT_SparseArray", expected)
    rperm <- c(4, 5, 2, 1, 3)  # reverse permutation
    expect_identical(aperm(current, rperm), svt)

    expect_identical(aperm(svt, perm[perm]), aperm(current, perm))
})

