
test_that(".aperm_SVT()", {
    .aperm_SVT <- SparseArray:::.aperm_SVT

    a <- array(1:198000, c(55, 8, 10, 3, 15))
    svt <- SparseArray(a)
    perm <- c(4, 3, 5, 1, 2)
    expected <- aperm(a, perm)
    current <- aperm(svt, perm)
    expect_identical(as.array(current), expected)
    rperm <- c(4, 5, 2, 1, 3)  # reverse permutation
    expect_identical(aperm(current, rperm), svt)

    expected <- aperm(a)
    current <- aperm(svt)
    expect_identical(as.array(current), expected)
    expect_identical(aperm(current), svt)

    a <- a * 0.1
    svt <- SparseArray(a) 
    perm <- c(4, 3, 5, 1, 2)
    expected <- aperm(a, perm)
    current <- aperm(svt, perm)
    expect_identical(as.array(current), expected)
    rperm <- c(4, 5, 2, 1, 3)  # reverse permutation
    expect_identical(aperm(current, rperm), svt)

    expected <- aperm(a)
    current <- aperm(svt)
    expect_identical(as.array(current), expected)
    expect_identical(aperm(current), svt)
})

