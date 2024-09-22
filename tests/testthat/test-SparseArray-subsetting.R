
test_that("SVT_SparseArray subsetting by an Mindex or Lindex", {
    ## --- 3D ---
    a0 <- array(0L, c(7, 10, 3),
                dimnames=list(NULL, letters[1:10], LETTERS[1:3]))
    a0[ , 2, 1] <- a0[c(1:2, 6), 4, 1] <- 1L
    a0[ , 8, 1] <- 81:87
    a0[ , -1, 3] <- 308:370
    Mindex2 <- rbind(c(7,  9), c(7, 10), c(6, 4), c(2, 4), c(1, 10),
                     c(7, 10), c(1,  1), c(5, 4), c(2, 4))
    Mindex3 <- rbind(cbind(Mindex2, 1),
                     cbind(Mindex2, 2),
                     cbind(Mindex2, 3))
    svt0 <- as(a0, "SVT_SparseArray")
    test_linear_subsetting(a0, svt0, Mindex3)

    ## --- 2D ---
    m0 <- a0[ , , 1]
    svt0 <- as(m0, "SVT_SparseArray")
    test_linear_subsetting(m0, svt0, Mindex2)

    ## --- 1D ---
    x0 <- as.array(m0[1, ])
    Mindex1 <- Mindex2[ , -2, drop=FALSE]
    svt0 <- as(x0, "SVT_SparseArray")
    test_linear_subsetting(x0, svt0, Mindex1)
})

test_that("SVT_SparseArray subsetting by an Nindex", {
    ## --- 3D ---
    a0 <- array(0L, c(7, 10, 3),
                dimnames=list(NULL, letters[1:10], LETTERS[1:3]))
    a0[ , 2, 1] <- a0[c(1:2, 6), 4, 1] <- 1L
    a0[ , 8, 1] <- 81:87
    a0[ , -1, 3] <- 308:370
    svt0 <- as(a0, "SVT_SparseArray")

    expect_identical(svt0[ , , ], svt0)

    a   <- a0  [ , c(4:3, 8), 1, drop=FALSE]
    svt <- svt0[ , c(4:3, 8), 1, drop=FALSE]
    check_SVT_SparseArray_object(svt, a)
    m   <- a0  [ , c(4:3, 8), 1, drop=TRUE]
    svt <- svt0[ , c(4:3, 8), 1, drop=TRUE]
    check_SVT_SparseArray_object(svt, m)

    a   <- a0  [7, , , drop=FALSE]
    svt <- svt0[7, , , drop=FALSE]
    check_SVT_SparseArray_object(svt, a)
    m   <- a0  [7, , , drop=TRUE]
    svt <- svt0[7, , , drop=TRUE]
    check_SVT_SparseArray_object(svt, m)

    ## --- 2D ---
    m0 <- a0[ , , 1]
    svt0 <- as(m0, "SVT_SparseMatrix")

    expect_identical(svt0[ , ], svt0)

    m   <- m0  [-5 , c(4:3, 8)]
    svt <- svt0[-5 , c(4:3, 8)]
    check_SVT_SparseArray_object(svt, m)

    expect_identical(svt0[ , 4], m0[ , 4])
    m   <- m0  [ , 4, drop=FALSE]
    svt <- svt0[ , 4, drop=FALSE]
    check_SVT_SparseArray_object(svt, m)

    expect_identical(svt0[6 , ], m0[6 , ])
    m   <- m0  [6, -1, drop=FALSE]
    svt <- svt0[6, -1, drop=FALSE]
    check_SVT_SparseArray_object(svt, m)

    ## --- 1D ---
    x0 <- as.array(m0[6, ])
    svt0 <- as(x0, "SVT_SparseArray")

    expect_identical(svt0[ ], svt0)

    x   <- x0  [c(8:4, 1, 4), drop=FALSE]
    svt <- svt0[c(8:4, 1, 4), drop=FALSE]
    check_SVT_SparseArray_object(svt, x)
    x   <- x0  [-4, drop=FALSE]
    svt <- svt0[-4, drop=FALSE]
    check_SVT_SparseArray_object(svt, x)
    x   <- x0  [-4, drop=FALSE]
    svt <- svt0[-4, drop=FALSE]
    check_SVT_SparseArray_object(svt, x)

    subscript <- c("d", "j", "j", "h")
    x   <- x0  [subscript, drop=FALSE]  # 'drop=TRUE' would do the same thing!
    svt <- svt0[subscript, drop=FALSE]
    check_SVT_SparseArray_object(svt, x)
    expect_identical(svt0[subscript], S4Arrays:::drop_even_if_1D(x))
})

