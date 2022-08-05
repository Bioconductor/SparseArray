.TEST_matrices <- list(
    matrix(1:15, nrow=3, ncol=5,
           dimnames=list(NULL, paste0("M1y", 1:5))),
    matrix(101:135, nrow=7, ncol=5,
           dimnames=list(paste0("M2x", 1:7), paste0("M2y", 1:5))),
    matrix(1001:1025, nrow=5, ncol=5,
           dimnames=list(paste0("M3x", 1:5), NULL))
)

.TEST_arrays <- list(
    array(1:60, c(3, 5, 4),
           dimnames=list(NULL, paste0("M1y", 1:5), NULL)),
    array(101:240, c(7, 5, 4),
           dimnames=list(paste0("M2x", 1:7), paste0("M2y", 1:5), NULL)),
    array(10001:10100, c(5, 5, 4),
           dimnames=list(paste0("M3x", 1:5), NULL, paste0("M3z", 1:4)))
)

test_that("rbind_SVT_SparseMatrix_objects", {
    m1 <- .TEST_matrices[[1]]
    m2 <- .TEST_matrices[[2]]
    m3 <- .TEST_matrices[[3]]
    svt1 <- as(m1, "SVT_SparseMatrix")
    svt2 <- as(m2, "SVT_SparseMatrix")
    svt3 <- as(m3, "SVT_SparseMatrix")

    m <- rbind(a=m1, b=m2, c=m3)
    svt <- rbind(a=svt1, b=svt2, c=svt3)
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    expect_identical(svt, as(m, "SVT_SparseMatrix"))

    ## Unary form.
    expect_identical(rbind(a=svt1), svt1)

    ## Zero rows.
    m1 <- matrix(nrow=0, ncol=3, dimnames=list(NULL, letters[1:3]))
    m2 <- matrix(1:15, ncol=3, dimnames=list(NULL, LETTERS[1:3]))
    svt1 <- as(m1, "SVT_SparseMatrix")
    svt2 <- as(m2, "SVT_SparseMatrix")

    m <- rbind(a=m1, b=m2)
    svt <- rbind(a=svt1, b=svt2)
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    expect_identical(svt, as(m, "SVT_SparseMatrix"))

    m <- rbind(a=m2, b=m1)
    svt <- rbind(a=svt2, b=svt1)
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    expect_identical(svt, as(m, "SVT_SparseMatrix"))

    m <- rbind(a=m1, b=m1)
    svt <- rbind(a=svt1, b=svt1)
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    expect_identical(svt, as(m, "SVT_SparseMatrix"))
})

test_that("cbind_SVT_SparseMatrix_objects", {
    m1 <- t(.TEST_matrices[[1]])
    m2 <- t(.TEST_matrices[[2]])
    m3 <- t(.TEST_matrices[[3]])
    svt1 <- as(m1, "SVT_SparseMatrix")
    svt2 <- as(m2, "SVT_SparseMatrix")
    svt3 <- as(m3, "SVT_SparseMatrix")

    m <- cbind(a=m1, b=m2, c=m3)
    svt <- cbind(a=svt1, b=svt2, c=svt3)
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    expect_identical(svt, as(m, "SVT_SparseMatrix"))

    ## Unary form.
    expect_identical(cbind(a=svt1), svt1)

    ## Zero cols.
    m1 <- matrix(nrow=3, ncol=0, dimnames=list(letters[1:3], NULL))
    m2 <- matrix(1:15, nrow=3, dimnames=list(LETTERS[1:3], NULL))
    svt1 <- as(m1, "SVT_SparseMatrix")
    svt2 <- as(m2, "SVT_SparseMatrix")

    m <- cbind(a=m1, b=m2)
    svt <- cbind(a=svt1, b=svt2)
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    expect_identical(svt, as(m, "SVT_SparseMatrix"))

    m <- cbind(a=m2, b=m1)
    svt <- cbind(a=svt2, b=svt1)
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    expect_identical(svt, as(m, "SVT_SparseMatrix"))

    m <- cbind(a=m1, b=m1)
    svt <- cbind(a=svt1, b=svt1)
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    expect_identical(svt, as(m, "SVT_SparseMatrix"))
})

test_that("arbind_3D_SVT_SparseArray_objects", {
    a1 <- .TEST_arrays[[1]]
    a2 <- .TEST_arrays[[2]]
    a3 <- .TEST_arrays[[3]]
    svt1 <- as(a1, "SVT_SparseArray")
    svt2 <- as(a2, "SVT_SparseArray")
    svt3 <- as(a3, "SVT_SparseArray")

    a <- arbind(a=a1, b=a2, c=a3)
    svt <- arbind(a=svt1, b=svt2, c=svt3)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    expect_identical(svt, as(a, "SVT_SparseArray"))

    ## Unary form.
    expect_identical(arbind(a=svt1), svt1)

    ## Zero ROWS.
    a1 <- a1[0, , ]
    svt1 <- as(a1, "SVT_SparseArray")

    a <- arbind(a=a1, b=a2)
    svt <- arbind(a=svt1, b=svt2)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    expect_identical(svt, as(a, "SVT_SparseArray"))

    a <- arbind(a=a2, b=a1)
    svt <- arbind(a=svt2, b=svt1)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    expect_identical(svt, as(a, "SVT_SparseArray"))

    a <- arbind(a=a1, b=a1)
    svt <- arbind(a=svt1, b=svt1)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    expect_identical(svt, as(a, "SVT_SparseArray"))
})

test_that("acbind_3D_SVT_SparseArray_objects", {
    a1 <- aperm(.TEST_arrays[[1]], c(2:1, 3))
    a2 <- aperm(.TEST_arrays[[2]], c(2:1, 3))
    a3 <- aperm(.TEST_arrays[[3]], c(2:1, 3))
    svt1 <- as(a1, "SVT_SparseArray")
    svt2 <- as(a2, "SVT_SparseArray")
    svt3 <- as(a3, "SVT_SparseArray")

    a <- acbind(a=a1, b=a2, c=a3)
    svt <- acbind(a=svt1, b=svt2, c=svt3)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    expect_identical(svt, as(a, "SVT_SparseArray"))

    ## Unary form.
    expect_identical(acbind(a=svt1), svt1)

    ## Zero COLS.
    a1 <- a1[ , 0, ]
    svt1 <- as(a1, "SVT_SparseArray")

    a <- acbind(a=a1, b=a2)
    svt <- acbind(a=svt1, b=svt2)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    expect_identical(svt, as(a, "SVT_SparseArray"))

    a <- acbind(a=a2, b=a1)
    svt <- acbind(a=svt2, b=svt1)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    expect_identical(svt, as(a, "SVT_SparseArray"))

    a <- acbind(a=a1, b=a1)
    svt <- acbind(a=svt1, b=svt1)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    expect_identical(svt, as(a, "SVT_SparseArray"))
})

