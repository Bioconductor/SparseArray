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
           dimnames=list(NULL, paste0("A1y", 1:5), NULL)),
    array(101:240, c(7, 5, 4),
           dimnames=list(paste0("A2x", 1:7), paste0("A2y", 1:5), NULL)),
    array(10001:10100, c(5, 5, 4),
           dimnames=list(paste0("A3x", 1:5), NULL, paste0("A3z", 1:4)))
)

test_that("rbind() on SVT_SparseMatrix objects", {
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

test_that("cbind() on SVT_SparseMatrix objects", {
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

test_that("arbind() on 3D SVT_SparseArray objects", {
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

test_that("acbind() on 3D SVT_SparseArray objects", {
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

if (SparseArray:::.SVT_VERSION != 0L) {

test_that("handling of lacunar leaves in .abind_SVT_SparseArray_objects()", {
    m1 <- matrix(c(0:2, 0L, 0L, 0L, 0L, 1L, 0L), ncol=3)
    m2 <- matrix(c(1L, 1L, 0L, 1L, 0L, 1L, 3:2, 0L), ncol=3)
    svt1 <- as(m1, "SVT_SparseMatrix")
    svt2 <- as(m2, "SVT_SparseMatrix")

    m <- rbind(m1, m1)
    svt <- rbind(svt1, svt1)
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    expect_identical(svt@SVT[[3L]], make_lacunar_leaf("integer", c(1L, 4L)))
    expect_identical(svt, as(m, "SVT_SparseMatrix"))

    m <- rbind(m2, m2)
    svt <- rbind(svt2, svt2)
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    expect_identical(svt@SVT[[1L]], make_lacunar_leaf("integer",
                                                      c(0:1, 3:4)))
    expect_identical(svt@SVT[[2L]], make_lacunar_leaf("integer",
                                                      c(0L, 2:3, 5L)))
    expect_identical(svt, as(m, "SVT_SparseMatrix"))

    m <- rbind(m1, m2)
    svt <- rbind(svt1, svt2)
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    expect_identical(svt@SVT[[2L]], make_lacunar_leaf("integer", c(3L, 5L)))
    expect_identical(svt, as(m, "SVT_SparseMatrix"))

    m <- rbind(m2, m1)
    svt <- rbind(svt2, svt1)
    check_SparseArray_object(svt, "SVT_SparseMatrix", m)
    expect_identical(svt@SVT[[2L]], make_lacunar_leaf("integer", c(0L, 2L)))
    expect_identical(svt, as(m, "SVT_SparseMatrix"))
})

}  # ----- end if (SparseArray:::.SVT_VERSION != 0L) -----

test_that("abind() default method on SparseArray objects", {
    a1 <- .TEST_arrays[[1]]
    a2 <- .TEST_arrays[[2]]
    a3 <- .TEST_arrays[[3]]
    svt_objects <- lapply(.TEST_arrays, as, "SVT_SparseArray")
    svt1 <- svt_objects[[1]]
    svt2 <- svt_objects[[2]]
    svt3 <- svt_objects[[3]]

    ## The default abind() method is defined in the S4Arrays package. It will
    ## be called if the input is a list object or if the supplied objects are
    ## a mix of SparseArray objects and ordinary arrays.

    ## --- Input is a list ---

    a <- abind(.TEST_arrays, along=1)
    svt <- abind(svt_objects, along=1)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    expect_identical(svt, as(a, "SVT_SparseArray"))

    expected_words <- c("all", "objects", "must", "be",
                        "supplied", "via", "list")
    regexp <- paste0("\\b", expected_words, "\\b", collapse=".*")
    expect_error(abind(svt1, svt_objects), regexp, ignore.case=TRUE)
    expect_error(abind(svt_objects, svt1), regexp, ignore.case=TRUE)

    ## --- Input is a mix of SparseArray objects and ordinary arrays ---

    expected <- abind(svt_objects, along=1)
    current <- abind(a1, svt2, a3, along=1)
    expect_identical(current, expected)

    m2 <- .TEST_matrices[[2]]
    a <- abind(m2, a2)
    svt <- abind(m2, svt2)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    expect_identical(svt, as(a, "SVT_SparseArray"))
})

test_that("abind(..., rev.along=0) on SparseArray objects", {
    a1 <- .TEST_arrays[[1]]
    a2 <- .TEST_arrays[[2]][1:3, , ]
    a3 <- .TEST_arrays[[3]][1:3, , ]
    svt1 <- as(a1, "SVT_SparseArray")
    svt2 <- as(a2, "SVT_SparseArray")
    svt3 <- as(a3, "SVT_SparseArray")

    a <- abind(a1, a2, a3, rev.along=0)
    svt <- abind(svt1, svt2, svt3, rev.along=0)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    expect_identical(svt, as(a, "SVT_SparseArray"))
    expect_identical(abind(svt1, svt2, svt3, along=4), svt)
})

