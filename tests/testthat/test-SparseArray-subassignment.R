.test_SparseArray_subassignment_by_Mindex_and_Lindex <-
    function(a0, Mindex, vals, expected_class)
{
    object0 <- as(a0, expected_class)
    Lindex <- Mindex2Lindex(Mindex, dim(a0))

    a <- `[<-`(a0, Mindex, value=vals)
    object <- `[<-`(object0, Mindex, value=vals)
    check_array_like_object(object, expected_class, a)
    object <- `[<-`(object0, Lindex, value=vals)
    check_array_like_object(object, expected_class, a)
    object <- `[<-`(object0, as.double(Lindex), value=vals)
    check_array_like_object(object, expected_class, a)
    object <- `[<-`(object0, Lindex + 0.5, value=vals)
    check_array_like_object(object, expected_class, a)
}

test_that("subassign an SVT_SparseArray object by an Mindex or Lindex", {
    ## Only zeros.
    a0 <- array(0L, c(7, 10, 3),
                dimnames=list(NULL, letters[1:10], LETTERS[1:3]))
    Mindex3 <- rbind(c(7,  9, 3), c(7, 10, 3), c(6, 4, 3), c(2, 4, 3),
                     c(1, 10, 3), c(7, 10, 3), c(1, 1, 3), c(5, 4, 3),
                     c(2,  4, 3))
    vals <- c(11:18, 0L)
    .test_SparseArray_subassignment_by_Mindex_and_Lindex(a0, Mindex3, vals,
                                                         "SVT_SparseArray")
    m0 <- a0[ , , 1]  # 2D
    Mindex2 <- Mindex3[ , -3]
    .test_SparseArray_subassignment_by_Mindex_and_Lindex(m0, Mindex2, vals,
                                                         "SVT_SparseMatrix")
    x0 <- as.array(m0[1, ])  # 1D
    Mindex1 <- Mindex2[ , -2, drop=FALSE]
    .test_SparseArray_subassignment_by_Mindex_and_Lindex(x0, Mindex1, vals,
                                                         "SVT_SparseArray")

    ## Add some nonzero elements.
    a0 <- make_3D_double_array()
    Mindex23 <- rbind(cbind(Mindex2, 1L), Mindex3)
    vals2 <- c(vals, vals)
    Mindex0 <- nzwhich(a0, arr.ind=TRUE)
    .test_SparseArray_subassignment_by_Mindex_and_Lindex(a0, Mindex23, vals2,
                                                         "SVT_SparseArray")
    .test_SparseArray_subassignment_by_Mindex_and_Lindex(a0, Mindex0, 0,
                                                         "SVT_SparseArray")
    m0 <- a0[ , , 1]  # 2D
    .test_SparseArray_subassignment_by_Mindex_and_Lindex(m0, Mindex2, vals,
                                                         "SVT_SparseMatrix")
    x0 <- as.array(m0[1, ])  # 1D
    .test_SparseArray_subassignment_by_Mindex_and_Lindex(x0, Mindex1, vals,
                                                         "SVT_SparseArray")

    ## Integer array.
    a0 <- make_3D_double_array()
    suppressWarnings(storage.mode(a0) <- "integer")
    .test_SparseArray_subassignment_by_Mindex_and_Lindex(a0, Mindex23, vals2,
                                                         "SVT_SparseArray")
    .test_SparseArray_subassignment_by_Mindex_and_Lindex(a0, Mindex0, 0L,
                                                         "SVT_SparseArray")

    ## Array type changed by subassignment.
    a0 <- make_3D_double_array()
    vals2 <- complex(real=vals2, imaginary=-0.75)
    .test_SparseArray_subassignment_by_Mindex_and_Lindex(a0, Mindex23, vals2,
                                                         "SVT_SparseArray")
    .test_SparseArray_subassignment_by_Mindex_and_Lindex(a0, Mindex0, -9.99,
                                                         "SVT_SparseArray")

    ## Assign random values to random array locations.
    set.seed(123)
    Mindex <- Lindex2Mindex(sample(length(a0)), dim(a0))
    vals <- sample(0:5, length(a0), replace=TRUE)
    .test_SparseArray_subassignment_by_Mindex_and_Lindex(a0, Mindex, vals,
                                                         "SVT_SparseArray")
    Mindex <- Lindex2Mindex(sample(length(a0), 5000, replace=TRUE), dim(a0))
    vals <- sample(-99:99, 5000, replace=TRUE)
    .test_SparseArray_subassignment_by_Mindex_and_Lindex(a0, Mindex, vals,
                                                         "SVT_SparseArray")
})

test_that("SparseArray:::.subassign_SVT_with_short_Rvector()", {
    subassign_SVT_with_short_Rvector <-
        SparseArray:::.subassign_SVT_with_short_Rvector

    test_subassign_SVT_with_short_Rvector <-
        function(svt0, Nindex, value, expected_type=type(value)) {
            svt <- subassign_SVT_with_short_Rvector(svt0, Nindex, value)
            expect_identical(type(svt), expected_type)
            a0 <- as.array(`type<-`(svt0, expected_type))
            a <- S4Arrays:::subassign_by_Nindex(a0, Nindex, value)
            expected_class <-
                if (is.matrix(a)) "SVT_SparseMatrix" else "SVT_SparseArray"
            check_array_like_object(svt, expected_class, a)
        }

    ## --- 1D objects ---

    svt0 <- SVT_SparseArray(dim=10, type="raw", dimnames=list(LETTERS[1:10]))
    x0 <- as.array(svt0)
    Nindex1 <- list(c(6:9, 2L))
    Nindex2 <- list(NULL)
    Nindex3 <- list(c(10L, 3:5, 3L))

    svt1 <- subassign_SVT_with_short_Rvector(svt0, Nindex1, as.raw(0:4))
    a1 <- S4Arrays:::subassign_by_Nindex(x0, Nindex1, as.raw(0:4))
    check_array_like_object(svt1, "SVT_SparseArray", a1)

    for (Nindex in list(Nindex1, Nindex2, Nindex3)) {
        value <- c(TRUE, FALSE, TRUE, TRUE, FALSE)
        test_subassign_SVT_with_short_Rvector(svt0, Nindex, value)
        test_subassign_SVT_with_short_Rvector(svt1, Nindex, value)
        value <- -2:2
        test_subassign_SVT_with_short_Rvector(svt0, Nindex, value)
        test_subassign_SVT_with_short_Rvector(svt1, Nindex, value)
        value <- c(-pi, NaN, 0, -Inf, NA)
        test_subassign_SVT_with_short_Rvector(svt0, Nindex, value)
        test_subassign_SVT_with_short_Rvector(svt1, Nindex, value)
        value <- 2.44 - value * 8i
        test_subassign_SVT_with_short_Rvector(svt0, Nindex, value)
        test_subassign_SVT_with_short_Rvector(svt1, Nindex, value)
        value <- c("hello", "", "world", "", "!")
        test_subassign_SVT_with_short_Rvector(svt0, Nindex, value)
        test_subassign_SVT_with_short_Rvector(svt1, Nindex, value)
        value <- list(NULL, 11:14, NULL, factor(), list("a", complex(0), NA))
        test_subassign_SVT_with_short_Rvector(svt0, Nindex, value)
        test_subassign_SVT_with_short_Rvector(svt1, Nindex, value)
    }

    ## --- 2D objects ---

    svt0 <- SVT_SparseArray(dim=c(10, 7), type="raw",
                            dimnames=list(LETTERS[1:10], letters[1:7]))
    m0 <- as.matrix(svt0)
    Nindex1 <- list(c(6:9, 2L), NULL)
    Nindex2 <- list(NULL, 2:4)
    Nindex3 <- list(c(10L, 3:5, 3L), c(6:3, 5L, 1L))
    Nindex4 <- list(c(3L, 5:2), c(1L, 7L))

    svt1 <- subassign_SVT_with_short_Rvector(svt0, Nindex1, as.raw(0:4))
    m1 <- S4Arrays:::subassign_by_Nindex(m0, Nindex1, as.raw(0:4))
    check_array_like_object(svt1, "SVT_SparseMatrix", m1)

    for (Nindex in list(Nindex1, Nindex2, Nindex3, Nindex4)) {
        value <- c(TRUE, FALSE, TRUE, TRUE, FALSE)
        test_subassign_SVT_with_short_Rvector(svt0, Nindex, value)
        test_subassign_SVT_with_short_Rvector(svt1, Nindex, value)
        value <- -2:2
        test_subassign_SVT_with_short_Rvector(svt0, Nindex, value)
        test_subassign_SVT_with_short_Rvector(svt1, Nindex, value)
        value <- c(-pi, NaN, 0, -Inf, NA)
        test_subassign_SVT_with_short_Rvector(svt0, Nindex, value)
        test_subassign_SVT_with_short_Rvector(svt1, Nindex, value)
        value <- 2.44 - value * 8i
        test_subassign_SVT_with_short_Rvector(svt0, Nindex, value)
        test_subassign_SVT_with_short_Rvector(svt1, Nindex, value)
        value <- c("hello", "", "world", "", "!")
        test_subassign_SVT_with_short_Rvector(svt0, Nindex, value)
        test_subassign_SVT_with_short_Rvector(svt1, Nindex, value)
        value <- list(NULL, 11:14, NULL, factor(), list("a", complex(0), NA))
        test_subassign_SVT_with_short_Rvector(svt0, Nindex, value)
        test_subassign_SVT_with_short_Rvector(svt1, Nindex, value)
    }

    ## --- 3D objects ---

    svt0 <- SVT_SparseArray(dim=c(10, 2, 7), type="raw",
                            dimnames=list(LETTERS[1:10], NULL, letters[1:7]))
    a0 <- as.array(svt0)
    Nindex1 <- list(c(6:9, 2L), NULL, NULL)
    Nindex2 <- list(NULL, 2L, 2:4)
    Nindex3 <- list(c(10L, 3:5, 3L), 2L, c(6:3, 5L, 1L))
    Nindex4 <- list(c(3L, 5:2), 2:1, c(1L, 7L))
    Nindex5 <- list(NULL, NULL, c(7L, 7L))

    svt1 <- subassign_SVT_with_short_Rvector(svt0, Nindex1, as.raw(0:4))
    a1 <- S4Arrays:::subassign_by_Nindex(a0, Nindex1, as.raw(0:4))
    check_array_like_object(svt1, "SVT_SparseArray", a1)

    for (Nindex in list(Nindex1, Nindex2, Nindex3, Nindex4, Nindex5)) {
        value <- c(TRUE, FALSE, TRUE, TRUE, FALSE)
        test_subassign_SVT_with_short_Rvector(svt0, Nindex, value)
        test_subassign_SVT_with_short_Rvector(svt1, Nindex, value)
        value <- -2:2
        test_subassign_SVT_with_short_Rvector(svt0, Nindex, value)
        test_subassign_SVT_with_short_Rvector(svt1, Nindex, value)
        value <- c(-pi, NaN, 0, -Inf, NA)
        test_subassign_SVT_with_short_Rvector(svt0, Nindex, value)
        test_subassign_SVT_with_short_Rvector(svt1, Nindex, value)
        value <- 2.44 - value * 8i
        test_subassign_SVT_with_short_Rvector(svt0, Nindex, value)
        test_subassign_SVT_with_short_Rvector(svt1, Nindex, value)
        value <- c("hello", "", "world", "", "!")
        test_subassign_SVT_with_short_Rvector(svt0, Nindex, value)
        test_subassign_SVT_with_short_Rvector(svt1, Nindex, value)
        value <- list(NULL, 11:14, NULL, factor(), list("a", complex(0), NA))
        test_subassign_SVT_with_short_Rvector(svt0, Nindex, value)
        test_subassign_SVT_with_short_Rvector(svt1, Nindex, value)
    }
})

test_that(".subassign_SVT_with_Rarray() and .subassign_SVT_with_SVT()", {
    test_subassign_SVT_with_Rarray_or_SVT <-
        function(svt0, Nindex, value, expected_type=type(value)) {
            svt <- SparseArray:::.subassign_SVT_with_Rarray(svt0, Nindex, value)
            expect_identical(type(svt), expected_type)
            a0 <- as.array(`type<-`(svt0, expected_type))
            a <- S4Arrays:::subassign_by_Nindex(a0, Nindex, value)
            expected_class <-
                if (is.matrix(a)) "SVT_SparseMatrix" else "SVT_SparseArray"
            check_array_like_object(svt, expected_class, a)
            value <- SVT_SparseArray(value)
            svt2 <- SparseArray:::.subassign_SVT_with_SVT(svt0, Nindex, value)
            expect_identical(svt2, svt)
        }

    ## --- 1D objects ---

    svt0 <- SVT_SparseArray(dim=10, type="integer")
    svt1 <- SVT_SparseArray(array(c(0L, 0L, 103:104, 0L, 0L, 107:109, 0L)))
    Nindex1 <- list(c(2:5, 8L))
    Nindex2 <- list(c(9L, 3:5, 3L))

    for (Nindex in list(Nindex1, Nindex2)) {
        value <- array(integer(5))
        test_subassign_SVT_with_Rarray_or_SVT(svt0, Nindex, value)
        test_subassign_SVT_with_Rarray_or_SVT(svt1, Nindex, value)
        value <- array(-2:2)
        test_subassign_SVT_with_Rarray_or_SVT(svt0, Nindex, value)
        test_subassign_SVT_with_Rarray_or_SVT(svt1, Nindex, value)
        value <- array(c(-pi, NaN, 0, -Inf, NA))
        test_subassign_SVT_with_Rarray_or_SVT(svt0, Nindex, value)
        test_subassign_SVT_with_Rarray_or_SVT(svt1, Nindex, value)
        value <- array(2.44 - value * 8i)
        test_subassign_SVT_with_Rarray_or_SVT(svt0, Nindex, value)
        test_subassign_SVT_with_Rarray_or_SVT(svt1, Nindex, value)
        value <- array(c("hello", "", "world", "", "!"))
        test_subassign_SVT_with_Rarray_or_SVT(svt0, Nindex, value)
        test_subassign_SVT_with_Rarray_or_SVT(svt1, Nindex, value)
        value <- array(list(NULL, 11:14, NULL, factor(),
                            list("a", complex(0), NA)))
        test_subassign_SVT_with_Rarray_or_SVT(svt0, Nindex, value)
        test_subassign_SVT_with_Rarray_or_SVT(svt1, Nindex, value)
    }

    value <- array(c(0L, -1:1, 0L, 2L, 0L, 0L, 3:4))
    test_subassign_SVT_with_Rarray_or_SVT(svt0, list(NULL), value)
    test_subassign_SVT_with_Rarray_or_SVT(svt1, list(NULL), value)

    ## --- 2D objects ---

    svt0 <- SVT_SparseArray(dim=c(10, 6), type="integer")
    svt1 <- `[<-`(svt0, (1:20)*3, value=1:20)
    svt1 <- `[<-`(svt1, , 2, value=0L)

    Nindex1 <- list(c(2:5, 8L), NULL)
    value1 <- -as.matrix(svt1)[6:10, ]
    test_subassign_SVT_with_Rarray_or_SVT(svt0, Nindex1, value1)
    test_subassign_SVT_with_Rarray_or_SVT(svt1, Nindex1, value1)

    Nindex2 <- list(NULL, 3:6)
    value2 <- -as.matrix(svt1)[ , 2:5]
    test_subassign_SVT_with_Rarray_or_SVT(svt0, Nindex2, value2)
    test_subassign_SVT_with_Rarray_or_SVT(svt1, Nindex2, value2)

    Nindex3 <- list(c(9L, 3:5, 3L), c(6L, 2L))
    value3 <- -as.matrix(svt1)[6:10, 5:6]
    test_subassign_SVT_with_Rarray_or_SVT(svt0, Nindex3, value3)
    test_subassign_SVT_with_Rarray_or_SVT(svt1, Nindex3, value3)

    ## --- 3D objects ---

    svt0 <- SVT_SparseArray(dim=c(6, 10, 2), type="integer")
    svt1 <- `[<-`(svt0, (1:24)*5, value=1:24)
    svt1 <- `[<-`(svt1, , 2:3, , value=0L)
    svt1 <- `[<-`(svt1, c(2, 4:6), 4, 2, value=1L)

    Nindex1 <- list(NULL, c(2:5, 8L), NULL)
    value1 <- -as.array(svt1)[ , 1:5, ]
    test_subassign_SVT_with_Rarray_or_SVT(svt0, Nindex1, value1)
    test_subassign_SVT_with_Rarray_or_SVT(svt1, Nindex1, value1)

    Nindex2 <- list(c(6:5, 1:2), 5L, NULL)
    value2 <- as.array(svt1)[1:4, 4, ]
    test_subassign_SVT_with_Rarray_or_SVT(svt0, Nindex2, value2)
    test_subassign_SVT_with_Rarray_or_SVT(svt1, Nindex2, value2)

    Nindex3 <- list(c(6:3, 2:5), 5:4, 2L)
    value3 <- matrix(101:116, ncol=2)
    test_subassign_SVT_with_Rarray_or_SVT(svt0, Nindex3, value3)
    test_subassign_SVT_with_Rarray_or_SVT(svt1, Nindex3, value3)
    test_subassign_SVT_with_Rarray_or_SVT(svt0, Nindex3, value3 + 0.5)
    test_subassign_SVT_with_Rarray_or_SVT(svt1, Nindex3, value3 + 0.5)
})

test_that("subassign an SVT_SparseArray object by an Nindex", {

    ## --- with an ordinary array on the right ---

    svt0 <- SVT_SparseArray(dim=c(4, 6), type="integer",
                            dimnames=list(letters[1:4], LETTERS[1:6]))
    m0 <- as.matrix(svt0)

    Rarray <- array(101:103, dim=c(1, 3))
    m <- `[<-`(m0, 2, 3:5, value=Rarray)
    svt <- `[<-`(svt0, 2, 3:5, value=Rarray)
    check_array_like_object(svt, "SVT_SparseMatrix", m)

    ## --- with an ordinary vector on the right that does not      ---
    ## --- get recycled along the **first** dimension of the array ---

    Rvector <- 201:202
    m <- `[<-`(m0, 4, 1:4, value=Rvector)
    svt <- `[<-`(svt0, 4, 1:4, value=Rvector)
    check_array_like_object(svt, "SVT_SparseMatrix", m)

    ## --- with a "short vector" on the right (gets recycled ---
    ## --- along the first dimension of the array)           ---

    set.seed(123)
    a0 <- array(0L, c(180, 400, 50))
    a0[sample(length(a0), 1e6)] <- sample(10L, 1e6, replace=TRUE)
    svt0 <- as(a0, "SVT_SparseArray")

    ## Wipe out all nonzeros:
    a <- `[<-`(a0, , , , value=0L)
    svt <- `[<-`(svt0, , , , value=0L)
    check_array_like_object(svt, "SVT_SparseArray", a)
    expect_null(svt@SVT)

    ## Wipe out all nonzeros in a column:
    a <- `[<-`(a0, , 8, 1, value=0L)
    svt <- `[<-`(svt0, , 8, 1, value=0L)
    check_array_like_object(svt, "SVT_SparseArray", a)
    expect_null(svt@SVT[[1L]][[8L]])
    i0 <- nzwhich(a0[ , 8, 1])
    svt2 <- `[<-`(svt0, i0, 8, 1, value=0L)
    expect_identical(svt2, svt)

    ## Wipe out all nonzeros in a row:
    a <- `[<-`(a0, 17, , 1, value=0L)
    svt <- `[<-`(svt0, 17, , 1, value=0L)
    check_array_like_object(svt, "SVT_SparseArray", a)
    j0 <- nzwhich(a0[17, , 1])
    svt2 <- `[<-`(svt0, 17, j0, 1, value=0L)
    expect_identical(svt2, svt)

    ## Inject zeros at random positions in a column:
    i <- sample(180L, 20L)
    a <- `[<-`(a0, i, 8, 1, value=0L)
    svt <- `[<-`(svt0, i, 8, 1, value=0L)
    check_array_like_object(svt, "SVT_SparseArray", a)

    ## Inject zeros at random positions in a row:
    j <- sample(400L, 50L)
    a <- `[<-`(a0, 17, j, 1, value=0L)
    svt <- `[<-`(svt0, 17, j, 1, value=0L)
    check_array_like_object(svt, "SVT_SparseArray", a)

    ## Inject zeros in a random set of rows:
    a <- `[<-`(a0, i, , 1, value=0L)
    svt <- `[<-`(svt0, i, , 1, value=0L)
    check_array_like_object(svt, "SVT_SparseArray", a)

    ## Inject zeros in a random set of columns:
    a <- `[<-`(a0, , j, 1, value=0L)
    svt <- `[<-`(svt0, , j, 1, value=0L)
    check_array_like_object(svt, "SVT_SparseArray", a)

    ## Inject zeros at random positions:
    a <- `[<-`(a0, i, j, 1, value=0L)
    svt <- `[<-`(svt0, i, j, 1, value=0L)
    check_array_like_object(svt, "SVT_SparseArray", a)
    a <- `[<-`(a0, i, j, , value=0L)
    svt <- `[<-`(svt0, i, j, , value=0L)
    check_array_like_object(svt, "SVT_SparseArray", a)

    ## Inject fixed nonzero value at random positions in a column:
    a <- `[<-`(a0, i, 8, 1, value=-555L)
    svt <- `[<-`(svt0, i, 8, 1, value=-555L)
    check_array_like_object(svt, "SVT_SparseArray", a)

    ## Inject fixed nonzero value at random positions in a row:
    a <- `[<-`(a0, 17, j, 1, value=-555L)
    svt <- `[<-`(svt0, 17, j, 1, value=-555L)
    check_array_like_object(svt, "SVT_SparseArray", a)

    ## Inject fixed nonzero value at random positions:
    a <- `[<-`(a0, i, j, 1, value=-555L)
    svt <- `[<-`(svt0, i, j, 1, value=-555L)
    check_array_like_object(svt, "SVT_SparseArray", a)
    a <- `[<-`(a0, i, j, , value=-555L)
    svt <- `[<-`(svt0, i, j, , value=-555L)
    check_array_like_object(svt, "SVT_SparseArray", a)

    ## Inject short vector with recycling:
    value <- c(-(101:104), 0L)
    a <- `[<-`(a0, i, j, 1, value=value)
    svt <- `[<-`(svt0, i, j, 1, value=value)
    check_array_like_object(svt, "SVT_SparseArray", a)
    a <- `[<-`(a0, i, , , value=value)
    svt <- `[<-`(svt0, i, , , value=value)
    check_array_like_object(svt, "SVT_SparseArray", a)
})

if (SparseArray:::SVT_VERSION != 0L) {

test_that("handling of lacunar leaves in SVT_SparseArray subassignment", {

    run_tests <- function(type) {

        check_svt123_leaves <- function(expected_leaf) {
            leaf_nzvals <- expected_leaf[[1L]]
            if (!is.null(leaf_nzvals)) {
                type(leaf_nzvals) <- type
                expected_leaf[[1L]] <- leaf_nzvals
            }
            expect_identical(svt1@SVT, expected_leaf)
            expect_identical(svt2@SVT[[2L]], expected_leaf)
            expect_identical(svt3@SVT[[2L]], expected_leaf)
            expect_identical(svt3@SVT[[4L]], expected_leaf)
            expect_identical(svt3@SVT[[5L]], expected_leaf)
        }

        svt1[2:3] <- -99L
        svt2[6:7] <- -99L
        svt3[2:3, c(2L, 4:5)] <- -99L
          m3[2:3, c(2L, 4:5)] <- -99L
        check_array_like_object(svt3, "SVT_SparseMatrix", m3)
        check_svt123_leaves(list(c(-99L, -99L), c(1L, 2L)))
        expect_identical(as(m3, "SVT_SparseMatrix"), svt3)

        svt1[2:3] <- 0L
        svt2[6:7] <- 0L
        svt3[2:3, c(2L, 4:5)] <- 0L
          m3[2:3, c(2L, 4:5)] <- 0L
        check_array_like_object(svt3, "SVT_SparseMatrix", m3)
        check_svt123_leaves(NULL)
        expect_identical(as(m3, "SVT_SparseMatrix"), svt3)

        svt1[1:4] <- 101:104
        svt2[5:8] <- 101:104
        svt3[1:4, c(2L, 4:5)] <- 101:104
          m3[1:4, c(2L, 4:5)] <- 101:104
        check_array_like_object(svt3, "SVT_SparseMatrix", m3)
        check_svt123_leaves(list(101:104, 0:3))
        expect_identical(as(m3, "SVT_SparseMatrix"), svt3)

        svt1[2:4] <- 1L
        svt2[6:8] <- 1L
        svt3[2:4, c(2L, 4:5)] <- 1L
          m3[2:4, c(2L, 4:5)] <- 1L
        check_array_like_object(svt3, "SVT_SparseMatrix", m3)
        check_svt123_leaves(list(c(101L, 1L, 1L, 1L), 0:3))
        expect_identical(as(m3, "SVT_SparseMatrix"), svt3)

        svt1[1:2] <- 0L
        svt2[5:6] <- 0L
        svt3[1:2, c(2L, 4:5)] <- 0L
          m3[1:2, c(2L, 4:5)] <- 0L
        check_array_like_object(svt3, "SVT_SparseMatrix", m3)
        check_svt123_leaves(make_lacunar_leaf(type, 2:3))
        expect_identical(as(m3, "SVT_SparseMatrix"), svt3)

        svt1[1L] <- NA
        svt2[5L] <- NA
        svt3[1L, c(2L, 4:5)] <- NA
          m3[1L, c(2L, 4:5)] <- NA
        check_array_like_object(svt3, "SVT_SparseMatrix", m3)
        check_svt123_leaves(list(c(NA, 1L, 1L), c(0L, 2L, 3L)))
        expect_identical(as(m3, "SVT_SparseMatrix"), svt3)

        svt1[1L] <- 1L
        svt2[5L] <- 1L
        svt3[1L, c(2L, 4:5)] <- 1L
          m3[1L, c(2L, 4:5)] <- 1L
        check_array_like_object(svt3, "SVT_SparseMatrix", m3)
        check_svt123_leaves(make_lacunar_leaf(type, c(0L, 2L, 3L)))
        expect_identical(as(m3, "SVT_SparseMatrix"), svt3)

        svt1[3L] <- 11L
        svt2[7L] <- 11L
        svt3[3L, c(2L, 4:5)] <- 11L
          m3[3L, c(2L, 4:5)] <- 11L
        check_array_like_object(svt3, "SVT_SparseMatrix", m3)
        check_svt123_leaves(list(c(1L, 11L, 1L), c(0L, 2L, 3L)))
        expect_identical(as(m3, "SVT_SparseMatrix"), svt3)

        svt1[2:3] <- 0L
        svt2[6:7] <- 0L
        svt3[2:3, c(2L, 4:5)] <- 0L
          m3[2:3, c(2L, 4:5)] <- 0L
        check_array_like_object(svt3, "SVT_SparseMatrix", m3)
        check_svt123_leaves(make_lacunar_leaf(type, c(0L, 3L)))
        expect_identical(as(m3, "SVT_SparseMatrix"), svt3)

        svt1[1:4] <- 1:0
        svt2[5:8] <- 1:0
        svt3[1:4, c(2L, 4:5)] <- 1:0
          m3[1:4, c(2L, 4:5)] <- 1:0
        check_array_like_object(svt3, "SVT_SparseMatrix", m3)
        check_svt123_leaves(make_lacunar_leaf(type, c(0L, 2L)))
        expect_identical(as(m3, "SVT_SparseMatrix"), svt3)

        svt1[c(1L, 4L)] <- 1L
        svt2[c(5L, 8L)] <- 1L
        svt3[c(1L, 4L), c(2L, 4:5)] <- 1L
          m3[c(1L, 4L), c(2L, 4:5)] <- 1L
        check_array_like_object(svt3, "SVT_SparseMatrix", m3)
        check_svt123_leaves(make_lacunar_leaf(type, c(0L, 2L, 3L)))
        expect_identical(as(m3, "SVT_SparseMatrix"), svt3)

        svt1[c(1L, 4:3)] <- c(0L, 1L, 0L)
        svt2[c(5L, 8:7)] <- c(0L, 1L, 0L)
        svt3[c(1L, 4:3), c(2L, 4:5)] <- c(0L, 1L, 0L)
          m3[c(1L, 4:3), c(2L, 4:5)] <- c(0L, 1L, 0L)
        check_array_like_object(svt3, "SVT_SparseMatrix", m3)
        check_svt123_leaves(make_lacunar_leaf(type, 3L))
        expect_identical(as(m3, "SVT_SparseMatrix"), svt3)

        svt1[1:2] <- 2:1
        svt2[5:6] <- 2:1
        svt3[1:2, c(2L, 4:5)] <- 2:1
          m3[1:2, c(2L, 4:5)] <- 2:1
        check_array_like_object(svt3, "SVT_SparseMatrix", m3)
        check_svt123_leaves(list(c(2L, 1L, 1L), c(0L, 1L, 3L)))
        expect_identical(as(m3, "SVT_SparseMatrix"), svt3)

        svt1[c(3L, 1L)] <- 1L
        svt2[c(7L, 5L)] <- 1L
        svt3[c(3L, 1L), c(2L, 4:5)] <- 1L
          m3[c(3L, 1L), c(2L, 4:5)] <- 1L
        check_array_like_object(svt3, "SVT_SparseMatrix", m3)
        check_svt123_leaves(make_lacunar_leaf(type, 0:3))
        expect_identical(as(m3, "SVT_SparseMatrix"), svt3)
    }

    svt1 <- as(array(0L, dim=4L), "SVT_SparseArray")  # 1D
    m3 <- matrix(0L, nrow=4, ncol=5)
    svt2 <- svt3 <- as(m3, "SVT_SparseMatrix")        # 2D
    run_tests("integer")

    type(svt1) <- "double"
    m3[ , 3L] <- 0.1
    svt2 <- svt3 <- as(m3, "SVT_SparseMatrix")
    run_tests("double")

})

}  # ----- end if (SparseArray:::SVT_VERSION != 0L) -----
