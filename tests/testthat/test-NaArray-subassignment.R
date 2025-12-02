.test_NaArray_subassignment_by_Mindex_and_Lindex <-
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

test_that("subassign an NaArray object by an Mindex or Lindex", {
    ## Only NAs.
    a0 <- array(NA_integer_, c(7, 10, 3),
                dimnames=list(NULL, letters[1:10], LETTERS[1:3]))
    Mindex3 <- rbind(c(7,  9, 3), c(7, 10, 3), c(6, 4, 3), c(2, 4, 3),
                     c(1, 10, 3), c(7, 10, 3), c(1, 1, 3), c(5, 4, 3),
                     c(2,  4, 3))
    vals <- c(0L, 12:18, NA)
    .test_NaArray_subassignment_by_Mindex_and_Lindex(a0, Mindex3, vals,
                                                     "NaArray")
    m0 <- a0[ , , 1]  # 2D
    Mindex2 <- Mindex3[ , -3]
    .test_NaArray_subassignment_by_Mindex_and_Lindex(m0, Mindex2, vals,
                                                     "NaMatrix")
    x0 <- as.array(m0[1, ])  # 1D
    Mindex1 <- Mindex2[ , -2, drop=FALSE]
    .test_NaArray_subassignment_by_Mindex_and_Lindex(x0, Mindex1, vals,
                                                     "NaArray")

    ## Add some non-NA elements.
    a0 <- make_3D_double_array(NA_real_)
    Mindex23 <- rbind(cbind(Mindex2, 1L), Mindex3)
    vals2 <- c(vals, vals)
    Mindex0 <- nnawhich(a0, arr.ind=TRUE)
    .test_NaArray_subassignment_by_Mindex_and_Lindex(a0, Mindex23, vals2,
                                                     "NaArray")
    .test_NaArray_subassignment_by_Mindex_and_Lindex(a0, Mindex0, 0,
                                                     "NaArray")
    m0 <- a0[ , , 1]  # 2D
    .test_NaArray_subassignment_by_Mindex_and_Lindex(m0, Mindex2, vals,
                                                     "NaMatrix")
    x0 <- as.array(m0[1, ])  # 1D
    .test_NaArray_subassignment_by_Mindex_and_Lindex(x0, Mindex1, vals,
                                                     "NaArray")

    ## Integer array.
    a0 <- make_3D_double_array(NA_real_)
    suppressWarnings(storage.mode(a0) <- "integer")
    .test_NaArray_subassignment_by_Mindex_and_Lindex(a0, Mindex23, vals2,
                                                     "NaArray")
    .test_NaArray_subassignment_by_Mindex_and_Lindex(a0, Mindex0, 0L,
                                                     "NaArray")

    ## Array type changed by subassignment.
    a0 <- make_3D_double_array(NA_real_)
    vals2 <- complex(real=vals2, imaginary=-0.75)
    .test_NaArray_subassignment_by_Mindex_and_Lindex(a0, Mindex23, vals2,
                                                     "NaArray")
    .test_NaArray_subassignment_by_Mindex_and_Lindex(a0, Mindex0, -9.99,
                                                     "NaArray")

    ## Assign random values to random array locations.
    set.seed(123)
    Mindex <- Lindex2Mindex(sample(length(a0)), dim(a0))
    vals <- sample(c(0:5, NA), length(a0), replace=TRUE)
    .test_NaArray_subassignment_by_Mindex_and_Lindex(a0, Mindex, vals,
                                                     "NaArray")
    Mindex <- Lindex2Mindex(sample(length(a0), 5000, replace=TRUE), dim(a0))
    vals <- sample(c(-99:99, NA), 5000, replace=TRUE)
    .test_NaArray_subassignment_by_Mindex_and_Lindex(a0, Mindex, vals,
                                                     "NaArray")
})

test_that("SparseArray:::.subassign_NaSVT_with_short_Rvector()", {
    subassign_NaSVT_with_short_Rvector <-
        SparseArray:::.subassign_NaSVT_with_short_Rvector

    test_subassign_NaSVT_with_short_Rvector <-
        function(naa0, Nindex, value, expected_type=type(value)) {
            naa <- subassign_NaSVT_with_short_Rvector(naa0, Nindex, value)
            expect_identical(type(naa), expected_type)
            a0 <- as.array(`type<-`(naa0, expected_type))
            a <- S4Arrays:::subassign_by_Nindex(a0, Nindex, value)
            expected_class <- if (is.matrix(a)) "NaMatrix" else "NaArray"
            check_array_like_object(naa, expected_class, a)
        }

    ## --- 1D objects ---

    naa0 <- NaArray(dim=10, type="integer", dimnames=list(LETTERS[1:10]))
    x0 <- as.array(naa0)
    Nindex1 <- list(c(6:9, 2L))
    Nindex2 <- list(NULL)
    Nindex3 <- list(c(10L, 3:5, 3L))

    naa1 <- subassign_NaSVT_with_short_Rvector(naa0, Nindex1, c(0:3, NA))
    a1 <- S4Arrays:::subassign_by_Nindex(x0, Nindex1, c(0:3, NA))
    check_array_like_object(naa1, "NaArray", a1)

    for (Nindex in list(Nindex1, Nindex2, Nindex3)) {
        value <- c(TRUE, NA, TRUE, NA, FALSE)
        test_subassign_NaSVT_with_short_Rvector(naa0, Nindex, value, "integer")
        test_subassign_NaSVT_with_short_Rvector(naa1, Nindex, value, "integer")
        value <- c(-2:1, NA)
        test_subassign_NaSVT_with_short_Rvector(naa0, Nindex, value)
        test_subassign_NaSVT_with_short_Rvector(naa1, Nindex, value)
        value <- c(-pi, NaN, 0, -Inf, NA)
        test_subassign_NaSVT_with_short_Rvector(naa0, Nindex, value)
        test_subassign_NaSVT_with_short_Rvector(naa1, Nindex, value)
        value <- 2.44 - value * 8i
        test_subassign_NaSVT_with_short_Rvector(naa0, Nindex, value)
        test_subassign_NaSVT_with_short_Rvector(naa1, Nindex, value)
        value <- c("hello", "", "world", NA, "!")
        test_subassign_NaSVT_with_short_Rvector(naa0, Nindex, value)
        test_subassign_NaSVT_with_short_Rvector(naa1, Nindex, value)
    }

    ## --- 2D objects ---

    naa0 <- NaArray(dim=c(10, 7), type="integer",
                    dimnames=list(LETTERS[1:10], letters[1:7]))
    m0 <- as.matrix(naa0)
    Nindex1 <- list(c(6:9, 2L), NULL)
    Nindex2 <- list(NULL, 2:4)
    Nindex3 <- list(c(10L, 3:5, 3L), c(6:3, 5L, 1L))
    Nindex4 <- list(c(3L, 5:2), c(1L, 7L))

    naa1 <- subassign_NaSVT_with_short_Rvector(naa0, Nindex1, c(0:3, NA))
    m1 <- S4Arrays:::subassign_by_Nindex(m0, Nindex1, c(0:3, NA))
    check_array_like_object(naa1, "NaMatrix", m1)

    for (Nindex in list(Nindex1, Nindex2, Nindex3, Nindex4)) {
        value <- c(TRUE, NA, TRUE, NA, FALSE)
        test_subassign_NaSVT_with_short_Rvector(naa0, Nindex, value, "integer")
        test_subassign_NaSVT_with_short_Rvector(naa1, Nindex, value, "integer")
        value <- c(-2:1, NA)
        test_subassign_NaSVT_with_short_Rvector(naa0, Nindex, value)
        test_subassign_NaSVT_with_short_Rvector(naa1, Nindex, value)
        value <- c(-pi, NaN, 0, -Inf, NA)
        test_subassign_NaSVT_with_short_Rvector(naa0, Nindex, value)
        test_subassign_NaSVT_with_short_Rvector(naa1, Nindex, value)
        value <- 2.44 - value * 8i
        test_subassign_NaSVT_with_short_Rvector(naa0, Nindex, value)
        test_subassign_NaSVT_with_short_Rvector(naa1, Nindex, value)
        value <- c("hello", "", "world", NA, "!")
        test_subassign_NaSVT_with_short_Rvector(naa0, Nindex, value)
        test_subassign_NaSVT_with_short_Rvector(naa1, Nindex, value)
    }

    ## --- 3D objects ---

    naa0 <- NaArray(dim=c(10, 2, 7), type="integer",
                    dimnames=list(LETTERS[1:10], NULL, letters[1:7]))
    a0 <- as.array(naa0)
    Nindex1 <- list(c(6:9, 2L), NULL, NULL)
    Nindex2 <- list(NULL, 2L, 2:4)
    Nindex3 <- list(c(10L, 3:5, 3L), 2L, c(6:3, 5L, 1L))
    Nindex4 <- list(c(3L, 5:2), 2:1, c(1L, 7L))
    Nindex5 <- list(NULL, NULL, c(7L, 7L))

    naa1 <- subassign_NaSVT_with_short_Rvector(naa0, Nindex1, c(0:3, NA))
    a1 <- S4Arrays:::subassign_by_Nindex(a0, Nindex1, c(0:3, NA))
    check_array_like_object(naa1, "NaArray", a1)

    for (Nindex in list(Nindex1, Nindex2, Nindex3, Nindex4, Nindex5)) {
        value <- c(TRUE, NA, TRUE, NA, FALSE)
        test_subassign_NaSVT_with_short_Rvector(naa0, Nindex, value, "integer")
        test_subassign_NaSVT_with_short_Rvector(naa1, Nindex, value, "integer")
        value <- c(-2:1, NA)
        test_subassign_NaSVT_with_short_Rvector(naa0, Nindex, value)
        test_subassign_NaSVT_with_short_Rvector(naa1, Nindex, value)
        value <- c(-pi, NaN, 0, -Inf, NA)
        test_subassign_NaSVT_with_short_Rvector(naa0, Nindex, value)
        test_subassign_NaSVT_with_short_Rvector(naa1, Nindex, value)
        value <- 2.44 - value * 8i
        test_subassign_NaSVT_with_short_Rvector(naa0, Nindex, value)
        test_subassign_NaSVT_with_short_Rvector(naa1, Nindex, value)
        value <- c("hello", "", "world", NA, "!")
        test_subassign_NaSVT_with_short_Rvector(naa0, Nindex, value)
        test_subassign_NaSVT_with_short_Rvector(naa1, Nindex, value)
    }
})

test_that(".subassign_NaSVT_with_Rarray() and .subassign_NaSVT_with_NaSVT()", {
    test_subassign_NaSVT_with_Rarray_or_NaSVT <-
        function(naa0, Nindex, value, expected_type=type(value)) {
            naa <- SparseArray:::.subassign_NaSVT_with_Rarray(naa0, Nindex,
                                                              value)
            expect_identical(type(naa), expected_type)
            a0 <- as.array(`type<-`(naa0, expected_type))
            a <- S4Arrays:::subassign_by_Nindex(a0, Nindex, value)
            expected_class <- if (is.matrix(a)) "NaMatrix" else "NaArray"
            check_array_like_object(naa, expected_class, a)
            value <- NaArray(value)
            naa2 <- SparseArray:::.subassign_NaSVT_with_NaSVT(naa0, Nindex,
                                                              value)
            expect_identical(naa2, naa)
        }

    ## --- 1D objects ---

    naa0 <- NaArray(dim=10, type="integer")
    naa1 <- NaArray(array(c(0L, 0L, 103:104, 0L, 0L, 107:109, 0L)))
    Nindex1 <- list(c(2:5, 8L))
    Nindex2 <- list(c(9L, 3:5, 3L))

    for (Nindex in list(Nindex1, Nindex2)) {
        value <- array(rep.int(NA_integer_, 5))
        test_subassign_NaSVT_with_Rarray_or_NaSVT(naa0, Nindex, value)
        test_subassign_NaSVT_with_Rarray_or_NaSVT(naa1, Nindex, value)
        value <- array(c(-2:1, NA))
        test_subassign_NaSVT_with_Rarray_or_NaSVT(naa0, Nindex, value)
        test_subassign_NaSVT_with_Rarray_or_NaSVT(naa1, Nindex, value)
        value <- array(c(-pi, NaN, 0, -Inf, NA))
        test_subassign_NaSVT_with_Rarray_or_NaSVT(naa0, Nindex, value)
        test_subassign_NaSVT_with_Rarray_or_NaSVT(naa1, Nindex, value)
        value <- array(2.44 - value * 8i)
        test_subassign_NaSVT_with_Rarray_or_NaSVT(naa0, Nindex, value)
        test_subassign_NaSVT_with_Rarray_or_NaSVT(naa1, Nindex, value)
        value <- array(c("hello", "", "world", NA, "!"))
        test_subassign_NaSVT_with_Rarray_or_NaSVT(naa0, Nindex, value)
        test_subassign_NaSVT_with_Rarray_or_NaSVT(naa1, Nindex, value)
    }

    value <- array(c(NA, -1:1, 0L, 2L, NA, 0L, 3:4))
    test_subassign_NaSVT_with_Rarray_or_NaSVT(naa0, list(NULL), value)
    test_subassign_NaSVT_with_Rarray_or_NaSVT(naa1, list(NULL), value)

    ## --- 2D objects ---

    naa0 <- NaArray(dim=c(10, 6), type="integer")
    naa1 <- `[<-`(naa0, (1:20)*3, value=1:20)
    naa1 <- `[<-`(naa1, , 2, value=NA_integer_)

    Nindex1 <- list(c(2:5, 8L), NULL)
    value1 <- -as.matrix(naa1)[6:10, ]
    test_subassign_NaSVT_with_Rarray_or_NaSVT(naa0, Nindex1, value1)
    test_subassign_NaSVT_with_Rarray_or_NaSVT(naa1, Nindex1, value1)

    Nindex2 <- list(NULL, 3:6)
    value2 <- -as.matrix(naa1)[ , 2:5]
    test_subassign_NaSVT_with_Rarray_or_NaSVT(naa0, Nindex2, value2)
    test_subassign_NaSVT_with_Rarray_or_NaSVT(naa1, Nindex2, value2)

    Nindex3 <- list(c(9L, 3:5, 3L), c(6L, 2L))
    value3 <- -as.matrix(naa1)[6:10, 5:6]
    test_subassign_NaSVT_with_Rarray_or_NaSVT(naa0, Nindex3, value3)
    test_subassign_NaSVT_with_Rarray_or_NaSVT(naa1, Nindex3, value3)

    ## --- 3D objects ---

    naa0 <- NaArray(dim=c(6, 10, 2), type="integer")
    naa1 <- `[<-`(naa0, (1:24)*5, value=1:24)
    naa1 <- `[<-`(naa1, , 2:3, , value=NA_integer_)
    naa1 <- `[<-`(naa1, c(2, 4:6), 4, 2, value=1L)

    Nindex1 <- list(NULL, c(2:5, 8L), NULL)
    value1 <- -as.array(naa1)[ , 1:5, ]
    test_subassign_NaSVT_with_Rarray_or_NaSVT(naa0, Nindex1, value1)
    test_subassign_NaSVT_with_Rarray_or_NaSVT(naa1, Nindex1, value1)

    Nindex2 <- list(c(6:5, 1:2), 5L, NULL)
    value2 <- as.array(naa1)[1:4, 4, ]
    test_subassign_NaSVT_with_Rarray_or_NaSVT(naa0, Nindex2, value2)
    test_subassign_NaSVT_with_Rarray_or_NaSVT(naa1, Nindex2, value2)

    Nindex3 <- list(c(6:3, 2:5), 5:4, 2L)
    value3 <- matrix(101:116, ncol=2)
    test_subassign_NaSVT_with_Rarray_or_NaSVT(naa0, Nindex3, value3)
    test_subassign_NaSVT_with_Rarray_or_NaSVT(naa1, Nindex3, value3)
    test_subassign_NaSVT_with_Rarray_or_NaSVT(naa0, Nindex3, value3 + 0.5)
    test_subassign_NaSVT_with_Rarray_or_NaSVT(naa1, Nindex3, value3 + 0.5)
})

test_that("subassign an NaArray object by an Nindex", {

    ## --- with an ordinary array on the right ---

    naa0 <- NaArray(dim=c(4, 6), type="integer",
                    dimnames=list(letters[1:4], LETTERS[1:6]))
    m0 <- as.matrix(naa0)

    Rarray <- array(101:103, dim=c(1, 3))
    m <- `[<-`(m0, 2, 3:5, value=Rarray)
    naa <- `[<-`(naa0, 2, 3:5, value=Rarray)
    check_array_like_object(naa, "NaMatrix", m)

    ## --- with an ordinary vector on the right that does not      ---
    ## --- get recycled along the **first** dimension of the array ---

    Rvector <- 201:202
    m <- `[<-`(m0, 4, 1:4, value=Rvector)
    naa <- `[<-`(naa0, 4, 1:4, value=Rvector)
    check_array_like_object(naa, "NaMatrix", m)

    ## --- with a "short vector" on the right (gets recycled ---
    ## --- along the first dimension of the array)           ---

    set.seed(123)
    a0 <- array(NA_integer_, c(180, 400, 50))
    a0[sample(length(a0), 1e6)] <- sample(10L, 1e6, replace=TRUE)
    naa0 <- as(a0, "NaArray")

    ## Wipe out all non-NAs:
    a <- `[<-`(a0, , , , value=NA_integer_)
    naa <- `[<-`(naa0, , , , value=NA_integer_)
    check_array_like_object(naa, "NaArray", a)
    expect_null(naa@NaSVT)

    ## Wipe out all non-NAs in a column:
    a <- `[<-`(a0, , 8, 1, value=NA_integer_)
    naa <- `[<-`(naa0, , 8, 1, value=NA_integer_)
    check_array_like_object(naa, "NaArray", a)
    expect_null(naa@NaSVT[[1L]][[8L]])
    i0 <- nzwhich(a0[ , 8, 1])
    naa2 <- `[<-`(naa0, i0, 8, 1, value=NA_integer_)
    expect_identical(naa2, naa)

    ## Wipe out all non-NAs in a row:
    a <- `[<-`(a0, 17, , 1, value=NA_integer_)
    naa <- `[<-`(naa0, 17, , 1, value=NA_integer_)
    check_array_like_object(naa, "NaArray", a)
    j0 <- nzwhich(a0[17, , 1])
    naa2 <- `[<-`(naa0, 17, j0, 1, value=NA_integer_)
    expect_identical(naa2, naa)

    ## Inject NAs at random positions in a column:
    i <- sample(180L, 20L)
    a <- `[<-`(a0, i, 8, 1, value=NA_integer_)
    naa <- `[<-`(naa0, i, 8, 1, value=NA_integer_)
    check_array_like_object(naa, "NaArray", a)

    ## Inject NAs at random positions in a row:
    j <- sample(400L, 50L)
    a <- `[<-`(a0, 17, j, 1, value=NA_integer_)
    naa <- `[<-`(naa0, 17, j, 1, value=NA_integer_)
    check_array_like_object(naa, "NaArray", a)

    ## Inject NAs in a random set of rows:
    a <- `[<-`(a0, i, , 1, value=NA_integer_)
    naa <- `[<-`(naa0, i, , 1, value=NA_integer_)
    check_array_like_object(naa, "NaArray", a)

    ## Inject NAs in a random set of columns:
    a <- `[<-`(a0, , j, 1, value=NA_integer_)
    naa <- `[<-`(naa0, , j, 1, value=NA_integer_)
    check_array_like_object(naa, "NaArray", a)

    ## Inject NAs at random positions:
    a <- `[<-`(a0, i, j, 1, value=NA_integer_)
    naa <- `[<-`(naa0, i, j, 1, value=NA_integer_)
    check_array_like_object(naa, "NaArray", a)
    a <- `[<-`(a0, i, j, , value=0L)
    naa <- `[<-`(naa0, i, j, , value=0L)
    check_array_like_object(naa, "NaArray", a)

    ## Inject fixed non-NA value at random positions in a column:
    a <- `[<-`(a0, i, 8, 1, value=-555L)
    naa <- `[<-`(naa0, i, 8, 1, value=-555L)
    check_array_like_object(naa, "NaArray", a)

    ## Inject fixed non-NA value at random positions in a row:
    a <- `[<-`(a0, 17, j, 1, value=-555L)
    naa <- `[<-`(naa0, 17, j, 1, value=-555L)
    check_array_like_object(naa, "NaArray", a)

    ## Inject fixed non-NA value at random positions:
    a <- `[<-`(a0, i, j, 1, value=-555L)
    naa <- `[<-`(naa0, i, j, 1, value=-555L)
    check_array_like_object(naa, "NaArray", a)
    a <- `[<-`(a0, i, j, , value=-555L)
    naa <- `[<-`(naa0, i, j, , value=-555L)
    check_array_like_object(naa, "NaArray", a)

    ## Inject short vector with recycling:
    value <- c(-(101:104), NA_integer_)
    a <- `[<-`(a0, i, j, 1, value=value)
    naa <- `[<-`(naa0, i, j, 1, value=value)
    check_array_like_object(naa, "NaArray", a)
    a <- `[<-`(a0, i, , , value=value)
    naa <- `[<-`(naa0, i, , , value=value)
    check_array_like_object(naa, "NaArray", a)
})

