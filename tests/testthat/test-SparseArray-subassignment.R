.test_SparseArray_subassignment_by_Mindex_and_Lindex <-
    function(a0, Mindex, vals, expected_class)
{
    object0 <- as(a0, expected_class)
    Lindex <- Mindex2Lindex(Mindex, dim(a0))

    a <- `[<-`(a0, Mindex, value=vals)
    object <- `[<-`(object0, Mindex, value=vals)
    check_SparseArray_object(object, expected_class, a)
    object <- `[<-`(object0, Lindex, value=vals)
    check_SparseArray_object(object, expected_class, a)
    object <- `[<-`(object0, as.double(Lindex), value=vals)
    check_SparseArray_object(object, expected_class, a)
    object <- `[<-`(object0, Lindex + 0.5, value=vals)
    check_SparseArray_object(object, expected_class, a)
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

    ## Add some nonzero values.
    a0 <- make_3D_test_array()
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
    a0 <- make_3D_test_array()
    suppressWarnings(storage.mode(a0) <- "integer")
    .test_SparseArray_subassignment_by_Mindex_and_Lindex(a0, Mindex23, vals2,
                                                         "SVT_SparseArray")
    .test_SparseArray_subassignment_by_Mindex_and_Lindex(a0, Mindex0, 0L,
                                                         "SVT_SparseArray")

    ## Array type changed by subassignment.
    a0 <- make_3D_test_array()
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

.test_SparseArray_subassignment_by_Nindex <-
    function(a0, index, vals, expected_class)
{
    object0 <- as(a0, expected_class)

    a <- `[<-`(a0, index, value=vals)
    object <- `[<-`(object0, index, value=vals)
    check_SparseArray_object(object, expected_class, a)
}

test_that(paste("subassign an SVT_SparseArray object by an Nindex",
                "and with a short vector"), {
    set.seed(123)
    a0 <- array(0L, c(180, 400, 50))
    a0[sample(length(a0), 1e6)] <- sample(10L, 1e6, replace=TRUE)
    svt0 <- as(a0, "SVT_SparseArray")

    ## Wipe out all nonzeros:
    a <- `[<-`(a0, , , , value=0L)
    svt <- `[<-`(svt0, , , , value=0L)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    expect_null(svt@SVT)

    ## Wipe out all nonzeros in a column:
    a <- `[<-`(a0, , 8, 1, value=0L)
    svt <- `[<-`(svt0, , 8, 1, value=0L)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    expect_null(svt@SVT[[1L]][[8L]])
    i0 <- nzwhich(a0[ , 8, 1])
    svt2 <- `[<-`(svt0, i0, 8, 1, value=0L)
    expect_identical(svt2, svt)

    ## Wipe out all nonzeros in a row:
    a <- `[<-`(a0, 17, , 1, value=0L)
    svt <- `[<-`(svt0, 17, , 1, value=0L)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    j0 <- nzwhich(a0[17, , 1])
    svt2 <- `[<-`(svt0, 17, j0, 1, value=0L)
    expect_identical(svt2, svt)

    ## Inject zeros at random positions in a column:
    i <- sample(180L, 20L)
    a <- `[<-`(a0, i, 8, 1, value=0L)
    svt <- `[<-`(svt0, i, 8, 1, value=0L)
    check_SparseArray_object(svt, "SVT_SparseArray", a)

    ## Inject zeros at random positions in a row:
    j <- sample(400L, 50L)
    a <- `[<-`(a0, 17, j, 1, value=0L)
    svt <- `[<-`(svt0, 17, j, 1, value=0L)
    check_SparseArray_object(svt, "SVT_SparseArray", a)

    ## Inject zeros in a random set of rows:
    a <- `[<-`(a0, i, , 1, value=0L)
    svt <- `[<-`(svt0, i, , 1, value=0L)
    check_SparseArray_object(svt, "SVT_SparseArray", a)

    ## Inject zeros in a random set of columns:
    a <- `[<-`(a0, , j, 1, value=0L)
    svt <- `[<-`(svt0, , j, 1, value=0L)
    check_SparseArray_object(svt, "SVT_SparseArray", a)

    ## Inject zeros at random positions:
    a <- `[<-`(a0, i, j, 1, value=0L)
    svt <- `[<-`(svt0, i, j, 1, value=0L)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    a <- `[<-`(a0, i, j, , value=0L)
    svt <- `[<-`(svt0, i, j, , value=0L)
    check_SparseArray_object(svt, "SVT_SparseArray", a)

    ## Inject fixed nonzero at random positions in a column:
    a <- `[<-`(a0, i, 8, 1, value=-555L)
    svt <- `[<-`(svt0, i, 8, 1, value=-555L)
    check_SparseArray_object(svt, "SVT_SparseArray", a)

    ## Inject fixed nonzero at random positions in a row:
    a <- `[<-`(a0, 17, j, 1, value=-555L)
    svt <- `[<-`(svt0, 17, j, 1, value=-555L)
    check_SparseArray_object(svt, "SVT_SparseArray", a)

    ## Inject fixed nonzero val at random positions:
    a <- `[<-`(a0, i, j, 1, value=-555L)
    svt <- `[<-`(svt0, i, j, 1, value=-555L)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    a <- `[<-`(a0, i, j, , value=-555L)
    svt <- `[<-`(svt0, i, j, , value=-555L)
    check_SparseArray_object(svt, "SVT_SparseArray", a)

    ## Inject short vector with recycling:
    value <- c(-(101:104), 0L)
    a <- `[<-`(a0, i, j, 1, value=value)
    svt <- `[<-`(svt0, i, j, 1, value=value)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
    a <- `[<-`(a0, i, , , value=value)
    svt <- `[<-`(svt0, i, , , value=value)
    check_SparseArray_object(svt, "SVT_SparseArray", a)
})

