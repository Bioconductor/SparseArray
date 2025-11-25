
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

test_that(".subset_SVT_as_Rarray() and .subset_SVT_as_SVT()", {
    test_subset_SVT_as_Rarray_or_SVT <-
        function(svt0, Nindex_list, a0) {
            for (Nindex in Nindex_list) {
                a <- S4Arrays:::subset_by_Nindex(a0, Nindex, drop=FALSE)
                expect_identical(
                    SparseArray:::.subset_SVT_as_Rarray(svt0, Nindex), a)
                svt <- SparseArray:::.subset_SVT_as_SVT(svt0, Nindex)
                check_SVT_SparseArray_object(svt, a)
            }
        }

    ## --- 1D objects ---

    x0 <- array(as.raw(c(0L, 72:73, 0L, 75L, 0L, 0L, 78:80)),
                dimnames=list(LETTERS[1:10]))
    Nindex_list <- list(list(NULL),
                        list(c(6:9, 2L)),
                        list(c(10L, 3:5, 3L)),
                        list(integer(0)))

    svt0 <- SVT_SparseArray(x0)
    test_subset_SVT_as_Rarray_or_SVT(svt0, Nindex_list, x0)
    type(svt0) <- "integer"
    test_subset_SVT_as_Rarray_or_SVT(svt0, Nindex_list, as.array(svt0))
    type(svt0) <- "double"
    test_subset_SVT_as_Rarray_or_SVT(svt0, Nindex_list, as.array(svt0))
    type(svt0) <- "complex"
    test_subset_SVT_as_Rarray_or_SVT(svt0, Nindex_list, as.array(svt0))
    type(svt0) <- "character"
    test_subset_SVT_as_Rarray_or_SVT(svt0, Nindex_list, as.array(svt0))
    type(svt0) <- "list"
    test_subset_SVT_as_Rarray_or_SVT(svt0, Nindex_list, as.array(svt0))

    ## --- 2D objects ---

    m0 <- rbind(a=x0, b=as.raw(0), c=rev(x0))
    Nindex_list <- list(list(NULL, NULL),
                        list(2L, NULL),
                        list(NULL, c(10L, 3:5, 3L)),
                        list(integer(0), 6:9))

    svt0 <- SVT_SparseArray(m0)
    test_subset_SVT_as_Rarray_or_SVT(svt0, Nindex_list, m0)
    type(svt0) <- "integer"
    test_subset_SVT_as_Rarray_or_SVT(svt0, Nindex_list, as.array(svt0))
    type(svt0) <- "double"
    test_subset_SVT_as_Rarray_or_SVT(svt0, Nindex_list, as.array(svt0))
    type(svt0) <- "complex"
    test_subset_SVT_as_Rarray_or_SVT(svt0, Nindex_list, as.array(svt0))
    type(svt0) <- "character"
    test_subset_SVT_as_Rarray_or_SVT(svt0, Nindex_list, as.array(svt0))
    type(svt0) <- "list"
    test_subset_SVT_as_Rarray_or_SVT(svt0, Nindex_list, as.array(svt0))

    ## --- 3D objects ---

    a0 <- array(0L, c(7, 10, 3),
                dimnames=list(NULL, letters[1:10], LETTERS[1:3]))
    a0[ , 2, 1] <- a0[c(1:2, 6), 4, 1] <- 1L
    a0[ , 8, 1] <- 81:87
    a0[ , -1, 3] <- 308:370
    Nindex_list <- list(list(NULL, NULL, NULL),
                        list(NULL, c(4:3, 8L), 1L),
                        list(7L, NULL, NULL),
                        list(NULL, NULL, 1L),
                        list(1:5, c(4:3, 8L, 1:3), c(1L, 3:2)))

    svt0 <- SVT_SparseArray(a0)
    test_subset_SVT_as_Rarray_or_SVT(svt0, Nindex_list, a0)
    type(svt0) <- "double"
    test_subset_SVT_as_Rarray_or_SVT(svt0, Nindex_list, as.array(svt0))
    type(svt0) <- "complex"
    test_subset_SVT_as_Rarray_or_SVT(svt0, Nindex_list, as.array(svt0))
    # Can't test for this at the moment because coercing a lacunar leaf
    # to "character" or "list" is not supported yet!
    #type(svt0) <- "character"
    #test_subset_SVT_as_Rarray_or_SVT(svt0, Nindex_list, as.array(svt0))
    #type(svt0) <- "list"
    #test_subset_SVT_as_Rarray_or_SVT(svt0, Nindex_list, as.array(svt0))
})

