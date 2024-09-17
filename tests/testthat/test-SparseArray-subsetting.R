### Subsetting a 1D ordinary array always preserves the "dim" and "dimnames"
### attributes i.e. it always returns another 1D array. This is inconsistent
### with the multidimensional case where, for example linear, subsetting (i.e.
### subsetting by a numeric vector or matrix) returns an ordinary vector
### (atomic or list). At the root of the problem is the behavior of
### base::drop() on a 1D array (see drop_even_if_1D() function in S4Arrays,
### in file R/dim-tuning-utils.R). So we "fix" subsetting of a 1D ordinary
### array by passing the result of the subsetting operation thru
### S4Arrays:::drop_even_if_1D().
.test_SVT_subset_by_Lindex_and_Mindex <- function(a, Mindex0)
{
    stopifnot(is.array(a))
    svt <- as(a, "SVT_SparseArray")
    Lindex0 <- Mindex2Lindex(Mindex0, dim(a))

    expected <- S4Arrays:::drop_even_if_1D(a[Mindex0])
    expect_identical(svt[Lindex0], expected)
    expect_identical(svt[Mindex0], expected)
    expect_identical(svt[Lindex0 + 0.99], expected)
    expect_identical(svt[Mindex0 + 0.99], expected)

    revLindex <- rev(Lindex0)
    revMindex <- Lindex2Mindex(revLindex, dim(a))
    expected <- S4Arrays:::drop_even_if_1D(a[revLindex])
    expect_identical(svt[revLindex], expected)
    expect_identical(svt[revMindex], expected)
    type(revLindex) <- "double"
    type(revMindex) <- "double"
    expect_identical(svt[revLindex], expected)
    expect_identical(svt[revMindex], expected)
    expect_identical(svt[revLindex + 0.99], expected)
    expect_identical(svt[revMindex + 0.99], expected)

    Lindex <- c(NA, Lindex0, NA, NA, Lindex0, NA, rev(Lindex0))
    ## Lindex2Mindex() and 'svt[Mindex]' don't accept NAs at the moment.
    #Mindex <- Lindex2Mindex(Lindex, dim(a))
    expected <- S4Arrays:::drop_even_if_1D(a[Lindex])
    expect_identical(svt[Lindex], expected)
    #expect_identical(svt[Mindex], expected)
    Lindex[1L] <- NaN      # this coerces 'Lindex' to "numeric"
    #Mindex[1L, 1L] <- NaN  # this coerces 'Mindex' to "numeric"
    expect_identical(svt[Lindex], expected)
    #expect_identical(svt[Mindex], expected)
}

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
    .test_SVT_subset_by_Lindex_and_Mindex(a0, Mindex3)

    ## --- 2D ---
    m0 <- a0[ , , 1]
    .test_SVT_subset_by_Lindex_and_Mindex(m0, Mindex2)

    ## --- 1D ---
    x0 <- as.array(m0[1, ])
    Mindex1 <- Mindex2[ , -2, drop=FALSE]
    .test_SVT_subset_by_Lindex_and_Mindex(x0, Mindex1)
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
    check_array_like_object(svt, "SVT_SparseArray", a)
    expect_identical(as(a, "SVT_SparseArray"), svt)
    m   <- a0  [ , c(4:3, 8), 1, drop=TRUE]
    svt <- svt0[ , c(4:3, 8), 1, drop=TRUE]
    check_array_like_object(svt, "SVT_SparseMatrix", m)
    expect_identical(as(m, "SVT_SparseMatrix"), svt)

    a   <- a0  [7, , , drop=FALSE]
    svt <- svt0[7, , , drop=FALSE]
    check_array_like_object(svt, "SVT_SparseArray", a)
    expect_identical(as(a, "SVT_SparseArray"), svt)
    m   <- a0  [7, , , drop=TRUE]
    svt <- svt0[7, , , drop=TRUE]
    check_array_like_object(svt, "SVT_SparseMatrix", m)
    expect_identical(as(m, "SVT_SparseMatrix"), svt)

    ## --- 2D ---
    m0 <- a0[ , , 1]
    svt0 <- as(m0, "SVT_SparseMatrix")

    expect_identical(svt0[ , ], svt0)

    m   <- m0  [-5 , c(4:3, 8)]
    svt <- svt0[-5 , c(4:3, 8)]
    check_array_like_object(svt, "SVT_SparseMatrix", m)
    expect_identical(as(m, "SVT_SparseMatrix"), svt)

    expect_identical(svt0[ , 4], m0[ , 4])
    m   <- m0  [ , 4, drop=FALSE]
    svt <- svt0[ , 4, drop=FALSE]
    check_array_like_object(svt, "SVT_SparseMatrix", m)
    expect_identical(as(m, "SVT_SparseMatrix"), svt)

    expect_identical(svt0[6 , ], m0[6 , ])
    m   <- m0  [6, -1, drop=FALSE]
    svt <- svt0[6, -1, drop=FALSE]
    check_array_like_object(svt, "SVT_SparseMatrix", m)
    expect_identical(as(m, "SVT_SparseMatrix"), svt)

    ## --- 1D ---
    x0 <- as.array(m0[6, ])
    svt0 <- as(x0, "SVT_SparseArray")

    expect_identical(svt0[ ], svt0)

    x   <- x0  [c(8:4, 1, 4), drop=FALSE]
    svt <- svt0[c(8:4, 1, 4), drop=FALSE]
    check_array_like_object(svt, "SVT_SparseArray", x)
    expect_identical(as(x, "SVT_SparseArray"), svt)
    x   <- x0  [-4, drop=FALSE]
    svt <- svt0[-4, drop=FALSE]
    check_array_like_object(svt, "SVT_SparseArray", x)
    expect_identical(as(x, "SVT_SparseArray"), svt)
    x   <- x0  [-4, drop=FALSE]
    svt <- svt0[-4, drop=FALSE]
    check_array_like_object(svt, "SVT_SparseArray", x)
    expect_identical(as(x, "SVT_SparseArray"), svt)

    subscript <- c("d", "j", "j", "h")
    x   <- x0  [subscript, drop=FALSE]  # 'drop=TRUE' would do the same thing!
    svt <- svt0[subscript, drop=FALSE]
    check_array_like_object(svt, "SVT_SparseArray", x)
    expect_identical(as(x, "SVT_SparseArray"), svt)
    expect_identical(svt0[subscript], S4Arrays:::drop_even_if_1D(x))
})

