
test_that(".tune_SVT_SparseArray_dims()", {
    .tune_SVT_SparseArray_dims <- SparseArray:::.tune_SVT_SparseArray_dims

    a0 <- array(0L, dim=c(5, 4, 1, 3))
    dimnames(a0) <- list(letters[1:5], NULL, NULL, LETTERS[1:3])
    a0[c(1:2, 8, 10, 15:17, 20, 24, 40)] <- (1:10)*10L
    svt0 <- SparseArray(a0)

    dim_tuner <- c(0L, 0L, 0L, 0L)
    svt <- .tune_SVT_SparseArray_dims(svt0, dim_tuner)
    expect_identical(svt, svt0)

    dim_tuner <- c(0L, 0L, -1L, 0L)
    svt <- .tune_SVT_SparseArray_dims(svt0, dim_tuner)
    check_array_like_object(svt, "SVT_SparseArray", drop(a0))
    expect_identical(as(drop(a0), "SVT_SparseArray"), svt)
    svt2 <- .tune_SVT_SparseArray_dims(svt, -dim_tuner)
    expect_identical(svt2, svt0)
  
    dim_tuner <- c(0L, 0L, 1L, 0L, 1L, 0L)
    svt <- .tune_SVT_SparseArray_dims(svt0, dim_tuner)
    a <- `dim<-`(a0, c(5, 4, 1, 1, 1, 3))
    dimnames(a)[c(1, 2, 4, 6)] <- dimnames(a0)
    check_array_like_object(svt, "SVT_SparseArray", a)
    expect_identical(as(a, "SVT_SparseArray"), svt)
    svt2 <- .tune_SVT_SparseArray_dims(svt, -dim_tuner)
    expect_identical(svt2, svt0)

    dim_tuner <- c(1L, 1L, 0L, 1L, 0L, -1L, 0L, 1L)
    svt <- .tune_SVT_SparseArray_dims(svt0, dim_tuner)
    a <- `dim<-`(a0, c(1, 1, 5, 1, 4, 3, 1))
    dimnames(a)[c(3, 5, 6)] <- dimnames(a0)[c(1, 2, 4)]
    check_array_like_object(svt, "SVT_SparseArray", a)
    expect_identical(as(a, "SVT_SparseArray"), svt)
    svt2 <- .tune_SVT_SparseArray_dims(svt, -dim_tuner)
    expect_identical(svt2, svt0)
})

test_that("`dim<-` and drop() on an SVT_SparseArray object", {
    ## --- integer matrix ---
    m0 <- matrix(c(1:0, -99L, NA, 2:1), ncol=3,
                 dimnames=list(LETTERS[1:2], letters[1:3]))
    svt0 <- SparseArray(m0)

    a <- `dim<-`(m0, c(1L, dim(m0)))
    dimnames(a) <- c(list(NULL), dimnames(m0))
    svt <- `dim<-`(svt0, c(1L, dim(svt0)))
    check_array_like_object(svt, "SVT_SparseArray", a)
    expect_identical(as(a, "SVT_SparseArray"), svt)

    expect_identical(drop(svt), svt0)

    m1 <- m0[2, , drop=FALSE]
    expect_identical(drop(SparseArray(m1)), drop(m1))
    x1 <- as.array(drop(m1))  # 1D array
    ## base::drop() is kind of messed up in the 1D case (see drop_even_if_1D()
    ## function in S4Arrays, in file R/dim-tuning-utils.R). But drop() does
    ## the right thing on a 1D SparseArray.
    expect_identical(drop(SparseArray(x1)), S4Arrays:::drop_even_if_1D(x1))

    ## --- double matrix ---
    m0[1, 3] <- NaN
    svt0 <- SparseArray(m0)

    a <- `dim<-`(m0, c(1L, dim(m0)))
    dimnames(a) <- c(list(NULL), dimnames(m0))
    svt <- `dim<-`(svt0, c(1L, dim(svt0)))
    check_array_like_object(svt, "SVT_SparseArray", a)
    expect_identical(as(a, "SVT_SparseArray"), svt)

    expect_identical(drop(svt), svt0)

    m1 <- m0[2, , drop=FALSE]
    expect_identical(drop(SparseArray(m1)), drop(m1))
    x1 <- as.array(drop(m1))  # 1D array
    ## base::drop() is kind of messed up in the 1D case (see drop_even_if_1D()
    ## function in S4Arrays, in file R/dim-tuning-utils.R). But drop() does
    ## the right thing on a 1D SparseArray.
    expect_identical(drop(SparseArray(x1)), S4Arrays:::drop_even_if_1D(x1))
})

