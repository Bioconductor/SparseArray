
test_that(".tune_SVT_SparseArray_dims()", {
    .tune_SVT_SparseArray_dims <- SparseArray:::.tune_SVT_SparseArray_dims

    a <- array(0L, dim=c(5, 4, 1, 3))
    dimnames(a) <- list(letters[1:5], NULL, NULL, LETTERS[1:3])
    a[c(1:2, 8, 10, 15:17, 20, 24, 40)] <- (1:10)*10L
    svt <- SparseArray(a)

    dim_tuner <- c(0L, 0L, 0L, 0L)
    current <- .tune_SVT_SparseArray_dims(svt, dim_tuner)
    expect_identical(current, svt)

    dim_tuner <- c(0L, 0L, -1L, 0L)
    current <- .tune_SVT_SparseArray_dims(svt, dim_tuner)
    expect_identical(as.array(current), drop(a))
    svt2 <- .tune_SVT_SparseArray_dims(current, -dim_tuner)
    expect_identical(svt2, svt)
  
    dim_tuner <- c(0L, 0L, 1L, 0L, 1L, 0L)
    current <- .tune_SVT_SparseArray_dims(svt, dim_tuner)
    expected <- `dim<-`(a, c(5, 4, 1, 1, 1, 3))
    dimnames(expected)[c(1, 2, 4, 6)] <- dimnames(a)
    expect_identical(as.array(current), expected)
    svt2 <- .tune_SVT_SparseArray_dims(current, -dim_tuner)
    expect_identical(svt2, svt)

    dim_tuner <- c(1L, 1L, 0L, 1L, 0L, -1L, 0L, 1L)
    current <- .tune_SVT_SparseArray_dims(svt, dim_tuner)
    expected <- `dim<-`(a, c(1, 1, 5, 1, 4, 3, 1))
    dimnames(expected)[c(3, 5, 6)] <- dimnames(a)[c(1, 2, 4)]
    expect_identical(as.array(current), expected)
    svt2 <- .tune_SVT_SparseArray_dims(current, -dim_tuner)
    expect_identical(svt2, svt)
})

