
.test_summarize_op1 <- function(a, svt, op)
{
    FUN <- match.fun(op)
    expected <- FUN(as.vector(a))
    current <- FUN(svt)
    expect_identical(current, expected)
}

.test_summarize_op2 <- function(a, svt, op)
{
    FUN <- match.fun(op)
    if (op %in% c("var", "sd") ||
        is.double(a) &&  op %in% c("sum", "prod", "mean"))
    {
        EXPECT_FUN <- expect_equal
    } else {
        EXPECT_FUN <- expect_identical
    }
    expected <- FUN(as.vector(a))
    current <- FUN(svt)
    if (op == "prod" && is.integer(current))
        expected <- as.integer(expected)
    EXPECT_FUN(current, expected)
    expected <- FUN(as.vector(a), na.rm=TRUE)
    current <- FUN(svt, na.rm=TRUE)
    if (op == "prod" && is.integer(current))
        expected <- as.integer(expected)
    EXPECT_FUN(current, expected)
}

test_that("anyNA() method for SVT_SparseArray objects", {
    ## input of type() "integer"
    m1 <- matrix(c(0L, 0L, 155L,
                   0L, 8L,  -1L), nrow=2, byrow=TRUE)
    svt1 <- as(m1, "SVT_SparseArray")
    .test_summarize_op1(m1, svt1, "anyNA")
    m1[1, 2] <- NA
    svt1 <- as(m1, "SVT_SparseArray")
    .test_summarize_op1(m1, svt1, "anyNA")

    ## input of type() "logical"
    m2 <- matrix(c(FALSE, FALSE, TRUE,
                   FALSE,  TRUE, TRUE), nrow=2, byrow=TRUE)
    svt2 <- as(m2, "SVT_SparseArray")
    .test_summarize_op1(m2, svt2, "anyNA")
    m2[1, 2] <- NA
    svt2 <- as(m2, "SVT_SparseArray")
    .test_summarize_op1(m2, svt2, "anyNA")

    ## input of type() "double"
    m3 <- matrix(c(0,    0,  pi,
                   0, 0.25, 1e3), nrow=2, byrow=TRUE)
    svt3 <- as(m3, "SVT_SparseArray")
    .test_summarize_op1(m3, svt3, "anyNA")
    m3[1, 2] <- NaN
    svt3 <- as(m3, "SVT_SparseArray")
    .test_summarize_op1(m3, svt3, "anyNA")
    m3[1, 2] <- NA
    svt3 <- as(m3, "SVT_SparseArray")
    .test_summarize_op1(m3, svt3, "anyNA")

    ## input of type() "complex"
    m4 <- matrix(c(0,    0,  pi,
                   0, 2-5i, 1e3), nrow=2, byrow=TRUE)
    svt4 <- as(m4, "SVT_SparseArray")
    .test_summarize_op1(m4, svt4, "anyNA")
    m4[1, 2] <- NaN       # 1st type of "complex" NaN
    svt4 <- as(m4, "SVT_SparseArray")
    .test_summarize_op1(m4, svt4, "anyNA")
    m4[1, 2] <- NaN * 1i  # 2nd type of "complex" NaN
    svt4 <- as(m4, "SVT_SparseArray")
    .test_summarize_op1(m4, svt4, "anyNA")
    m4[1, 2] <- NA
    svt4 <- as(m4, "SVT_SparseArray")
    .test_summarize_op1(m4, svt4, "anyNA")

    ## input of type() "character"
    m5 <- matrix(c("",     "", "Hello",
                   "", "dear", "world"), nrow=2, byrow=TRUE)
    svt5 <- as(m5, "SVT_SparseArray")
    .test_summarize_op1(m5, svt5, "anyNA")
    m5[1, 2] <- NA
    svt5 <- as(m5, "SVT_SparseArray")
    .test_summarize_op1(m5, svt5, "anyNA")
})

test_that("other summarization methods for SVT_SparseArray objects", {
    ## input of type() "integer"
    m1 <- matrix(c( 0L, 0L,  NA, 0L, NA,
                    NA, 0L, -3L, 1L, NA,
                    0L, 0L,  0L, 0L, 0L,
                   15L, 0L,  0L, 0L, NA), nrow=4, byrow=TRUE)
    svt1 <- as(m1, "SVT_SparseArray")
    .test_summarize_op2(m1, svt1, "any")
    .test_summarize_op2(m1, svt1, "all")
    .test_summarize_op2(m1, svt1, "min")
    .test_summarize_op2(m1, svt1, "max")
    .test_summarize_op2(m1, svt1, "range")
    .test_summarize_op2(m1, svt1, "sum")
    .test_summarize_op2(m1, svt1, "prod")
    .test_summarize_op2(m1, svt1, "mean")
    .test_summarize_op2(m1, svt1, "var")
    .test_summarize_op2(m1, svt1, "sd")
    m0 <- m1[0, ]
    expect_warning(min(SparseArray(m0)), "NAs introduced")
    expect_identical(suppressWarnings(min(SparseArray(m0))), NA_integer_)
    expect_warning(max(SparseArray(m0)), "NAs introduced")
    expect_identical(suppressWarnings(max(SparseArray(m0))), NA_integer_)
    expect_warning(range(SparseArray(m0)), "NAs introduced")
    expect_identical(suppressWarnings(range(SparseArray(m0))),
                     rep.int(NA_integer_,2))

    ## input of type() "logical"
    m2 <- is.na(m1)
    svt2 <- as(m2, "SVT_SparseArray")
    .test_summarize_op2(m2, svt2, "any")
    .test_summarize_op2(m2, svt2, "all")
    .test_summarize_op2(m2, svt2, "min")
    .test_summarize_op2(m2, svt2, "max")
    .test_summarize_op2(m2, svt2, "range")
    .test_summarize_op2(m2, svt2, "sum")
    .test_summarize_op2(m2, svt2, "prod")
    .test_summarize_op2(m2, svt2, "mean")
    .test_summarize_op2(m2, svt2, "var")
    .test_summarize_op2(m2, svt2, "sd")
    m0 <- m2[0, ]
    expect_warning(min(SparseArray(m0)), "NAs introduced")
    expect_identical(suppressWarnings(min(SparseArray(m0))), NA_integer_)
    expect_warning(max(SparseArray(m0)), "NAs introduced")
    expect_identical(suppressWarnings(max(SparseArray(m0))), NA_integer_)
    expect_warning(range(SparseArray(m0)), "NAs introduced")
    expect_identical(suppressWarnings(range(SparseArray(m0))),
                     rep.int(NA_integer_,2))

    ## input of type() "double"
    a <- array(0, 6:4)
    a[1, , 2] <- c(1e12, -1234.55, -2.1, -1, -0.55)
    a[3, , 2] <- c(-0.55, 0, 1e-10, 0.88, 1)
    a[5, , 2] <- c(pi, 10.33, 3.4567895e8, 300, 2009.01)
    a[6, 3:4, 2] <- c(NA, NaN)
    svt3 <- as(a, "SVT_SparseArray")
    expect_error(any(svt3), "does not support")
    expect_error(all(svt3), "does not support")
    .test_summarize_op2(a, svt3, "min")
    .test_summarize_op2(a, svt3, "max")
    .test_summarize_op2(a, svt3, "range")
    .test_summarize_op2(a, svt3, "sum")
    .test_summarize_op2(a, svt3, "prod")
    .test_summarize_op2(a, svt3, "mean")
    .test_summarize_op2(a, svt3, "var")
    .test_summarize_op2(a, svt3, "sd")
})

