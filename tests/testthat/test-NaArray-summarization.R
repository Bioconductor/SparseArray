
test_that("anyNA() method for NaArray objects", {
    ## input of type() "integer"
    m1 <- matrix(c(0L, 0L, 155L,
                   0L, 8L,  -1L), nrow=2, byrow=TRUE)
    naa1 <- as(m1, "NaArray")
    test_summarize_op1(m1, naa1, "anyNA")
    m1[1, 2] <- NA
    naa1 <- as(m1, "NaArray")
    test_summarize_op1(m1, naa1, "anyNA")

    ## input of type() "logical"
    m2 <- matrix(c(FALSE, FALSE, TRUE,
                   FALSE,  TRUE, TRUE), nrow=2, byrow=TRUE)
    naa2 <- as(m2, "NaArray")
    test_summarize_op1(m2, naa2, "anyNA")
    m2[1, 2] <- NA
    naa2 <- as(m2, "NaArray")
    test_summarize_op1(m2, naa2, "anyNA")

    ## input of type() "double"
    m3 <- matrix(c(0,    0,  pi,
                   0, 0.25, 1e3), nrow=2, byrow=TRUE)
    naa3 <- as(m3, "NaArray")
    test_summarize_op1(m3, naa3, "anyNA")
    m3[1, 2] <- NaN
    naa3 <- as(m3, "NaArray")
    test_summarize_op1(m3, naa3, "anyNA")
    m3[1, 2] <- NA
    naa3 <- as(m3, "NaArray")
    test_summarize_op1(m3, naa3, "anyNA")

    ## input of type() "complex"
    m4 <- matrix(c(0,    0,  pi,
                   0, 2-5i, 1e3), nrow=2, byrow=TRUE)
    naa4 <- as(m4, "NaArray")
    test_summarize_op1(m4, naa4, "anyNA")
    m4[1, 2] <- NaN       # 1st type of "complex" NaN
    naa4 <- as(m4, "NaArray")
    test_summarize_op1(m4, naa4, "anyNA")
    m4[1, 2] <- NaN * 1i  # 2nd type of "complex" NaN
    naa4 <- as(m4, "NaArray")
    test_summarize_op1(m4, naa4, "anyNA")
    m4[1, 2] <- NA
    naa4 <- as(m4, "NaArray")
    test_summarize_op1(m4, naa4, "anyNA")

    ## input of type() "character"
    m5 <- matrix(c("",     "", "Hello",
                   "", "dear", "world"), nrow=2, byrow=TRUE)
    naa5 <- as(m5, "NaArray")
    test_summarize_op1(m5, naa5, "anyNA")
    m5[1, 2] <- NA
    naa5 <- as(m5, "NaArray")
    test_summarize_op1(m5, naa5, "anyNA")
})

test_that("other summarization methods for NaArray objects", {
    ## input of type() "integer"
    m1 <- matrix(c( 0L, 0L,  NA, 0L, NA,
                    NA, 0L, -3L, 1L, NA,
                    0L, 0L,  0L, 0L, 0L,
                   15L, 0L,  0L, 0L, NA), nrow=4, byrow=TRUE)
    naa1 <- as(m1, "NaArray")
    test_summarize_op2(m1, naa1, "any")
    test_summarize_op2(m1, naa1, "all")
    test_summarize_op2(m1, naa1, "min")
    test_summarize_op2(m1, naa1, "max")
    test_summarize_op2(m1, naa1, "range")
    test_summarize_op2(m1, naa1, "sum")
    test_summarize_op2(m1, naa1, "prod")
    test_summarize_op2(m1, naa1, "mean")
    test_summarize_op2(m1, naa1, "var")
    test_summarize_op2(m1, naa1, "sd")
    m0 <- m1[0, ]
    naa0 <- naa1[0, ]
    expect_warning(min(naa0), "NAs introduced")
    expect_warning(max(naa0), "NAs introduced")
    expect_warning(range(naa0), "NAs introduced")
    expect_identical(suppressWarnings(min(naa0)), NA_integer_)
    expect_identical(suppressWarnings(max(naa0)), NA_integer_)
    expect_identical(suppressWarnings(range(naa0)), rep(NA_integer_,2))

    ## input of type() "logical"
    m2 <- is.na(m1)
    naa2 <- as(m2, "NaArray")
    test_summarize_op2(m2, naa2, "any")
    test_summarize_op2(m2, naa2, "all")
    test_summarize_op2(m2, naa2, "min")
    test_summarize_op2(m2, naa2, "max")
    test_summarize_op2(m2, naa2, "range")
    test_summarize_op2(m2, naa2, "sum")
    test_summarize_op2(m2, naa2, "prod")
    test_summarize_op2(m2, naa2, "mean")
    test_summarize_op2(m2, naa2, "var")
    test_summarize_op2(m2, naa2, "sd")
    m0 <- m2[0, ]
    naa0 <- naa2[0, ]
    expect_warning(min(naa0), "NAs introduced")
    expect_warning(max(naa0), "NAs introduced")
    expect_warning(range(naa0), "NAs introduced")
    expect_identical(suppressWarnings(min(naa0)), NA_integer_)
    expect_identical(suppressWarnings(max(naa0)), NA_integer_)
    expect_identical(suppressWarnings(range(naa0)), rep(NA_integer_,2))
})

test_that("summarization methods for 3D NaArray objects", {
    ## input of type() "double"
    a <- array(0, 6:4)
    a[1, , 2] <- c(1e12, -1234.55, -2.1, -1, -0.55)
    a[3, , 2] <- c(-0.55, 0, 1e-10, 0.88, 1)
    a[5, , 2] <- c(pi, 10.33, 3.4567895e8, 300, 2009.01)
    naa3 <- as(a, "NaArray")
    test_summarize_op1(a, naa3, "anyNA")
    a[6, 3, 2] <- NA
    a[6, 4, 2] <- NaN
    naa3 <- as(a, "NaArray")
    test_summarize_op1(a, naa3, "anyNA")
    expect_error(any(naa3), "does not support")
    expect_error(all(naa3), "does not support")
    test_summarize_op2(a, naa3, "min")
    test_summarize_op2(a, naa3, "max")
    test_summarize_op2(a, naa3, "range")
    test_summarize_op2(a, naa3, "sum")
    test_summarize_op2(a, naa3, "prod")
    test_summarize_op2(a, naa3, "mean")
    test_summarize_op2(a, naa3, "var")
    test_summarize_op2(a, naa3, "sd")
})

