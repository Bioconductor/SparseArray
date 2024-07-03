.test_isFUN <- function(a, svt, isFUN)
{
    FUN <- match.fun(isFUN)
    expected <- FUN(a)
    current  <- FUN(svt)
    check_SparseArray_object(current, class(svt), expected)
    expect_identical(as(expected, "SVT_SparseArray"), current)
}

test_that("is.na() and is.nan() on SVT_SparseArray objects", {
    a <- array(0, 6:4)
    a[1, , 2] <- c(-Inf, -1234.55, -2.1, -1, -0.55)
    a[3, , 2] <- c(-0.55, 0, 1e-10, 0.88, 1)
    a[5, , 2] <- c(pi, 10.33, 3.4567895e8, 1e300, Inf)
    a[6, 3:4, 2] <- c(NA, NaN)
    svt <- as(a, "SVT_SparseArray")
    for (isFUN in c("is.na", "is.nan", "is.infinite"))
        .test_isFUN(a, svt, isFUN)

    a[1, , 4] <- complex(real=0, imag=c(-Inf, -1234.55, -2.1, -1, -0.55))
    a[2, , 4] <- complex(real=c(  0,   0,  NA,  NA, 0),
                         imag=c(  0,  NA,   0,  NA, 0))
    a[3, , 4] <- complex(real=c(  0,   0, NaN, NaN, 0),
                         imag=c(  0, NaN,   0, NaN, 0))
    a[4, , 4] <- complex(real=c( NA,  NA, NaN, NaN, 0),
                         imag=c( NA, NaN,  NA, NaN, 0))
    a[5, , 4] <- complex(real=c(Inf,  NA, Inf,  NA, 0),
                         imag=c(Inf, Inf,  NA,  NA, 0))
    a[6, , 4] <- complex(real=c(Inf, NaN, Inf, NaN, 0),
                         imag=c(Inf, Inf, NaN, NaN, 0))
    svt <- as(a, "SVT_SparseArray")
    for (isFUN in c("is.na", "is.nan", "is.infinite"))
        .test_isFUN(a, svt, isFUN)
})

