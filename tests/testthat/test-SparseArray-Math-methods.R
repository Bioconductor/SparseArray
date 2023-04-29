
.test_Math_op <- function(a, svt, op)
{
    FUN <- match.fun(op)
    ## Catch the warning if any. Unfortunately 'FUN(a)' is called twice.
    ## According to
    ##   https://stackoverflow.com/questions/3903157/how-can-i-check-whether-a-function-call-results-in-a-warning
    ## it's possible to avoid this by using handlers but I didn't try it.
    ## The handlers-based solution looks very ugly to start with.
    expected <- tryCatch(FUN(a), warning=identity)
    if (inherits(expected, "warning")) {
        expected_warning_msg <- expected$message
        expect_warning(FUN(svt), expected_warning_msg)
        expected <- suppressWarnings(FUN(a))
        current <- suppressWarnings(FUN(svt))
    } else {
        current <- FUN(svt)
    }
    expect_true(is(current, "SVT_SparseArray"))
    expect_true(validObject(current))
    expect_identical(as.array(current), expected)
}

test_that("'Math' ops on SVT_SparseArray objects", {
    a <- array(0, 6:4)
    a[1, , 2] <- c(-Inf, -1234.55, -2.1, -1, -0.55)
    a[3, , 2] <- c(-0.55, 0, 1e-10, 0.88, 1)
    a[5, , 2] <- c(pi, 10.33, 3.4567895e8, 1e300, Inf)
    a[6, 3:4, 2] <- c(NA, NaN)
    svt <- as(a, "SVT_SparseArray")

    ## 'Math' group (+ 'Math2' group, called with 'digits' argument missing).
    for (op in c(SparseArray:::.SUPPORTED_MATH_OPS, "round", "signif"))
        .test_Math_op(a, svt, op)

    ## 'Math2' group, called with various values of the 'digits' argument.
    for (op in c("round", "signif")) {
        for (digits in -6:7) {
            FUN <- match.fun(op)
            expected <- FUN(a, digits)
            current <- FUN(svt, digits)
            expect_true(is(current, "SVT_SparseArray"))
            expect_true(validObject(current))
            expect_identical(as.array(current), expected)
        }
    }
})

