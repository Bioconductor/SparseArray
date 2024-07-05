.test_isFUN <- function(a, svt, isFUN)
{
    FUN <- match.fun(isFUN)
    expected <- FUN(a)
    current  <- FUN(svt)
    check_SparseArray_object(current, class(svt), expected)
    expect_identical(as(expected, "SVT_SparseArray"), current)
}

test_that("is.na(), is.nan(), and is.infinite() on SVT_SparseArray objects", {
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

test_that("pmin() and pmax() on SVT_SparseArray objects", {
     set.seed(2009)

     ## --- no NAs ---
     svt1 <- poissonSparseMatrix(6, 5, density=0.5,
                                 dimnames=list(letters[1:6], NULL))
     svt2 <- poissonSparseMatrix(6, 5, density=0.6,
                                 dimnames=list(NULL, LETTERS[1:5])) * 0.77

     expect_identical(pmin(svt1), svt1)
     expect_identical(pmin(svt1, svt1), svt1)
     expect_identical(pmin(svt1, svt1, svt1), svt1)
     expect_identical(pmin(svt2), svt2)
     expect_identical(pmin(svt2, svt2), svt2)
     expect_identical(pmin(svt2, svt2, svt2), svt2)
     expect_identical(pmax(svt1), svt1)
     expect_identical(pmax(svt1, svt1), svt1)
     expect_identical(pmax(svt1, svt1, svt1), svt1)
     expect_identical(pmax(svt2), svt2)
     expect_identical(pmax(svt2, svt2), svt2)
     expect_identical(pmax(svt2, svt2, svt2), svt2)

     a1 <- as.array(svt1)
     a2 <- as.array(svt2)

     expected <- pmin(a1, a2)
     current <- pmin(svt1, svt2)
     check_SparseArray_object(current, "SVT_SparseMatrix", expected)
     expected <- pmax(a1, a2)
     current <- pmax(svt1, svt2)
     check_SparseArray_object(current, "SVT_SparseMatrix", expected)
     ## Swaping the order produces the same results **except** for the
     ## dimnames!
     expected <- pmin(a2, a1)
     current <- pmin(svt2, svt1)
     check_SparseArray_object(current, "SVT_SparseMatrix", expected)
     expected <- pmax(a2, a1)
     current <- pmax(svt2, svt1)
     check_SparseArray_object(current, "SVT_SparseMatrix", expected)

     ## --- with NAs ---

     svt1[cbind(6, 1:5)] <- c(NA, -2, NA, -4, NA)
     svt2[cbind(6, 1:5)] <- c(NA, NA, 0, NaN, NaN)
     a1 <- as.array(svt1)
     a2 <- as.array(svt2)

     ## Handling of NAs/NaNs differ slightly between base::pmin()/pmax()
     ## and our methods for SVT_SparseMatrix objects e.g. whether we get an
     ## NA or NaN in the output depends on the order of the arguments.
     ## Let's fix that:
     .fix_NAs_for_no_narm <- function(ans, x, y) {
         NA_idx <- which((is.na(x) & !is.nan(x)) | ((is.na(y) & !is.nan(y))))
         ans[NA_idx] <- NA
         ans
     }

     expected <- .fix_NAs_for_no_narm(pmin(a1, a2), a1, a2)
     current <- .fix_NAs_for_no_narm(pmin(svt1, svt2), a1, a2)
     check_SparseArray_object(current, "SVT_SparseMatrix", expected)
     expected <- .fix_NAs_for_no_narm(pmax(a1, a2), a1, a2)
     current <- .fix_NAs_for_no_narm(pmax(svt1, svt2), a1, a2)
     check_SparseArray_object(current, "SVT_SparseMatrix", expected)
     ## Swaping the order produces the same results **except** for the
     ## dimnames!
     expected <- .fix_NAs_for_no_narm(pmin(a2, a1), a2, a1)
     current <- .fix_NAs_for_no_narm(pmin(svt2, svt1), a2, a1)
     check_SparseArray_object(current, "SVT_SparseMatrix", expected)
     expected <- .fix_NAs_for_no_narm(pmax(a2, a1), a2, a1)
     current <- .fix_NAs_for_no_narm(pmax(svt2, svt1), a2, a1)
     check_SparseArray_object(current, "SVT_SparseMatrix", expected)

     .fix_NAs_for_narm <- function(ans, x, y) {
         NA_idx <- which(is.na(x) & is.nan(y) | is.nan(x) & is.na(y))
         ans[NA_idx] <- NaN
         ans
     }

     expected <- .fix_NAs_for_narm(pmin(a1, a2, na.rm=TRUE), a1, a2)
     current <- .fix_NAs_for_narm(pmin(svt1, svt2, na.rm=TRUE), a1, a2)
     check_SparseArray_object(current, "SVT_SparseMatrix", expected)
     expected <- .fix_NAs_for_narm(pmax(a1, a2, na.rm=TRUE), a1, a2)
     current <- .fix_NAs_for_narm(pmax(svt1, svt2, na.rm=TRUE), a1, a2)
     check_SparseArray_object(current, "SVT_SparseMatrix", expected)
     ## Swaping the order produces the same results **except** for the
     ## dimnames!
     expected <- .fix_NAs_for_narm(pmin(a2, a1, na.rm=TRUE), a2, a1)
     current <- .fix_NAs_for_narm(pmin(svt2, svt1, na.rm=TRUE), a2, a1)
     check_SparseArray_object(current, "SVT_SparseMatrix", expected)
     expected <- .fix_NAs_for_narm(pmax(a2, a1, na.rm=TRUE), a2, a1)
     current <- .fix_NAs_for_narm(pmax(svt2, svt1, na.rm=TRUE), a2, a1)
     check_SparseArray_object(current, "SVT_SparseMatrix", expected)
})

