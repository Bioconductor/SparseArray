
.test_matrixStats_method1 <- function(a, svt, generic)
{
    coo <- as(svt, "COO_SparseArray")
    FUN <- match.fun(generic)
    expected <- FUN(a, useNames=FALSE)
    current <- FUN(svt, useNames=FALSE)
    expect_identical(current, expected)
    expect_identical(FUN(coo, useNames=FALSE), current)
    expected <- FUN(a, useNames=TRUE)
    current <- FUN(svt, useNames=TRUE)
    expect_identical(current, expected)
    expect_identical(FUN(coo, useNames=TRUE), current)
}

.test_matrixStats_method2 <- function(a, svt, generic, dims)
{
    coo <- as(svt, "COO_SparseArray")
    FUN <- match.fun(generic)
    op <- sub("^(col|row)", "", generic)
    if (op %in% c("Vars", "Sds") ||
        is.double(a) && op %in% c("Sums", "Prods", "Means", "Sums2", "Means2"))
    {
        EXPECT_FUN <- expect_equal
    } else {
        EXPECT_FUN <- expect_identical
    }
    if (op %in% c("Sums", "Means")) {
        if (missing(dims)) {
            ## No 'useNames' arg.
            expected <- FUN(a)
            current <- FUN(svt)
            EXPECT_FUN(current, expected)
            expect_identical(FUN(coo), current)
            expected <- FUN(a, na.rm=TRUE)
            current <- FUN(svt, na.rm=TRUE)
            EXPECT_FUN(current, expected)
            expect_identical(FUN(coo, na.rm=TRUE), current)
        } else {
            expected <- FUN(a, dims=dims)
            current <- FUN(svt, dims=dims)
            EXPECT_FUN(current, expected)
            expect_identical(FUN(coo, dims=dims), current)
            expected <- FUN(a, na.rm=TRUE, dims=dims)
            current <- FUN(svt, na.rm=TRUE, dims=dims)
            EXPECT_FUN(current, expected)
            expect_identical(FUN(coo, na.rm=TRUE, dims=dims), current)
        }
    } else {
        expected <- FUN(a, useNames=FALSE)
        current <- FUN(svt, useNames=FALSE)
        EXPECT_FUN(current, expected)
        expect_identical(FUN(coo, useNames=FALSE), current)
        expected <- FUN(a, na.rm=TRUE, useNames=FALSE)
        current <- FUN(svt, na.rm=TRUE, useNames=FALSE)
        EXPECT_FUN(current, expected)
        expect_identical(FUN(coo, na.rm=TRUE, useNames=FALSE), current)
        expected <- FUN(a, useNames=TRUE)
        current <- FUN(svt, useNames=TRUE)
        EXPECT_FUN(current, expected)
        expect_identical(FUN(coo, useNames=TRUE), current)
        expected <- FUN(a, na.rm=TRUE, useNames=TRUE)
        current <- FUN(svt, na.rm=TRUE, useNames=TRUE)
        EXPECT_FUN(current, expected)
        expect_identical(FUN(coo, na.rm=TRUE, useNames=TRUE), current)
    }
}

### Tests *Mins() and *Maxs() methods on a 3D array-like object.
.test_colrowMinsMaxs_3D <- function(object)
{
    a <- as.array(object)

    ## Base R does NOT allow an ordinay array with dimensions of extent
    ## zero to carry a character(0) in its dimnames, only a NULL. This means
    ## that, if 'object' is an Array derivative, then 'dimnames(object)'
    ## and 'dimnames(as.array(object))' are not guaranteed to be identical
    ## because a character(0) in the former will be replaced with a NULL in
    ## the latter.
    ## As a consequence, calling colMins/Maxs() or rowMins/Maxs()
    ## on 'object' won't necessarily produce the exact same result as
    ## calling simple_colMins3D/Maxs3D() or simple_rowMins3D/Maxs3D()
    ## on 'as.array(object)' when the result is a vector of length 0.
    ## More precisely, one can be named while the other is not.

    ## Does not look at the names if 'a1' and 'a2' are vectors of length 0.
    expect_almost_identical <- function(a1, a2) {
        if (is.vector(a1) && length(a1) == 0L) {
            a1 <- unname(a1)
            a2 <- unname(a2)
        }
        expect_identical(a1, a2)
    }

    ## dims == 1 (default)

    expected <- simple_colMins3D(a)
    expect_almost_identical(colMins(object), expected)
    expected <- simple_colMaxs3D(a)
    expect_almost_identical(colMaxs(object), expected)

    expected <- simple_rowMins3D(a)
    expect_almost_identical(rowMins(object), expected)
    expected <- simple_rowMaxs3D(a)
    expect_almost_identical(rowMaxs(object), expected)

    expected <- simple_colMins3D(a, na.rm=TRUE)
    expect_almost_identical(colMins(object, na.rm=TRUE), expected)
    expected <- simple_colMaxs3D(a, na.rm=TRUE)
    expect_almost_identical(colMaxs(object, na.rm=TRUE), expected)

    expected <- simple_rowMins3D(a, na.rm=TRUE)
    expect_almost_identical(rowMins(object, na.rm=TRUE), expected)
    expected <- simple_rowMaxs3D(a, na.rm=TRUE)
    expect_almost_identical(rowMaxs(object, na.rm=TRUE), expected)

    ## dims == 2

    expected <- simple_colMins3D(a, dims=2)
    expect_almost_identical(colMins(object, dims=2), expected)
    expected <- simple_colMaxs3D(a, dims=2)
    expect_almost_identical(colMaxs(object, dims=2), expected)

    expected <- simple_rowMins3D(a, dims=2)
    expect_almost_identical(rowMins(object, dims=2), expected)
    expected <- simple_rowMaxs3D(a, dims=2)
    expect_almost_identical(rowMaxs(object, dims=2), expected)

    expected <- simple_colMins3D(a, na.rm=TRUE, dims=2)
    expect_almost_identical(colMins(object, na.rm=TRUE, dims=2), expected)
    expected <- simple_colMaxs3D(a, na.rm=TRUE, dims=2)
    expect_almost_identical(colMaxs(object, na.rm=TRUE, dims=2), expected)

    expected <- simple_rowMins3D(a, na.rm=TRUE, dims=2)
    expect_almost_identical(rowMins(object, na.rm=TRUE, dims=2), expected)
    expected <- simple_rowMaxs3D(a, na.rm=TRUE, dims=2)
    expect_almost_identical(rowMaxs(object, na.rm=TRUE, dims=2), expected)
}

test_that("colAnyNAs/rowAnyNAs() methods for 2D SparseArray objects", {
    ## input of type() "integer"
    m1 <- matrix(c(0L, 0L, 155L,
                   0L, 8L,  -1L), nrow=2, byrow=TRUE,
                 dimnames=list(LETTERS[1:2], letters[1:3]))
    svt1 <- as(m1, "SVT_SparseArray")
    .test_matrixStats_method1(m1, svt1, "colAnyNAs")
    .test_matrixStats_method1(m1, svt1, "rowAnyNAs")
    m1[1, 2] <- NA
    svt1 <- as(m1, "SVT_SparseArray")
    .test_matrixStats_method1(m1, svt1, "colAnyNAs")
    .test_matrixStats_method1(m1, svt1, "rowAnyNAs")

    ## input of type() "logical"
    m2 <- matrix(c(FALSE, FALSE, TRUE,
                   FALSE,  TRUE, TRUE), nrow=2, byrow=TRUE,
                 dimnames=list(LETTERS[1:2], letters[1:3]))
    svt2 <- as(m2, "SVT_SparseArray")
    .test_matrixStats_method1(m2, svt2, "colAnyNAs")
    .test_matrixStats_method1(m2, svt2, "rowAnyNAs")
    m2[1, 2] <- NA
    svt2 <- as(m2, "SVT_SparseArray")
    .test_matrixStats_method1(m2, svt2, "colAnyNAs")
    .test_matrixStats_method1(m2, svt2, "rowAnyNAs")

    ## input of type() "double"
    m3 <- matrix(c(0,    0,  pi,
                   0, 0.25, 1e3), nrow=2, byrow=TRUE,
                 dimnames=list(LETTERS[1:2], letters[1:3]))
    svt3 <- as(m3, "SVT_SparseArray")
    .test_matrixStats_method1(m3, svt3, "colAnyNAs")
    .test_matrixStats_method1(m3, svt3, "rowAnyNAs")
    m3[1, 2] <- svt3[1, 2] <- NaN
    .test_matrixStats_method1(m3, svt3, "colAnyNAs")
    .test_matrixStats_method1(m3, svt3, "rowAnyNAs")
    m3[1, 2] <- svt3[1, 2] <- NA
    .test_matrixStats_method1(m3, svt3, "colAnyNAs")
    .test_matrixStats_method1(m3, svt3, "rowAnyNAs")

    ## input of type() "complex"
    m4 <- matrix(c(0,    0,  pi,
                   0, 2-5i, 1e3), nrow=2, byrow=TRUE,
                 dimnames=list(LETTERS[1:2], letters[1:3]))
    svt4 <- as(m4, "SVT_SparseArray")
    .test_matrixStats_method1(m4, svt4, "colAnyNAs")
    .test_matrixStats_method1(m4, svt4, "rowAnyNAs")
    m4[1, 2] <- NaN       # 1st type of "complex" NaN
    svt4 <- as(m4, "SVT_SparseArray")
    .test_matrixStats_method1(m4, svt4, "colAnyNAs")
    .test_matrixStats_method1(m4, svt4, "rowAnyNAs")
    m4[1, 2] <- NaN * 1i  # 2nd type of "complex" NaN
    svt4 <- as(m4, "SVT_SparseArray")
    .test_matrixStats_method1(m4, svt4, "colAnyNAs")
    .test_matrixStats_method1(m4, svt4, "rowAnyNAs")
    m4[1, 2] <- NA
    svt4 <- as(m4, "SVT_SparseArray")
    .test_matrixStats_method1(m4, svt4, "colAnyNAs")
    .test_matrixStats_method1(m4, svt4, "rowAnyNAs")

    ## input of type() "character"
    m5 <- matrix(c("",     "", "Hello",
                   "", "dear", "world"), nrow=2, byrow=TRUE,
                 dimnames=list(LETTERS[1:2], letters[1:3]))
    svt5 <- as(m5, "SVT_SparseArray")
    .test_matrixStats_method1(m5, svt5, "colAnyNAs")
    .test_matrixStats_method1(m5, svt5, "rowAnyNAs")
    m5[1, 2] <- NA
    svt5 <- as(m5, "SVT_SparseArray")
    .test_matrixStats_method1(m5, svt5, "colAnyNAs")
    .test_matrixStats_method1(m5, svt5, "rowAnyNAs")
})

test_that("other matrixStats methods for 2D SparseArray objects", {
    ## input of type() "integer"
    m1 <- matrix(c( 0L, 0L,  NA, 0L, NA,
                    NA, 0L, -3L, 1L, NA,
                    0L, 0L,  0L, 0L, 0L,
                   15L, 0L,  0L, 0L, NA), nrow=4, byrow=TRUE,
                 dimnames=list(LETTERS[1:4], letters[1:5]))
    svt1 <- as(m1, "SVT_SparseArray")
    .test_matrixStats_method2(m1, svt1, "colAnys")
    .test_matrixStats_method2(m1, svt1, "rowAnys")
    .test_matrixStats_method2(m1, svt1, "colAlls")
    .test_matrixStats_method2(m1, svt1, "rowAlls")
    .test_matrixStats_method2(m1, svt1, "colMins")
    .test_matrixStats_method2(m1, svt1, "rowMins")
    .test_matrixStats_method2(m1, svt1, "colMaxs")
    .test_matrixStats_method2(m1, svt1, "rowMaxs")
    .test_matrixStats_method2(m1, svt1, "colRanges")
    .test_matrixStats_method2(m1, svt1, "rowRanges")
    .test_matrixStats_method2(m1, svt1, "colSums")
    .test_matrixStats_method2(m1, svt1, "rowSums")
    .test_matrixStats_method2(m1, svt1, "colProds")
    .test_matrixStats_method2(m1, svt1, "rowProds")
    .test_matrixStats_method2(m1, svt1, "colMeans")
    .test_matrixStats_method2(m1, svt1, "rowMeans")
    .test_matrixStats_method2(m1, svt1, "colSums2")
    .test_matrixStats_method2(m1, svt1, "rowSums2")
    .test_matrixStats_method2(m1, svt1, "colMeans2")
    .test_matrixStats_method2(m1, svt1, "rowMeans2")
    .test_matrixStats_method2(m1, svt1, "colVars")
    .test_matrixStats_method2(m1, svt1, "rowVars")
    .test_matrixStats_method2(m1, svt1, "colSds")
    .test_matrixStats_method2(m1, svt1, "rowSds")
    m0 <- m1[0, ]
    svt0 <- svt1[0, ]
    expected <- rep(NA_integer_, 5L)
    expect_warning(colMins(svt0), "NAs introduced")
    expect_warning(colMaxs(svt0), "NAs introduced")
    expect_warning(colRanges(svt0), "NAs introduced")
    expect_identical(suppressWarnings(colMins(svt0, useNames=FALSE)), expected)
    expect_identical(suppressWarnings(colMaxs(svt0, useNames=FALSE)), expected)
    expect_identical(suppressWarnings(colRanges(svt0, useNames=FALSE)),
                     cbind(expected, expected, deparse.level=0))
    expect_identical(suppressWarnings(colMins(svt0)),
                     setNames(expected, colnames(m0)))
    expect_identical(suppressWarnings(colMaxs(svt0)),
                     setNames(expected, colnames(m0)))
    expect_identical(suppressWarnings(colRanges(svt0)),
                     cbind(setNames(expected, colnames(m0)), expected,
                           deparse.level=0))
    expect_identical(unname(rowMins(svt0)), rowMins(m0))
    expect_identical(unname(rowMaxs(svt0)), rowMaxs(m0))
    expect_identical(rowRanges(svt0), rowRanges(m0))

    ## input of type() "logical"
    m2 <- is.na(m1)
    svt2 <- as(m2, "SVT_SparseArray")
    .test_matrixStats_method2(m2, svt2, "colAnys")
    .test_matrixStats_method2(m2, svt2, "rowAnys")
    .test_matrixStats_method2(m2, svt2, "colAlls")
    .test_matrixStats_method2(m2, svt2, "rowAlls")
    storage.mode(m2) <- "integer"
    .test_matrixStats_method2(m2, svt2, "colMins")
    .test_matrixStats_method2(m2, svt2, "rowMins")
    .test_matrixStats_method2(m2, svt2, "colMaxs")
    .test_matrixStats_method2(m2, svt2, "rowMaxs")
    .test_matrixStats_method2(m2, svt2, "colRanges")
    .test_matrixStats_method2(m2, svt2, "rowRanges")
    .test_matrixStats_method2(m2, svt2, "colSums")
    .test_matrixStats_method2(m2, svt2, "rowSums")
    .test_matrixStats_method2(m2, svt2, "colProds")
    .test_matrixStats_method2(m2, svt2, "rowProds")
    .test_matrixStats_method2(m2, svt2, "colMeans")
    .test_matrixStats_method2(m2, svt2, "rowMeans")
    .test_matrixStats_method2(m2, svt2, "colSums2")
    .test_matrixStats_method2(m2, svt2, "rowSums2")
    .test_matrixStats_method2(m2, svt2, "colMeans2")
    .test_matrixStats_method2(m2, svt2, "rowMeans2")
    .test_matrixStats_method2(m2, svt2, "colVars")
    .test_matrixStats_method2(m2, svt2, "rowVars")
    .test_matrixStats_method2(m2, svt2, "colSds")
    .test_matrixStats_method2(m2, svt2, "rowSds")
    m0 <- m2[0, ]
    svt0 <- svt2[0, ]
    expected <- rep(NA_integer_, 5L)
    expect_warning(colMins(svt0), "NAs introduced")
    expect_warning(colMaxs(svt0), "NAs introduced")
    expect_warning(colRanges(svt0), "NAs introduced")
    expect_identical(suppressWarnings(colMins(svt0, useNames=FALSE)), expected)
    expect_identical(suppressWarnings(colMaxs(svt0, useNames=FALSE)), expected)
    expect_identical(suppressWarnings(colRanges(svt0, useNames=FALSE)),
                     cbind(expected, expected, deparse.level=0))
    expect_identical(suppressWarnings(colMins(svt0)),
                     setNames(expected, colnames(m0)))
    expect_identical(suppressWarnings(colMaxs(svt0)),
                     setNames(expected, colnames(m0)))
    expect_identical(suppressWarnings(colRanges(svt0)),
                     cbind(setNames(expected, colnames(m0)), expected,
                           deparse.level=0))
    expect_identical(unname(rowMins(svt0)), rowMins(m0))
    expect_identical(unname(rowMaxs(svt0)), rowMaxs(m0))
    expect_identical(rowRanges(svt0), rowRanges(m0))
})

test_that("matrixStats methods for 3D SparseArray objects", {
    ## input of type() "double"
    a <- array(0, 6:4,
               dimnames=list(letters[1:6], letters[22:26], LETTERS[1:4]))
    a[1, , 2] <- c(1e12, -1234.55, -2.1, -1, -0.55)
    a[3, , 2] <- c(-0.55, 0, 1e-10, 0.88, 1)
    a[5, , 2] <- c(pi, 10.33, 3.4567895e8, 300, 2009.01)
    a[6, 3:4, 2] <- c(NA, NaN)
    svt3 <- as(a, "SVT_SparseArray")
    coo3 <- as(svt3, "COO_SparseArray")

    .test_colrowMinsMaxs_3D(svt3)
    .test_colrowMinsMaxs_3D(coo3)

    ## dims == 1 (default)
    .test_matrixStats_method2(a, svt3, "colSums")
    .test_matrixStats_method2(a, svt3, "rowSums")
    .test_matrixStats_method2(a, svt3, "colMeans")
    .test_matrixStats_method2(a, svt3, "rowMeans")

    ## dims == 2
    .test_matrixStats_method2(a, svt3, "colSums", dims=2)
    .test_matrixStats_method2(a, svt3, "rowSums", dims=2)
    .test_matrixStats_method2(a, svt3, "colMeans", dims=2)
    .test_matrixStats_method2(a, svt3, "rowMeans", dims=2)

})

test_that("more torturing of the *Mins() and *Maxs() methods for SparseArray", {
    ## We use a 3D SVT_SparseArray object for the torture.

    ## input of type() "integer"
    svt0 <- SVT_SparseArray(dim=5:3,
                dimnames=list(letters[1:5], letters[23:26], LETTERS[1:3]),
                type="integer")
    svt0[c(1, 6, 16, 20:22, 36, 39:40, 60)] <-
                c(2L, -5L, NA, NA, -11L, 99L, -8L, NA, NA, NA)

    suppressWarnings(.test_colrowMinsMaxs_3D(svt0))
    expect_warning(rowMins(svt0, na.rm=TRUE, dims=2), "NAs introduced")
    expect_warning(rowMaxs(svt0, na.rm=TRUE, dims=2), "NAs introduced")

    svt <- svt0[ , , 0]
    suppressWarnings(.test_colrowMinsMaxs_3D(svt))
    expect_warning(rowMins(svt), "NAs introduced")
    expect_warning(rowMaxs(svt), "NAs introduced")
    expect_warning(rowMins(svt, dims=2), "NAs introduced")
    expect_warning(rowMaxs(svt, dims=2), "NAs introduced")

    svt <- svt0[ , 0, ]
    suppressWarnings(.test_colrowMinsMaxs_3D(svt))
    expect_warning(rowMins(svt), "NAs introduced")
    expect_warning(rowMaxs(svt), "NAs introduced")
    expect_warning(colMins(svt, dims=2), "NAs introduced")
    expect_warning(colMaxs(svt, dims=2), "NAs introduced")

    svt <- svt0[ 0, , ]
    suppressWarnings(.test_colrowMinsMaxs_3D(svt))
    expect_warning(colMins(svt), "NAs introduced")
    expect_warning(colMaxs(svt), "NAs introduced")
    expect_warning(colMins(svt, dims=2), "NAs introduced")
    expect_warning(colMaxs(svt, dims=2), "NAs introduced")

    ## input of type() "double"
    svt0[39:40] <- c(NaN, NaN)
    .test_colrowMinsMaxs_3D(svt0)
    .test_colrowMinsMaxs_3D(svt0[ , , 0])
    .test_colrowMinsMaxs_3D(svt0[ , 0, ])
    .test_colrowMinsMaxs_3D(svt0[ 0, , ])
})

