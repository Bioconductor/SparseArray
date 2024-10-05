
.test_NaArray_matrixStats_method1 <- function(a, naa, generic)
{
    FUN <- match.fun(generic)
    expected <- FUN(a, useNames=FALSE)
    current <- FUN(naa, useNames=FALSE)
    expect_identical(current, expected)
    expected <- FUN(a, useNames=TRUE)
    current <- FUN(naa, useNames=TRUE)
    expect_identical(current, expected)
}

.test_NaArray_matrixStats_method2 <- function(a, naa, generic, dims)
{
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
            current <- FUN(naa)
            EXPECT_FUN(current, expected)
            expected <- FUN(a, na.rm=TRUE)
            current <- FUN(naa, na.rm=TRUE)
            EXPECT_FUN(current, expected)
        } else {
            expected <- FUN(a, dims=dims)
            current <- FUN(naa, dims=dims)
            EXPECT_FUN(current, expected)
            expected <- FUN(a, na.rm=TRUE, dims=dims)
            current <- FUN(naa, na.rm=TRUE, dims=dims)
            EXPECT_FUN(current, expected)
        }
    } else {
        expected <- FUN(a, useNames=FALSE)
        current <- FUN(naa, useNames=FALSE)
        EXPECT_FUN(current, expected)
        expected <- FUN(a, na.rm=TRUE, useNames=FALSE)
        current <- FUN(naa, na.rm=TRUE, useNames=FALSE)
        EXPECT_FUN(current, expected)
        expected <- FUN(a, useNames=TRUE)
        current <- FUN(naa, useNames=TRUE)
        EXPECT_FUN(current, expected)
        expected <- FUN(a, na.rm=TRUE, useNames=TRUE)
        current <- FUN(naa, na.rm=TRUE, useNames=TRUE)
        EXPECT_FUN(current, expected)
    }
}

test_that("colAnyNAs()/rowAnyNAs() methods for 2D NaArray objects", {
    ## input of type() "integer"
    m1 <- matrix(c(0L, 0L, 155L,
                   0L, 8L,  -1L), nrow=2, byrow=TRUE,
                 dimnames=list(LETTERS[1:2], letters[1:3]))
    nam1 <- as(m1, "NaArray")
    .test_NaArray_matrixStats_method1(m1, nam1, "colAnyNAs")
    .test_NaArray_matrixStats_method1(m1, nam1, "rowAnyNAs")
    m1[1, 2] <- NA
    nam1 <- as(m1, "NaArray")
    .test_NaArray_matrixStats_method1(m1, nam1, "colAnyNAs")
    .test_NaArray_matrixStats_method1(m1, nam1, "rowAnyNAs")

    ## input of type() "logical"
    m2 <- matrix(c(FALSE, FALSE, TRUE,
                   FALSE,  TRUE, TRUE), nrow=2, byrow=TRUE,
                 dimnames=list(LETTERS[1:2], letters[1:3]))
    nam2 <- as(m2, "NaArray")
    .test_NaArray_matrixStats_method1(m2, nam2, "colAnyNAs")
    .test_NaArray_matrixStats_method1(m2, nam2, "rowAnyNAs")
    m2[1, 2] <- NA
    nam2 <- as(m2, "NaArray")
    .test_NaArray_matrixStats_method1(m2, nam2, "colAnyNAs")
    .test_NaArray_matrixStats_method1(m2, nam2, "rowAnyNAs")

    ## input of type() "double"
    m3 <- matrix(c(0,    0,  pi,
                   0, 0.25, 1e3), nrow=2, byrow=TRUE,
                 dimnames=list(LETTERS[1:2], letters[1:3]))
    nam3 <- as(m3, "NaArray")
    .test_NaArray_matrixStats_method1(m3, nam3, "colAnyNAs")
    .test_NaArray_matrixStats_method1(m3, nam3, "rowAnyNAs")
    m3[1, 2] <- NaN
    nam3 <- as(m3, "NaArray")
    .test_NaArray_matrixStats_method1(m3, nam3, "colAnyNAs")
    .test_NaArray_matrixStats_method1(m3, nam3, "rowAnyNAs")
    m3[1, 2] <- NA
    nam3 <- as(m3, "NaArray")
    .test_NaArray_matrixStats_method1(m3, nam3, "colAnyNAs")
    .test_NaArray_matrixStats_method1(m3, nam3, "rowAnyNAs")

    ## input of type() "complex"
    m4 <- matrix(c(0,    0,  pi,
                   0, 2-5i, 1e3), nrow=2, byrow=TRUE,
                 dimnames=list(LETTERS[1:2], letters[1:3]))
    nam4 <- as(m4, "NaArray")
    .test_NaArray_matrixStats_method1(m4, nam4, "colAnyNAs")
    .test_NaArray_matrixStats_method1(m4, nam4, "rowAnyNAs")
    m4[1, 2] <- NaN       # 1st type of "complex" NaN
    nam4 <- as(m4, "NaArray")
    .test_NaArray_matrixStats_method1(m4, nam4, "colAnyNAs")
    .test_NaArray_matrixStats_method1(m4, nam4, "rowAnyNAs")
    m4[1, 2] <- NaN * 1i  # 2nd type of "complex" NaN
    nam4 <- as(m4, "NaArray")
    .test_NaArray_matrixStats_method1(m4, nam4, "colAnyNAs")
    .test_NaArray_matrixStats_method1(m4, nam4, "rowAnyNAs")
    m4[1, 2] <- NA
    nam4 <- as(m4, "NaArray")
    .test_NaArray_matrixStats_method1(m4, nam4, "colAnyNAs")
    .test_NaArray_matrixStats_method1(m4, nam4, "rowAnyNAs")

    ## input of type() "character"
    m5 <- matrix(c("",     "", "Hello",
                   "", "dear", "world"), nrow=2, byrow=TRUE,
                 dimnames=list(LETTERS[1:2], letters[1:3]))
    nam5 <- as(m5, "NaArray")
    .test_NaArray_matrixStats_method1(m5, nam5, "colAnyNAs")
    .test_NaArray_matrixStats_method1(m5, nam5, "rowAnyNAs")
    m5[1, 2] <- NA
    nam5 <- as(m5, "NaArray")
    .test_NaArray_matrixStats_method1(m5, nam5, "colAnyNAs")
    .test_NaArray_matrixStats_method1(m5, nam5, "rowAnyNAs")
})

test_that("other matrixStats methods for 2D NaArray objects", {
    ## input of type() "integer"
    m1 <- matrix(c( 0L, 0L,  NA, 0L, NA,
                    NA, 0L, -3L, 1L, NA,
                    0L, 0L,  0L, 0L, 0L,
                   15L, 0L,  0L, 0L, NA), nrow=4, byrow=TRUE,
                 dimnames=list(LETTERS[1:4], letters[1:5]))
    nam1 <- as(m1, "NaArray")
    .test_NaArray_matrixStats_method2(m1, nam1, "colAnys")
    #.test_NaArray_matrixStats_method2(m1, nam1, "rowAnys")
    .test_NaArray_matrixStats_method2(m1, nam1, "colAlls")
    #.test_NaArray_matrixStats_method2(m1, nam1, "rowAlls")
    .test_NaArray_matrixStats_method2(m1, nam1, "colMins")
    .test_NaArray_matrixStats_method2(m1, nam1, "rowMins")
    .test_NaArray_matrixStats_method2(m1, nam1, "colMaxs")
    .test_NaArray_matrixStats_method2(m1, nam1, "rowMaxs")
    .test_NaArray_matrixStats_method2(m1, nam1, "colRanges")
    .test_NaArray_matrixStats_method2(m1, nam1, "rowRanges")
    .test_NaArray_matrixStats_method2(m1, nam1, "colSums")
    .test_NaArray_matrixStats_method2(m1, nam1, "rowSums")
    .test_NaArray_matrixStats_method2(m1, nam1, "colProds")
    #.test_NaArray_matrixStats_method2(m1, nam1, "rowProds")
    .test_NaArray_matrixStats_method2(m1, nam1, "colMeans")
    #.test_NaArray_matrixStats_method2(m1, nam1, "rowMeans")
    .test_NaArray_matrixStats_method2(m1, nam1, "colSums2")
    .test_NaArray_matrixStats_method2(m1, nam1, "rowSums2")
    .test_NaArray_matrixStats_method2(m1, nam1, "colMeans2")
    #.test_NaArray_matrixStats_method2(m1, nam1, "rowMeans2")
    .test_NaArray_matrixStats_method2(m1, nam1, "colVars")
    #.test_NaArray_matrixStats_method2(m1, nam1, "rowVars")
    .test_NaArray_matrixStats_method2(m1, nam1, "colSds")
    #.test_NaArray_matrixStats_method2(m1, nam1, "rowSds")
    m0 <- m1[0, ]
    nam0 <- nam1[0, ]
    expected <- rep(NA_integer_, 5L)
    expect_warning(colMins(nam0), "NAs introduced")
    expect_warning(colMaxs(nam0), "NAs introduced")
    expect_warning(colRanges(nam0), "NAs introduced")
    expect_identical(suppressWarnings(colMins(nam0, useNames=FALSE)), expected)
    expect_identical(suppressWarnings(colMaxs(nam0, useNames=FALSE)), expected)
    expect_identical(suppressWarnings(colRanges(nam0, useNames=FALSE)),
                     cbind(expected, expected, deparse.level=0))
    expect_identical(suppressWarnings(colMins(nam0)),
                     setNames(expected, colnames(m0)))
    expect_identical(suppressWarnings(colMaxs(nam0)),
                     setNames(expected, colnames(m0)))
    expect_identical(suppressWarnings(colRanges(nam0)),
                     cbind(setNames(expected, colnames(m0)), expected,
                           deparse.level=0))
    expect_identical(unname(rowMins(nam0)), rowMins(m0))
    expect_identical(unname(rowMaxs(nam0)), rowMaxs(m0))
    expect_identical(rowRanges(nam0), rowRanges(m0))

    ## input of type() "logical"
    m2 <- is.na(m1)
    nam2 <- as(m2, "NaArray")
    .test_NaArray_matrixStats_method2(m2, nam2, "colAnys")
    #.test_NaArray_matrixStats_method2(m2, nam2, "rowAnys")
    .test_NaArray_matrixStats_method2(m2, nam2, "colAlls")
    #.test_NaArray_matrixStats_method2(m2, nam2, "rowAlls")
    storage.mode(m2) <- "integer"
    .test_NaArray_matrixStats_method2(m2, nam2, "colMins")
    .test_NaArray_matrixStats_method2(m2, nam2, "rowMins")
    .test_NaArray_matrixStats_method2(m2, nam2, "colMaxs")
    .test_NaArray_matrixStats_method2(m2, nam2, "rowMaxs")
    .test_NaArray_matrixStats_method2(m2, nam2, "colRanges")
    .test_NaArray_matrixStats_method2(m2, nam2, "rowRanges")
    .test_NaArray_matrixStats_method2(m2, nam2, "colSums")
    .test_NaArray_matrixStats_method2(m2, nam2, "rowSums")
    .test_NaArray_matrixStats_method2(m2, nam2, "colProds")
    #.test_NaArray_matrixStats_method2(m2, nam2, "rowProds")
    .test_NaArray_matrixStats_method2(m2, nam2, "colMeans")
    #.test_NaArray_matrixStats_method2(m2, nam2, "rowMeans")
    .test_NaArray_matrixStats_method2(m2, nam2, "colSums2")
    .test_NaArray_matrixStats_method2(m2, nam2, "rowSums2")
    .test_NaArray_matrixStats_method2(m2, nam2, "colMeans2")
    #.test_NaArray_matrixStats_method2(m2, nam2, "rowMeans2")
    .test_NaArray_matrixStats_method2(m2, nam2, "colVars")
    #.test_NaArray_matrixStats_method2(m2, nam2, "rowVars")
    .test_NaArray_matrixStats_method2(m2, nam2, "colSds")
    #.test_NaArray_matrixStats_method2(m2, nam2, "rowSds")
    m0 <- m2[0, ]
    nam0 <- nam2[0, ]
    expected <- rep(NA_integer_, 5L)
    expect_warning(colMins(nam0), "NAs introduced")
    expect_warning(colMaxs(nam0), "NAs introduced")
    expect_warning(colRanges(nam0), "NAs introduced")
    expect_identical(suppressWarnings(colMins(nam0, useNames=FALSE)), expected)
    expect_identical(suppressWarnings(colMaxs(nam0, useNames=FALSE)), expected)
    expect_identical(suppressWarnings(colRanges(nam0, useNames=FALSE)),
                     cbind(expected, expected, deparse.level=0))
    expect_identical(suppressWarnings(colMins(nam0)),
                     setNames(expected, colnames(m0)))
    expect_identical(suppressWarnings(colMaxs(nam0)),
                     setNames(expected, colnames(m0)))
    expect_identical(suppressWarnings(colRanges(nam0)),
                     cbind(setNames(expected, colnames(m0)), expected,
                           deparse.level=0))
    expect_identical(unname(rowMins(nam0)), rowMins(m0))
    expect_identical(unname(rowMaxs(nam0)), rowMaxs(m0))
    expect_identical(rowRanges(nam0), rowRanges(m0))
})

test_that("matrixStats methods for 3D NaArray objects", {
    ## input of type() "double"
    a <- array(NA, 6:4,
               dimnames=list(letters[1:6], letters[22:26], LETTERS[1:4]))
    a[1, , 2] <- c(1e12, -1234.55, -2.1, -1, -0.55)
    a[3, , 2] <- c(-0.55, 0, 1e-10, 0.88, 1)
    a[5, , 2] <- c(pi, 10.33, 3.4567895e8, 300, 2009.01)
    a[6, 3:4, 2] <- c(0, NaN)
    naa3 <- as(a, "NaArray")

    test_3D_colrowMinsMaxs(naa3)

    ## dims == 1 (default)
    .test_NaArray_matrixStats_method2(a, naa3, "colSums")
    .test_NaArray_matrixStats_method2(a, naa3, "rowSums")
    .test_NaArray_matrixStats_method2(a, naa3, "colMeans")
    #.test_NaArray_matrixStats_method2(a, naa3, "rowMeans")

    ## dims == 2
    .test_NaArray_matrixStats_method2(a, naa3, "colSums", dims=2)
    .test_NaArray_matrixStats_method2(a, naa3, "rowSums", dims=2)
    .test_NaArray_matrixStats_method2(a, naa3, "colMeans", dims=2)
    #.test_NaArray_matrixStats_method2(a, naa3, "rowMeans", dims=2)
})

test_that("more torturing of the *Mins()/*Maxs() methods for NaArray", {

    ## --- 2D objects ---

    m1 <- rbind(c(0L, -8L, NA), c(NA, NA, 1L))
    m2 <- rbind(c(NA, 0L, NA, NA), c(8L, 9L, 1L, 1L), -(8:11))
    for (m in list(m1, m2)) {
        nam <- NaArray(m)
        expect_identical(rowMins(nam), rowMins(m))
        expect_identical(rowMaxs(nam), rowMaxs(m))
        expect_identical(colMins(nam), colMins(m))
        expect_identical(colMaxs(nam), colMaxs(m))
        expect_identical(rowMins(nam, na.rm=TRUE), rowMins(m, na.rm=TRUE))
        expect_identical(rowMaxs(nam, na.rm=TRUE), rowMaxs(m, na.rm=TRUE))
        expect_identical(colMins(nam, na.rm=TRUE), colMins(m, na.rm=TRUE))
        expect_identical(colMaxs(nam, na.rm=TRUE), colMaxs(m, na.rm=TRUE))
    }

    ## --- 3D objects ---

    ## input of type() "integer"
    naa0 <- NaArray(dim=5:3,
                dimnames=list(letters[1:5], letters[23:26], LETTERS[1:3]),
                type="integer")
    naa0[c(1, 6, 16, 20:22, 36, 39:40, 60)] <-
                c(2L, -5L, 0L, 0L, -11L, 99L, -8L, 0L, 0L, 0L)

    suppressWarnings(test_3D_colrowMinsMaxs(naa0))
    expect_warning(rowMins(naa0, na.rm=TRUE, dims=2), "NAs introduced")
    expect_warning(rowMaxs(naa0, na.rm=TRUE, dims=2), "NAs introduced")

    naa <- naa0[ , , 0]
    suppressWarnings(test_3D_colrowMinsMaxs(naa))
    expect_warning(rowMins(naa), "NAs introduced")
    expect_warning(rowMaxs(naa), "NAs introduced")
    expect_warning(rowMins(naa, dims=2), "NAs introduced")
    expect_warning(rowMaxs(naa, dims=2), "NAs introduced")

    naa <- naa0[ , 0, ]
    suppressWarnings(test_3D_colrowMinsMaxs(naa))
    expect_warning(rowMins(naa), "NAs introduced")
    expect_warning(rowMaxs(naa), "NAs introduced")
    expect_warning(colMins(naa, dims=2), "NAs introduced")
    expect_warning(colMaxs(naa, dims=2), "NAs introduced")

    naa <- naa0[ 0, , ]
    suppressWarnings(test_3D_colrowMinsMaxs(naa))
    expect_warning(colMins(naa), "NAs introduced")
    expect_warning(colMaxs(naa), "NAs introduced")
    expect_warning(colMins(naa, dims=2), "NAs introduced")
    expect_warning(colMaxs(naa, dims=2), "NAs introduced")

    ## input of type() "double"
    naa0[39:40] <- c(NaN, NaN)
    test_3D_colrowMinsMaxs(naa0)
    test_3D_colrowMinsMaxs(naa0[ , , 0])
    test_3D_colrowMinsMaxs(naa0[ , 0, ])
    test_3D_colrowMinsMaxs(naa0[ 0, , ])
})

