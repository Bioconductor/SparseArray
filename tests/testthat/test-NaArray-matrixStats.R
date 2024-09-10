
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
    naa1 <- as(m1, "NaArray")
    .test_NaArray_matrixStats_method1(m1, naa1, "colAnyNAs")
    #.test_NaArray_matrixStats_method1(m1, naa1, "rowAnyNAs")
    m1[1, 2] <- NA
    naa1 <- as(m1, "NaArray")
    .test_NaArray_matrixStats_method1(m1, naa1, "colAnyNAs")
    #.test_NaArray_matrixStats_method1(m1, naa1, "rowAnyNAs")

    ## input of type() "logical"
    m2 <- matrix(c(FALSE, FALSE, TRUE,
                   FALSE,  TRUE, TRUE), nrow=2, byrow=TRUE,
                 dimnames=list(LETTERS[1:2], letters[1:3]))
    naa2 <- as(m2, "NaArray")
    .test_NaArray_matrixStats_method1(m2, naa2, "colAnyNAs")
    #.test_NaArray_matrixStats_method1(m2, naa2, "rowAnyNAs")
    m2[1, 2] <- NA
    naa2 <- as(m2, "NaArray")
    .test_NaArray_matrixStats_method1(m2, naa2, "colAnyNAs")
    #.test_NaArray_matrixStats_method1(m2, naa2, "rowAnyNAs")

    ## input of type() "double"
    m3 <- matrix(c(0,    0,  pi,
                   0, 0.25, 1e3), nrow=2, byrow=TRUE,
                 dimnames=list(LETTERS[1:2], letters[1:3]))
    naa3 <- as(m3, "NaArray")
    .test_NaArray_matrixStats_method1(m3, naa3, "colAnyNAs")
    #.test_NaArray_matrixStats_method1(m3, naa3, "rowAnyNAs")
    m3[1, 2] <- NaN
    naa3 <- as(m3, "NaArray")
    .test_NaArray_matrixStats_method1(m3, naa3, "colAnyNAs")
    #.test_NaArray_matrixStats_method1(m3, naa3, "rowAnyNAs")
    m3[1, 2] <- NA
    naa3 <- as(m3, "NaArray")
    .test_NaArray_matrixStats_method1(m3, naa3, "colAnyNAs")
    #.test_NaArray_matrixStats_method1(m3, naa3, "rowAnyNAs")

    ## input of type() "complex"
    m4 <- matrix(c(0,    0,  pi,
                   0, 2-5i, 1e3), nrow=2, byrow=TRUE,
                 dimnames=list(LETTERS[1:2], letters[1:3]))
    naa4 <- as(m4, "NaArray")
    .test_NaArray_matrixStats_method1(m4, naa4, "colAnyNAs")
    #.test_NaArray_matrixStats_method1(m4, naa4, "rowAnyNAs")
    m4[1, 2] <- NaN       # 1st type of "complex" NaN
    naa4 <- as(m4, "NaArray")
    .test_NaArray_matrixStats_method1(m4, naa4, "colAnyNAs")
    #.test_NaArray_matrixStats_method1(m4, naa4, "rowAnyNAs")
    m4[1, 2] <- NaN * 1i  # 2nd type of "complex" NaN
    naa4 <- as(m4, "NaArray")
    .test_NaArray_matrixStats_method1(m4, naa4, "colAnyNAs")
    #.test_NaArray_matrixStats_method1(m4, naa4, "rowAnyNAs")
    m4[1, 2] <- NA
    naa4 <- as(m4, "NaArray")
    .test_NaArray_matrixStats_method1(m4, naa4, "colAnyNAs")
    #.test_NaArray_matrixStats_method1(m4, naa4, "rowAnyNAs")

    ## input of type() "character"
    m5 <- matrix(c("",     "", "Hello",
                   "", "dear", "world"), nrow=2, byrow=TRUE,
                 dimnames=list(LETTERS[1:2], letters[1:3]))
    naa5 <- as(m5, "NaArray")
    .test_NaArray_matrixStats_method1(m5, naa5, "colAnyNAs")
    #.test_NaArray_matrixStats_method1(m5, naa5, "rowAnyNAs")
    m5[1, 2] <- NA
    naa5 <- as(m5, "NaArray")
    .test_NaArray_matrixStats_method1(m5, naa5, "colAnyNAs")
    #.test_NaArray_matrixStats_method1(m5, naa5, "rowAnyNAs")
})

test_that("other matrixStats methods for 2D NaArray objects", {
    ## input of type() "integer"
    m1 <- matrix(c( 0L, 0L,  NA, 0L, NA,
                    NA, 0L, -3L, 1L, NA,
                    0L, 0L,  0L, 0L, 0L,
                   15L, 0L,  0L, 0L, NA), nrow=4, byrow=TRUE,
                 dimnames=list(LETTERS[1:4], letters[1:5]))
    naa1 <- as(m1, "NaArray")
    .test_NaArray_matrixStats_method2(m1, naa1, "colAnys")
    #.test_NaArray_matrixStats_method2(m1, naa1, "rowAnys")
    .test_NaArray_matrixStats_method2(m1, naa1, "colAlls")
    #.test_NaArray_matrixStats_method2(m1, naa1, "rowAlls")
    .test_NaArray_matrixStats_method2(m1, naa1, "colMins")
    #.test_NaArray_matrixStats_method2(m1, naa1, "rowMins")
    .test_NaArray_matrixStats_method2(m1, naa1, "colMaxs")
    #.test_NaArray_matrixStats_method2(m1, naa1, "rowMaxs")
    .test_NaArray_matrixStats_method2(m1, naa1, "colRanges")
    #.test_NaArray_matrixStats_method2(m1, naa1, "rowRanges")
    .test_NaArray_matrixStats_method2(m1, naa1, "colSums")
    .test_NaArray_matrixStats_method2(m1, naa1, "rowSums")
    .test_NaArray_matrixStats_method2(m1, naa1, "colProds")
    #.test_NaArray_matrixStats_method2(m1, naa1, "rowProds")
    .test_NaArray_matrixStats_method2(m1, naa1, "colMeans")
    #.test_NaArray_matrixStats_method2(m1, naa1, "rowMeans")
    .test_NaArray_matrixStats_method2(m1, naa1, "colSums2")
    .test_NaArray_matrixStats_method2(m1, naa1, "rowSums2")
    .test_NaArray_matrixStats_method2(m1, naa1, "colMeans2")
    #.test_NaArray_matrixStats_method2(m1, naa1, "rowMeans2")
    .test_NaArray_matrixStats_method2(m1, naa1, "colVars")
    #.test_NaArray_matrixStats_method2(m1, naa1, "rowVars")
    .test_NaArray_matrixStats_method2(m1, naa1, "colSds")
    #.test_NaArray_matrixStats_method2(m1, naa1, "rowSds")
    m0 <- m1[0, ]
    naa0 <- naa1[0, ]
    expected <- rep(NA_integer_, 5L)
    expect_warning(colMins(naa0), "NAs introduced")
    expect_warning(colMaxs(naa0), "NAs introduced")
    expect_warning(colRanges(naa0), "NAs introduced")
    expect_identical(suppressWarnings(colMins(naa0, useNames=FALSE)), expected)
    expect_identical(suppressWarnings(colMaxs(naa0, useNames=FALSE)), expected)
    expect_identical(suppressWarnings(colRanges(naa0, useNames=FALSE)),
                     cbind(expected, expected, deparse.level=0))
    expect_identical(suppressWarnings(colMins(naa0)),
                     setNames(expected, colnames(m0)))
    expect_identical(suppressWarnings(colMaxs(naa0)),
                     setNames(expected, colnames(m0)))
    expect_identical(suppressWarnings(colRanges(naa0)),
                     cbind(setNames(expected, colnames(m0)), expected,
                           deparse.level=0))
    #expect_identical(unname(rowMins(naa0)), rowMins(m0))
    #expect_identical(unname(rowMaxs(naa0)), rowMaxs(m0))
    #expect_identical(rowRanges(naa0), rowRanges(m0))

    ## input of type() "logical"
    m2 <- is.na(m1)
    naa2 <- as(m2, "NaArray")
    .test_NaArray_matrixStats_method2(m2, naa2, "colAnys")
    #.test_NaArray_matrixStats_method2(m2, naa2, "rowAnys")
    .test_NaArray_matrixStats_method2(m2, naa2, "colAlls")
    #.test_NaArray_matrixStats_method2(m2, naa2, "rowAlls")
    storage.mode(m2) <- "integer"
    .test_NaArray_matrixStats_method2(m2, naa2, "colMins")
    #.test_NaArray_matrixStats_method2(m2, naa2, "rowMins")
    .test_NaArray_matrixStats_method2(m2, naa2, "colMaxs")
    #.test_NaArray_matrixStats_method2(m2, naa2, "rowMaxs")
    .test_NaArray_matrixStats_method2(m2, naa2, "colRanges")
    #.test_NaArray_matrixStats_method2(m2, naa2, "rowRanges")
    .test_NaArray_matrixStats_method2(m2, naa2, "colSums")
    .test_NaArray_matrixStats_method2(m2, naa2, "rowSums")
    .test_NaArray_matrixStats_method2(m2, naa2, "colProds")
    #.test_NaArray_matrixStats_method2(m2, naa2, "rowProds")
    .test_NaArray_matrixStats_method2(m2, naa2, "colMeans")
    #.test_NaArray_matrixStats_method2(m2, naa2, "rowMeans")
    .test_NaArray_matrixStats_method2(m2, naa2, "colSums2")
    .test_NaArray_matrixStats_method2(m2, naa2, "rowSums2")
    .test_NaArray_matrixStats_method2(m2, naa2, "colMeans2")
    #.test_NaArray_matrixStats_method2(m2, naa2, "rowMeans2")
    .test_NaArray_matrixStats_method2(m2, naa2, "colVars")
    #.test_NaArray_matrixStats_method2(m2, naa2, "rowVars")
    .test_NaArray_matrixStats_method2(m2, naa2, "colSds")
    #.test_NaArray_matrixStats_method2(m2, naa2, "rowSds")
    m0 <- m2[0, ]
    naa0 <- naa2[0, ]
    expected <- rep(NA_integer_, 5L)
    expect_warning(colMins(naa0), "NAs introduced")
    expect_warning(colMaxs(naa0), "NAs introduced")
    expect_warning(colRanges(naa0), "NAs introduced")
    expect_identical(suppressWarnings(colMins(naa0, useNames=FALSE)), expected)
    expect_identical(suppressWarnings(colMaxs(naa0, useNames=FALSE)), expected)
    expect_identical(suppressWarnings(colRanges(naa0, useNames=FALSE)),
                     cbind(expected, expected, deparse.level=0))
    expect_identical(suppressWarnings(colMins(naa0)),
                     setNames(expected, colnames(m0)))
    expect_identical(suppressWarnings(colMaxs(naa0)),
                     setNames(expected, colnames(m0)))
    expect_identical(suppressWarnings(colRanges(naa0)),
                     cbind(setNames(expected, colnames(m0)), expected,
                           deparse.level=0))
    #expect_identical(unname(rowMins(naa0)), rowMins(m0))
    #expect_identical(unname(rowMaxs(naa0)), rowMaxs(m0))
    #expect_identical(rowRanges(naa0), rowRanges(m0))
})

test_that("matrixStats methods for 3D NaArray objects", {
    ## input of type() "double"
    a <- array(0, 6:4,
               dimnames=list(letters[1:6], letters[22:26], LETTERS[1:4]))
    a[1, , 2] <- c(1e12, -1234.55, -2.1, -1, -0.55)
    a[3, , 2] <- c(-0.55, 0, 1e-10, 0.88, 1)
    a[5, , 2] <- c(pi, 10.33, 3.4567895e8, 300, 2009.01)
    a[6, 3:4, 2] <- c(NA, NaN)
    naa3 <- as(a, "NaArray")

    ## dims == 1 (default)
    expected <- apply(a, MARGIN=3, colMins, useNames=TRUE)
    expect_identical(colMins(naa3), expected)
    expected <- apply(a, MARGIN=1, min)
    #expect_identical(rowMins(naa3), expected)
    expected <- apply(a, MARGIN=3, colMaxs, useNames=TRUE)
    expect_identical(colMaxs(naa3), expected)
    expected <- apply(a, MARGIN=1, max)
    #expect_identical(rowMaxs(naa3), expected)
    .test_NaArray_matrixStats_method2(a, naa3, "colSums")
    .test_NaArray_matrixStats_method2(a, naa3, "rowSums")
    .test_NaArray_matrixStats_method2(a, naa3, "colMeans")
    #.test_NaArray_matrixStats_method2(a, naa3, "rowMeans")

    ## dims == 2
    expected <- apply(a, MARGIN=3, min)
    expect_identical(colMins(naa3, dims=2), expected)
    expected <- apply(a, MARGIN=2, rowMins, useNames=TRUE)
    #expect_identical(rowMins(naa3, dims=2), expected)
    expected <- apply(a, MARGIN=3, max)
    expect_identical(colMaxs(naa3, dims=2), expected)
    expected <- apply(a, MARGIN=2, rowMaxs, useNames=TRUE)
    #expect_identical(rowMaxs(naa3, dims=2), expected)
    .test_NaArray_matrixStats_method2(a, naa3, "colSums", dims=2)
    .test_NaArray_matrixStats_method2(a, naa3, "rowSums", dims=2)
    .test_NaArray_matrixStats_method2(a, naa3, "colMeans", dims=2)
    #.test_NaArray_matrixStats_method2(a, naa3, "rowMeans", dims=2)
})

