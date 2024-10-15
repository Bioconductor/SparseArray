
test_that("nzwhich(), nzvals(), `nzvals<-`() on COO_SparseArray objects", {

    dim <- c(5L, 8L, 2L)
    nzcoo  <- Lindex2Mindex(c(5:12, 15, 3:1), dim)
    nzdata <- c(101:107, 0, 109:110, 0, 112)
    coo0 <- COO_SparseArray(dim, nzcoo, nzdata)
    a0 <- as.array(coo0)

    check_array_like_object(coo0, "COO_SparseArray", a0)
    expect_identical(nzwhich(coo0), nzwhich(a0))
    expect_identical(nzvals(coo0), nzvals(a0))

    coo <- coo0
    a <- a0

    nzvals(coo) <- nzvals(a) <- nzvals(a) + 0.99
    check_array_like_object(coo, "COO_SparseArray", a)
    expect_identical(nzwhich(coo0), nzwhich(a0))
    expect_identical(nzvals(coo0), nzvals(a0))

    nzvals(coo) <- nzvals(a) <- -0.1 * (11:15)
    check_array_like_object(coo, "COO_SparseArray", a)
    expect_identical(nzwhich(coo0), nzwhich(a0))
    expect_identical(nzvals(coo0), nzvals(a0))
})

test_that("COO_SparseMatrix <==> [C|R|T]sparseMatrix coercions", {
    CsparseMatrix <- SparseArray:::CsparseMatrix
    RsparseMatrix <- SparseArray:::RsparseMatrix
    TsparseMatrix <- SparseArray:::TsparseMatrix
    normalize_COO_SparseArray <- SparseArray:::.normalize_COO_SparseArray

    ## Construct a COO_SparseMatrix object 'coo1' with an unnormalized 'nzcoo'
    ## slot and a normalized 'nzdata' slot:
    nzcoo1 <- rbind(c(4L,4L), c(4L,3L), c(3L,3L), c(1L,2L), c(3L,1L),
                    c(4L,1L), c(5L,1L), c(1L,2L), c(1L,4L))
    nzdata1 <- 11:19
    coo1 <- COO_SparseArray(5:4, nzcoo=nzcoo1, nzdata=nzdata1,
                            dimnames=list(letters[1:5], LETTERS[1:4]))

    ## COO_SparseMatrix <==> dgCMatrix.
    dgc1 <- CsparseMatrix(dim(coo1),
                          i=nzcoo1[ , 1], j=nzcoo1[ , 2], nzdata=nzdata1,
                          dimnames=dimnames(coo1))
    expect_identical(as(coo1, "dgCMatrix"), dgc1)
    expect_identical(as(coo1, "CsparseMatrix"), dgc1)
    target <- normalize_COO_SparseArray(`type<-`(coo1, type(dgc1)))
    expect_identical(as(dgc1, "COO_SparseMatrix"), target)
    expect_identical(as(dgc1, "COO_SparseArray"), target)
    expect_identical(as(target, "dgCMatrix"), dgc1)

    ## COO_SparseMatrix <==> dgRMatrix.
    dgr1 <- RsparseMatrix(dim(coo1),
                          i=nzcoo1[ , 1], j=nzcoo1[ , 2], nzdata=nzdata1,
                          dimnames=dimnames(coo1))
    expect_identical(as(coo1, "dgRMatrix"), dgr1)
    expect_identical(as(coo1, "RsparseMatrix"), dgr1)
    target <- t(normalize_COO_SparseArray(t(`type<-`(coo1, type(dgr1)))))
    expect_identical(as(dgr1, "COO_SparseMatrix"), target)
    expect_identical(as(dgr1, "COO_SparseArray"), target)
    expect_identical(as(dgr1, "SparseMatrix"), target)
    expect_identical(as(dgr1, "SparseArray"), target)
    expect_identical(as(target, "dgRMatrix"), dgr1)

    ## COO_SparseMatrix <==> dgTMatrix.
    dgt1 <- TsparseMatrix(dim(coo1),
                          i=nzcoo1[ , 1], j=nzcoo1[ , 2], nzdata=nzdata1,
                          dimnames=dimnames(coo1))
    expect_identical(as(coo1, "dgTMatrix"), dgt1)
    expect_identical(as(coo1, "TsparseMatrix"), dgt1)
    expect_identical(as(coo1, "sparseMatrix"), dgt1)
    target <- normalize_COO_SparseArray(`type<-`(coo1, type(dgc1)))
    expect_identical(as(dgt1, "COO_SparseMatrix"), target)
    expect_identical(as(dgt1, "COO_SparseArray"), target)

    dgt1b <- TsparseMatrix(dim(coo1),
                           i=nzcoo1[ , 1], j=nzcoo1[ , 2], nzdata=nzdata1,
                           dimnames=dimnames(coo1), drop.dups=FALSE)
    target <- as(as.matrix(dgt1b), "COO_SparseMatrix")
    expect_identical(as(dgt1b, "COO_SparseMatrix"), target)
    expect_identical(as(dgt1b, "COO_SparseArray"), target)

    ## COO_SparseMatrix <==> ngCMatrix.
    ngc1 <- CsparseMatrix(dim(coo1),
                          i=nzcoo1[ , 1], j=nzcoo1[ , 2],
                          dimnames=dimnames(coo1))
    expect_identical(as(coo1, "ngCMatrix"), ngc1)
    target <- as(as.matrix(ngc1), "COO_SparseMatrix")
    expect_identical(as(ngc1, "COO_SparseMatrix"), target)
    expect_identical(as(ngc1, "COO_SparseArray"), target)
    expect_identical(as(target, "ngCMatrix"), ngc1)

    ## COO_SparseMatrix <==> ngRMatrix.
    ngr1 <- RsparseMatrix(dim(coo1),
                          i=nzcoo1[ , 1], j=nzcoo1[ , 2],
                          dimnames=dimnames(coo1))
    expect_identical(as(coo1, "ngRMatrix"), ngr1)
    target <- as(as.matrix(ngr1), "COO_SparseMatrix")
    target <- t(normalize_COO_SparseArray(t(target)))
    expect_identical(as(ngr1, "COO_SparseMatrix"), target)
    expect_identical(as(ngr1, "COO_SparseArray"), target)
    expect_identical(as(target, "ngRMatrix"), ngr1)

    ## COO_SparseMatrix <==> ngTMatrix.
    ngt1 <- TsparseMatrix(dim(coo1),
                          i=nzcoo1[ , 1], j=nzcoo1[ , 2],
                          dimnames=dimnames(coo1))
    expect_identical(as(coo1, "ngTMatrix"), ngt1)
    target <- as(as.matrix(ngt1), "COO_SparseMatrix")
    expect_identical(as(ngt1, "COO_SparseMatrix"), target)
    expect_identical(as(ngt1, "COO_SparseArray"), target)
    object <- as(target, "ngTMatrix")
    expect_true(is(object, "ngTMatrix"))
    expect_true(validObject(object))
    expect_true(S4Vectors:::sortedIntegerPairs(object@j, object@i,
                                               strictly=TRUE))

    ngt1b <- TsparseMatrix(dim(coo1),
                           i=nzcoo1[ , 1], j=nzcoo1[ , 2],
                           dimnames=dimnames(coo1), drop.dups=FALSE)
    target <- as(as.matrix(ngt1b), "COO_SparseMatrix")
    expect_identical(as(ngt1b, "COO_SparseMatrix"), target)
    expect_identical(as(ngt1b, "COO_SparseArray"), target)
})

