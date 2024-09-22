
test_that("NaArray subsetting by an Mindex or Lindex", {
    for (background in c(NA_integer_, 0L)) {
        ## --- 3D ---
        a0 <- array(background, c(7, 10, 3),
                    dimnames=list(NULL, letters[1:10], LETTERS[1:3]))
        a0[ , 2, 1] <- a0[c(1:2, 6), 4, 1] <- 1L
        a0[ , 8, 1] <- 81:87
        a0[ , -1, 3] <- 308:370
        Mindex2 <- rbind(c(7,  9), c(7, 10), c(6, 4), c(2, 4), c(1, 10),
                         c(7, 10), c(1,  1), c(5, 4), c(2, 4))
        Mindex3 <- rbind(cbind(Mindex2, 1),
                         cbind(Mindex2, 2),
                         cbind(Mindex2, 3))
        naa0 <- as(a0, "NaArray")
        test_linear_subsetting(a0, naa0, Mindex3)

        ## --- 2D ---
        m0 <- a0[ , , 1]
        naa0 <- as(m0, "NaArray")
        test_linear_subsetting(m0, naa0, Mindex2)

        ## --- 1D ---
        x0 <- as.array(m0[1, ])
        Mindex1 <- Mindex2[ , -2, drop=FALSE]
        naa0 <- as(x0, "NaArray")
        test_linear_subsetting(x0, naa0, Mindex1)
    }
})

test_that("NaArray subsetting by an Nindex", {
    a1 <- a2 <- array(NA_integer_, c(7, 10, 3),
                      dimnames=list(NULL, letters[1:10], LETTERS[1:3]))
    a2[ , 2, 1] <- a2[c(1:2, 6), 4, 1] <- 1L
    a2[ , 8, 1] <- 81:87
    a2[ , -1, 3] <- 308:370
    a3 <- array(0L, dim(a1), dimnames=dimnames(a1))
    a4 <- a2
    a4[is.na(a2)] <- 0L

    for (a0 in list(a1, a2, a3, a4)) {
        ## --- 3D ---
        naa0 <- as(a0, "NaArray")

        expect_identical(naa0[ , , ], naa0)

        a   <- a0  [ , c(4:3, 8), 1, drop=FALSE]
        naa <- naa0[ , c(4:3, 8), 1, drop=FALSE]
        check_NaArray_object(naa, a)
        m   <- a0  [ , c(4:3, 8), 1, drop=TRUE]
        naa <- naa0[ , c(4:3, 8), 1, drop=TRUE]
        check_NaArray_object(naa, m)

        a   <- a0  [7, , , drop=FALSE]
        naa <- naa0[7, , , drop=FALSE]
        check_NaArray_object(naa, a)
        m   <- a0  [7, , , drop=TRUE]
        naa <- naa0[7, , , drop=TRUE]
        check_NaArray_object(naa, m)

        ## --- 2D ---
        m0 <- a0[ , , 1]
        naa0 <- as(m0, "NaMatrix")

        expect_identical(naa0[ , ], naa0)

        m   <- m0  [-5 , c(4:3, 8)]
        naa <- naa0[-5 , c(4:3, 8)]
        check_NaArray_object(naa, m)

        expect_identical(naa0[ , 4], m0[ , 4])
        m   <- m0  [ , 4, drop=FALSE]
        naa <- naa0[ , 4, drop=FALSE]
        check_NaArray_object(naa, m)

        expect_identical(naa0[6 , ], m0[6 , ])
        m   <- m0  [6, -1, drop=FALSE]
        naa <- naa0[6, -1, drop=FALSE]
        check_NaArray_object(naa, m)

        ## --- 1D ---
        x0 <- as.array(m0[6, ])
        naa0 <- as(x0, "NaArray")

        expect_identical(naa0[ ], naa0)

        x   <- x0  [c(8:4, 1, 4), drop=FALSE]
        naa <- naa0[c(8:4, 1, 4), drop=FALSE]
        check_NaArray_object(naa, x)
        x   <- x0  [-4, drop=FALSE]
        naa <- naa0[-4, drop=FALSE]
        check_NaArray_object(naa, x)
        x   <- x0  [-4, drop=FALSE]
        naa <- naa0[-4, drop=FALSE]
        check_NaArray_object(naa, x)

        subscript <- c("d", "j", "j", "h")
        x   <- x0  [subscript, drop=FALSE]  # 'drop=TRUE' would do the same!
        naa <- naa0[subscript, drop=FALSE]
        check_NaArray_object(naa, x)
        expect_identical(naa0[subscript], S4Arrays:::drop_even_if_1D(x))
    }
})

