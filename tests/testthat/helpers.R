make_3D_test_array <- function()
{
    set.seed(123)
    a0 <- array(0.0, c(7, 10, 3),
                dimnames=list(NULL, letters[1:10], LETTERS[1:3]))
    a0[5*(1:26)] <- runif(26, min=-5, max=10)
    a0[2, c(1:4, 7:9), 1] <- c(NA, NaN, Inf, 3e9, 256, -0.999, -1)
    a0
}

check_SparseArray_object <- function(object, expected_class, a0)
{
    expect_true(class(object) == expected_class)
    expect_true(validObject(object))
    expect_identical(dim(object), dim(a0))
    expect_identical(dimnames(object), dimnames(a0))
    expect_identical(type(object), type(a0))
    a <- as.array(object)
    expect_identical(a, a0)
}

