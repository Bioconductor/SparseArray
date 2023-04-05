library(Matrix)
library(SparseArray)

summarize_array_like_objects <- function(objects)
{
    stopifnot(is.list(objects))
    object_classes <- sapply(objects, function(x) class(x)[[1L]])
    object_dims <- sapply(objects, function(x) paste(dim(x), collapse="x"))
    setNames(sprintf("<%s %s>", object_dims, object_classes), names(objects))
}

print_object_summaries <- function(object_summaries)
{
    cat("Objects:\n")
    object_names <- format(names(object_summaries))
    for (i in seq_along(object_summaries))
        cat("  ", object_names[[i]], " = ", object_summaries[[i]], "\n", sep="")
}

exec_and_print_scores <- function(FUN, list_of_args)
{
    FUN <- match.fun(FUN)
    stopifnot(is.list(list_of_args))
    lapply(seq_along(list_of_args),
        function(i) {
            args <- list_of_args[[i]]
            stopifnot(is.list(args))
            t <- system.time(ans <- do.call(FUN, args))[["user.self"]]
            if (i == 1L)
                t1 <<- t
            label <- paste0(names(list_of_args)[[i]], ":")
            tm <- sprintf("%6.3fs", t)
            score <- sprintf("%8.4g", t1/t)
            cat("  ", format(label, width=40), " ",
                      tm, " ", score, "\n", sep="")
            ans
        }
    )
}

## --- symetric crossprod ---

benchmark_sym_crossprod <- function(dgcm1)
{
    cat("--- Benchmark symetric crossprod ---\n")

    m1 <- as.matrix(dgcm1)
    svt1 <- as(dgcm1, "SVT_SparseMatrix")
    objects <- list(m1=m1, dgcm1=dgcm1, svt1=svt1)
    object_summaries <- summarize_array_like_objects(objects)
    print_object_summaries(object_summaries)
    cat("Times/Scores (higher score is faster):\n")

    ## --- Unary crossprod() ---

    list_of_args <- lapply(objects, list)  # make list of lists
    names(list_of_args) <- sprintf("crossprod(%s)", names(objects))
    res <- exec_and_print_scores(crossprod, list_of_args)

    ## Sanity checks:
    cp1 <- res[[1L]]
    cp2 <- res[[2L]]
    cp3 <- res[[3L]]
    stopifnot(identical(as.matrix(cp2), cp1))
    stopifnot(identical(cp3, cp1))

    ## --- Binary symetric crossprod() ---

    list_of_args <- lapply(objects, function(x) list(x, x))
    names(list_of_args) <- sprintf("crossprod(%s, %s)",
                                   names(objects), names(objects))
    res <- exec_and_print_scores(crossprod, list_of_args)

    ## Sanity checks:
    stopifnot(identical(res[[1L]], cp1))
    stopifnot(identical(res[[2L]], as(cp2, "generalMatrix")))
    stopifnot(identical(res[[3L]], cp3))
}

## --- asymetric crossprod ---

benchmark_asym_crossprod <- function(dgcm1, dgcm2)
{
    cat("--- Benchmark asymetric crossprod ---\n")

    m1 <- as.matrix(dgcm1)
    svt1 <- as(dgcm1, "SVT_SparseMatrix")
    m2 <- as.matrix(dgcm2)
    svt2 <- as(dgcm2, "SVT_SparseMatrix")
    objects1 <- list(m1=m1, dgcm1=dgcm1, svt1=svt1)
    objects2 <- list(m2=m2, dgcm2=dgcm2, svt2=svt2)
    object_summaries <- summarize_array_like_objects(c(objects1, objects2))
    print_object_summaries(object_summaries)
    cat("Times/Scores (higher score is faster):\n")

    list_of_args <- mapply(function(x1, x2) list(x1, x2),
                           objects1, objects2, SIMPLIFY=FALSE)
    names(list_of_args) <- sprintf("crossprod(%s, %s)",
                                   names(objects1), names(objects2))
    res <- exec_and_print_scores(crossprod, list_of_args)

    ## Sanity checks:
    cp1 <- res[[1L]]
    cp2 <- res[[2L]]
    cp3 <- res[[3L]]
    stopifnot(identical(as.matrix(cp2), cp1))
    stopifnot(identical(cp3, cp1))

    list_of_args <- mapply(function(x1, x2) list(x2, x1),
                           objects1, objects2, SIMPLIFY=FALSE)
    names(list_of_args) <- sprintf("crossprod(%s, %s)",
                                   names(objects2), names(objects1))
    res <- exec_and_print_scores(crossprod, list_of_args)

    ## Sanity checks:
    stopifnot(identical(cp1, t(res[[1L]])))
    stopifnot(identical(cp2, t(res[[2L]])))
    stopifnot(identical(cp3, t(res[[3L]])))
}

# Make big dgCMatrix objects with no NAs/NaNs or infinite values.

make_dgcm1 <- function()
{
    set.seed(333)
    rsparsematrix(25000, 400, density=0.07)
}

make_dgcm2 <- function()
{
    set.seed(333)
    rsparsematrix(25000, 650, density=0.20)
}

dgcm1 <- make_dgcm1()
benchmark_sym_crossprod(dgcm1)
# --- Benchmark symetric crossprod ---
# Objects:
#   m1    = <25000x400 matrix>
#   dgcm1 = <25000x400 dgCMatrix>
#   svt1  = <25000x400 SVT_SparseMatrix>
# Times/Scores (higher score is faster):
#   crossprod(m1):                            2.877s        1
#   crossprod(dgcm1):                         0.218s     13.2
#   crossprod(svt1):                          0.224s    12.84
#   crossprod(m1, m1):                        5.725s        1
#   crossprod(dgcm1, dgcm1):                  0.195s    29.36
#   crossprod(svt1, svt1):                    0.381s    15.03

dgcm2 <- make_dgcm2()
benchmark_asym_crossprod(dgcm1, dgcm2)
# --- Benchmark asymetric crossprod ---
# Objects:
#   m1    = <25000x400 matrix>
#   dgcm1 = <25000x400 dgCMatrix>
#   svt1  = <25000x400 SVT_SparseMatrix>
#   m2    = <25000x650 matrix>
#   dgcm2 = <25000x650 dgCMatrix>
#   svt2  = <25000x650 SVT_SparseMatrix>
# Times/Scores (higher score is faster):
#   crossprod(m1, m2):                        9.362s        1
#   crossprod(dgcm1, dgcm2):                  0.830s    11.28
#   crossprod(svt1, svt2):                    0.641s    14.61
#   crossprod(m2, m1):                        6.287s        1
#   crossprod(dgcm2, dgcm1):                  0.683s    9.205
#   crossprod(svt2, svt1):                    0.608s    10.34

