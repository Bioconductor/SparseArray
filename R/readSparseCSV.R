### =========================================================================
### readSparseCSV() and writeSparseCSV()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### writeSparseCSV()
###

.write_csv_line <- function(rowname, vals, filepath, sep, write.zeros)
{
    if (write.zeros) {
        line <- as.character(vals)
    } else {
        nzidx <- which_is_nonzero(vals)
        line <- character(length(vals))
        line[nzidx] <- as.character(vals[nzidx])
    }
    cat(rowname, paste0(sep, line), "\n", file=filepath, sep="", append=TRUE)
}

.write_csv_block <- function(rownames, block, filepath, sep, write.zeros)
{
    stopifnot(is.matrix(block))
    ## Write 'block' row by row. Not very efficient!
    for (i in seq_len(nrow(block))) {
        vals <- as.integer(block[i, ])
        .write_csv_line(rownames[[i]], vals, filepath, sep, write.zeros)
    }
}

### TODO: Maybe use a block-processing approach that is RealizationSink-based
### i.e. it would use:
### - A dedicated RealizationSink object for CSV files e.g. CSVRealizationSink
###   (to be implemented) and its required methods. See
###   DelayedArray/R/write_block.R in the DelayedArray package.
### - DelayedArray::BLOCK_write_to_sink()
### - See writeTENxMatrix() in the HDF5Array package for the details.
###   Note that writeHDF5Array() in the same package also uses a
###   RealizationSink-based approach.
### Major difference with writeTENxMatrix(): BLOCK_write_to_sink() must use
### a **row-oriented** grid obtained with rowAutoGrid(..., nrow=chunknrow).
### I'm not sure BLOCK_write_to_sink() will be able to handle the case
### when 'transpose' is TRUE case out-of-the-box though. Might require some
### tweaks, hopefully nothing really complicated.
writeSparseCSV <- function(x, filepath, sep=",", transpose=FALSE,
                              write.zeros=FALSE, chunknrow=250)
{
    ## Check 'x'.
    x_dim <- dim(x)
    if (length(x_dim) != 2L)
        stop("'x' must be a matrix-like object")
    x_rownames <- rownames(x)
    x_colnames <- colnames(x)
    if (is.null(x_rownames) || is.null(x_colnames))
        stop("'x' must have rownames and colnames")
    x_type <- type(x)
    if (!(x_type %in% c("logical", "integer", "double", "raw")))
        stop(wmsg("'x' must be of type \"integer\" (types \"logical\", ",
                  "\"double\", and \"raw\" are also supported via ",
                  "coercion to \"integer\")"))

    ## Check 'filepath', 'sep', 'transpose', and 'write.zeros'.
    if (!isSingleString(filepath))
        stop(wmsg("'filepath' must be a single string"))
    if (!(isSingleString(sep) && nchar(sep) == 1L))
        stop(wmsg("'sep' must be a single character"))
    if (!isTRUEorFALSE(transpose))
        stop(wmsg("'transpose' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(write.zeros))
        stop(wmsg("'write.zeros' must be TRUE or FALSE"))

    x_nrow <- x_dim[[1L]]
    x_ncol <- x_dim[[2L]]

    if (transpose) {
        chunks <- breakInChunks(x_ncol, chunksize=chunknrow)
        cat(paste0(sep, x_rownames), "\n", file=filepath, sep="")
        ## Note that walking on the columns of a SVT_SparseMatrix should be
        ## slightly more efficient than walking on its rows.
        for (chunk_id in seq_along(chunks)) {
            idx <- chunks[[chunk_id]]
            block <- t(as.matrix(x[ , idx, drop=FALSE]))
            .write_csv_block(x_colnames[idx], block, filepath, sep,
                             write.zeros)
        }
    } else {
        chunks <- breakInChunks(x_nrow, chunksize=chunknrow)
        cat(paste0(sep, x_colnames), "\n", file=filepath, sep="")
        for (chunk_id in seq_along(chunks)) {
            idx <- chunks[[chunk_id]]
            block <- as.matrix(x[idx, , drop=FALSE])
            .write_csv_block(x_rownames[idx], block, filepath, sep,
                             write.zeros)
        }
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### readSparseCSV()
###

.scan_first_two_lines <- function(filepath, sep=",")
{
    con <- file(filepath, "r")
    on.exit(close(con))
    line1 <- scan(con, what=character(), sep=sep, nlines=1L, quiet=TRUE)
    line2 <- scan(con, what=character(), sep=sep, nlines=1L, quiet=TRUE)
    list(line1, line2)
}

.looks_like_a_name <- function(x)
    nzchar(x) & is.na(suppressWarnings(as.numeric(x)))

### Guess by looking at the first 2 lines in the file.
### Not ready (this is a mess!)
.guess_dimnames_and_ncol <- function(line1, line2, rownames=NA, colnames=NA)
{
    if (!(is.logical(rownames) && length(rownames) == 1L))
        stop(wmsg("'rownames' must be a single logical value"))
    if (!(is.logical(colnames) && length(colnames) == 1L))
        stop(wmsg("'colnames' must be a single logical value"))
    n1 <- length(line1)
    n2 <- length(line2)
    if (n2 == 0L) {
        if (n1 == 0L)
            stop(wmsg("invalid file: first two lines are empty"))
        if (isTRUE(rownames) && isTRUE(colnames))
            stop(wmsg("file does not seem to contain both rownames and ",
                      "colnames (2nd line is empty)"))
        if (isTRUE(rownames))
            return(list(c(TRUE, FALSE), n1 - 1L))
        if (isTRUE(colnames))
            return(list(c(FALSE, TRUE), n1))
        return(list(c(FALSE, FALSE), n1))
    }
    if (n1 == n2 - 1L) {
        if (isFALSE(rownames) || isFALSE(colnames))
            stop(wmsg("file seems to contain both rownames and colnames ",
                      "(1st line contains one less item than 2nd line)"))
        return(list(c(TRUE, TRUE), n1))
    }
    if (n1 != n2)
        stop(wmsg("invalid file: nb of items in 2nd line is not equal to ",
                  "n2 or to n2-1 where n2 is the nb of items in 1st line"))
    item1_looks_like_a_name <- .looks_like_a_name(line1[[1L]])
    if (!item1_looks_like_a_name) {
        if (is.na(colnames))
            colnames <- any(.looks_like_a_name(line1[-1L]))
        if (is.na(rownames))
            rownames <- .looks_like_a_name(line2[[1L]])
        return(list(c(rownames, colnames), n2))
    }
    ## File contains either rownames or colnames, but not both.
    if (isTRUE(rownames) && isTRUE(colnames))
        stop(wmsg("file does not seem to contain both rownames and colnames ",
                  "(1st item in the file looks like a name)"))
    if (isFALSE(rownames) && isFALSE(colnames))
        stop(wmsg("file seems to contain either rownames or colnames ",
                  "(1st item in the file looks like a name)"))
    if (!is.na(rownames))
        return(c(rownames, !rownames))
    if (!is.na(colnames))
        return(c(!colnames, colnames))
    if (.looks_like_a_name(line2[[1L]]))
        return(c(TRUE, FALSE))
    if (any(.looks_like_a_name(line1[-1L])))
        return(c(FALSE, TRUE))
    ## Maybe file is more likely to have rownames than colnames but who knows,
    ## this is just a random guess.
    c(TRUE, FALSE)
}

.readSparseCSV_as_SVT_SparseMatrix <- function(con, sep, csv_colnames,
                                               transpose=FALSE)
{
    tmpenv <- new.env(parent=emptyenv())
    C_ans <- .Call2("C_readSparseCSV_as_SVT_SparseMatrix",
                    con, sep, transpose, length(csv_colnames), tmpenv,
                    PACKAGE="SparseArray")
    rm(tmpenv)

    ## Construct SVT_SparseMatrix object.
    csv_rownames <- C_ans[[1L]]
    ans_SVT <- C_ans[[2L]]
    if (transpose) {
        ans_rownames <- csv_colnames
        ans_colnames <- csv_rownames
    } else {
        ans_rownames <- csv_rownames
        ans_colnames <- csv_colnames
    }
    ans_dim <- c(length(ans_rownames), length(ans_colnames))
    ans_dimnames <- list(ans_rownames, ans_colnames)
    ans_type <- "integer"
    new_SVT_SparseArray(ans_dim, ans_dimnames, ans_type, ans_SVT, check=FALSE)
}

### Returns an SVT_SparseMatrix object by default.
readSparseCSV <- function(filepath, sep=",", transpose=FALSE)
{
    ## Check 'filepath', 'sep', and 'transpose'.
    if (!isSingleString(filepath))
        stop(wmsg("'filepath' must be a single string"))
    if (!(isSingleString(sep) && nchar(sep) == 1L))
        stop(wmsg("'sep' must be a single character"))
    if (!isTRUEorFALSE(transpose))
        stop(wmsg("'transpose' must be TRUE or FALSE"))

    first_two_lines <- .scan_first_two_lines(filepath, sep=sep)
    line1 <- first_two_lines[[1L]]
    line2 <- first_two_lines[[2L]]
    n1 <- length(line1)
    n2 <- length(line2)
    if (n1 < 2L)
        stop(wmsg("first line in the file must contain ",
                  "at least 2 items (found ", n1, ")"))
    if (n1 != n2)
        stop(wmsg("first two lines in the file must contain ",
                  "the same number of items"))
    #dimnames_and_ncol <- .guess_dimnames_and_ncol(line1, line2,
    #                                              rownames, colnames)
    #rownames <- dimnames_and_ncol[[1L]][1L]
    #colnames <- dimnames_and_ncol[[1L]][2L]
    #ncol <- dimnames_and_ncol[[2L]]

    #filexp <- open_input_files(filepath)[[1L]]
    con <- file(filepath, "r")
    on.exit(close(con))
    .readSparseCSV_as_SVT_SparseMatrix(con, sep, line1[-1L],
                                       transpose=transpose)
}

readSparseTable <- function(...)
{
    .Deprecated("readSparseCSV")
    readSparseCSV(...)
}

### Reproduces the mysterious "segfault from C stack overflow" error that
### we use to get when readSparseCSV() was called in the context of creating
### the vignette with 'R CMD build', but with C code that is much simpler.
### See src/test.c for more information.
### NOT exported.
test <- function() .Call2("C_test", PACKAGE="SparseArray")

