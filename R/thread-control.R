### =========================================================================
### OpenMP thread control
### -------------------------------------------------------------------------
###


.normarg_nthread <- function(nthread)
{
    if (!isSingleNumber(nthread))
        stop(wmsg("'nthread' must be a single number"))
    if (!is.integer(nthread))
        nthread <- as.integer(nthread)
    nthread
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### R wrappers to OpenMP thread control functions
###

### Wrapper to omp_get_num_procs().
### Returns 0 if OpenMP is not available (e.g. on macOS).
.get_num_procs <- function()
    .Call2("C_get_num_procs", PACKAGE="SparseArray")

### Wrapper to omp_get_max_threads().
### Default is controlled by environment variable OMP_NUM_THREADS.
### Returns 0 if OpenMP is not available (e.g. on macOS).
.get_max_threads <- function()
    .Call2("C_get_max_threads", PACKAGE="SparseArray")

### Wrapper to omp_set_num_threads().
### No-op if OpenMP is not available (e.g. on macOS).
### Returns previous omp_get_max_threads() value.
.set_max_threads <- function(nthread)
{
    nthread <- .normarg_nthread(nthread)
    .Call2("C_set_max_threads", nthread, PACKAGE="SparseArray")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get/set SparseArray option "nthread"
###

.default_SparseArray_nthread <- function()
{
    nthread <- .get_max_threads()
    if (nthread == 0L)
        return(nthread)
    n <- .get_num_procs() %/% 3L
    if (nthread > n)
        nthread <- n
    if (nthread == 0L)
        nthread <- 1L
    nthread
}

get_SparseArray_nthread <- function()
{
    default <- .default_SparseArray_nthread()
    nthread <- get_SparseArray_option("nthread", default=default)
    if (!isSingleNumber(nthread) || nthread < 0L)
        warning(wmsg("invalid 'getOption(\"SparseArray\")$nthread'"))
    nthread
}

set_SparseArray_nthread <- function(nthread=NULL)
{
    if (.get_max_threads() == 0L) {
        nthread <- 0L
    } else if (is.null(nthread)) {
        nthread <- .default_SparseArray_nthread()
    } else {
        nthread <- .normarg_nthread(nthread)
        if (nthread < 1L)
            stop(wmsg("'nthread' must be >= 1"))
    }
    set_SparseArray_option("nthread", nthread)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SparseArray.Call()
###

SparseArray.Call <- function(.NAME, ...)
{
    prev_max_threads <- .set_max_threads(get_SparseArray_nthread())
    on.exit(.set_max_threads(prev_max_threads))
    .Call2(.NAME, ..., PACKAGE="SparseArray")
}

