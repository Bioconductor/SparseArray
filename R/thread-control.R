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

### Wrapper to omp_get_initial_device().
### Returns 0 if OpenMP is not available (e.g. on macOS).
.get_initial_device <- function()
    .Call2("C_get_initial_device", PACKAGE="SparseArray")

### Wrapper to omp_pause_resource().
### No-op and returns 0 if OpenMP is not available (e.g. on macOS).
.pause_resource <- function(hard_pause=FALSE, device_num=.get_initial_device())
{
    ret <- .Call2("C_pause_resource", hard_pause, device_num,
                  PACKAGE="SparseArray")
    if (ret != 0L)
        warning(wmsg("omp_pause_resource() failed to relinquish ",
                     "resources on device ", device_num))
    ret
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

### IMPORTANT NOTE: Once we're done, we **must** call .pause_resource()
### to relinquish resources used by OpenMP. This protects us from running
### into a potential "OpenMP/fork deadlock" situation the next time
### OMP-parallelized code is executed. See this post on SO
### https://stackoverflow.com/questions/49049388 for details about
### the "OpenMP/fork deadlock" problem.
### For example, if we don't call .pause_resource(), then OpenMP resources
### are left in a state that is incompatible with parallelization via
### BiocParallel::MulticoreParam(). More precisely, R will hang the next
### time OMP-parallelized code is executed on the workers started by
### BiocParallel::MulticoreParam() (these workers are started with fork()).
### This can be reproduced with the following code:
if (FALSE) {
    library(SparseArray)
    omp_loops <- function(nloop)
        .Call("C_simple_omp_parallel_for_loop", as.integer(nloop),
              PACKAGE="SparseArray")
    omp_loops(30)  # first OMP-parallelized code execution
    library(BiocParallel)
    BPPARAM <- MulticoreParam(2)
    bplapply(1:3, function(i) omp_loops(30), BPPARAM=BPPARAM)  # hangs!
    ## Note that this does not happen if we replace MulticoreParam() with
    ## SerialParam() or if we don't call omp_loops() a first time before
    ## calling it again in the bplapply() loop.
}

SparseArray.Call <- function(.NAME, ...)
{
    prev_max_threads <- .set_max_threads(get_SparseArray_nthread())
    on.exit({.pause_resource(); .set_max_threads(prev_max_threads)})
    .Call2(.NAME, ..., PACKAGE="SparseArray")
}

