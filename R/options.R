### =========================================================================
### Handle SparseArray options
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###

### Must return a named list.
.get_SparseArray_options <- function()
{
    SparseArray_options <- getOption("SparseArray")
    if (is.null(SparseArray_options))
        return(setNames(list(), character(0)))
    if (!is.list(SparseArray_options) || is.null(names(SparseArray_options)))
        stop(wmsg("invalid 'getOption(\"SparseArray\")' ",
                  "(should be a named list)"))
    SparseArray_options
}

get_SparseArray_option <- function(name, default=NULL)
{
    stopifnot(isSingleString(name))
    SparseArray_options <- .get_SparseArray_options()
    if (name %in% names(SparseArray_options))
        return(SparseArray_options[[name]])
    default
}

set_SparseArray_option <- function(name, value)
{
    stopifnot(isSingleString(name))
    SparseArray_options <- .get_SparseArray_options()
    prev_value <- SparseArray_options[[name]]
    SparseArray_options[name] <- list(value)
    options(SparseArray=SparseArray_options)
    invisible(prev_value)
}

SparseArray_option_is_set <- function(name)
{
    stopifnot(isSingleString(name))
    SparseArray_options <- .get_SparseArray_options()
    name %in% names(SparseArray_options)
}

