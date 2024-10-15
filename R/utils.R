
### A silly trick used only to trigger an error when the function is called
### with no arguments.
check_unused_arguments <- function() NULL

vector_of_zeros <- function(mode="logical", length=0L)
{
    vector(mode=mode, length=length)
}

vector_of_ones <- function(mode="logical", length=0L)
{
    as.fun <- base::get(paste0("as.", mode), envir=asNamespace("base"),
                        mode="function")
    rep.int(as.fun(1L), length)
}

coercion_can_introduce_zeros <- function(from_type, to_type)
{
    if (!isSingleString(from_type))
        stop(wmsg("'from_type' must be a single string"))
    if (!isSingleString(to_type))
        stop(wmsg("'to_type' must be a single string"))
    if (!(to_type %in% c("double", "logical")))
        stop(wmsg("'to_type' must be \"double\" or \"logical\""))
    .Call2("C_coercion_can_introduce_zeros", from_type, to_type,
                                             PACKAGE="SparseArray")
}

coercion_can_introduce_NAs <- function(from_type, to_type)
{
    if (!isSingleString(from_type))
        stop(wmsg("'from_type' must be a single string"))
    if (!isSingleString(to_type))
        stop(wmsg("'to_type' must be a single string"))
    .Call2("C_coercion_can_introduce_NAs", from_type, to_type,
                                           PACKAGE="SparseArray")
}

lacunar_mode_is_on <- function()
    .Call2("C_lacunar_mode_is_on", PACKAGE="SparseArray")

