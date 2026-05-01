.STRINGS01 <- c("", "1")

.onLoad <- function(libname, pkgname)
{
    if (!SparseArray_option_is_set("nthread"))
        set_SparseArray_nthread()
    .Call2("C_init_character0_character1", .STRINGS01, PACKAGE="SparseArray")
}

.onUnload <- function(libpath)
{
    library.dynam.unload("SparseArray", libpath)
}

