.onLoad <- function(libname, pkgname)
{
    if (!SparseArray_option_is_set("nthread"))
        set_SparseArray_nthread()
}

.onUnload <- function(libpath)
{
    library.dynam.unload("SparseArray", libpath)
}

