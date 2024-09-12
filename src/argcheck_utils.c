/****************************************************************************
 *               Some utility functions for argument checking               *
 ****************************************************************************/
#include "argcheck_utils.h"

#include "Rvector_utils.h"  /* _get_Rtype_from_Rstring() */

#include <string.h>  /* for memcmp() */


SEXPTYPE _get_and_check_Rtype_from_Rstring(SEXP type,
		const char *fun, const char *argname)
{
	SEXPTYPE Rtype = _get_Rtype_from_Rstring(type);
	if (Rtype == 0)
		error("SparseArray internal error in %s():\n"
		      "    invalid '%s' value", fun, argname);
	return Rtype;
}

int _get_and_check_na_background(SEXP na_background,
		const char *fun, const char *argname)
{
	if (!(IS_LOGICAL(na_background) && LENGTH(na_background) == 1))
		error("SparseArray internal error in %s():\n"
		      "    '%s' must be TRUE or FALSE", fun, argname);
	return LOGICAL(na_background)[0];
}

void _check_array_conformability(SEXP x_dim, SEXP y_dim)
{
	int ndim = LENGTH(x_dim);
	if (ndim != LENGTH(y_dim) ||
	    memcmp(INTEGER(x_dim), INTEGER(y_dim), sizeof(int) * ndim) != 0)
		error("non-conformable arrays");
	return;
}

