#ifndef _ARGCHECK_UTILS_H_
#define _ARGCHECK_UTILS_H_

#include <Rdefines.h>

SEXPTYPE _get_and_check_Rtype_from_Rstring(
	SEXP type,
	const char *fun,
	const char *argname
);

int _get_and_check_na_background(
	SEXP na_background,
	const char *fun,
	const char *argname
);

void _check_array_conformability(
	SEXP x_dim,
	SEXP y_dim
);

#endif  /* _ARGCHECK_UTILS_H_ */

