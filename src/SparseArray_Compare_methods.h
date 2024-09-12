#ifndef _SPARSEARRAY_COMPARE_METHODS_H_
#define _SPARSEARRAY_COMPARE_METHODS_H_

#include <Rdefines.h>

SEXP C_Compare_SVT1_v2(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP x_na_background,
	SEXP v2,
	SEXP op
);

SEXP C_Compare_SVT1_SVT2(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP x_na_background,
	SEXP y_dim,
	SEXP y_type,
	SEXP y_SVT,
	SEXP y_na_background,
	SEXP op
);

#endif  /* _SPARSEARRAY_COMPARE_METHODS_H_ */

