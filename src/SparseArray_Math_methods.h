#ifndef _SPARSEARRAY_MATH_METHODS_H_
#define _SPARSEARRAY_MATH_METHODS_H_

#include <Rdefines.h>

SEXP C_Math_SVT(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP x_na_background,
	SEXP op,
	SEXP digits
);

#endif  /* _SPARSEARRAY_MATH_METHODS_H_ */

