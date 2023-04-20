#ifndef _SPARSEARRAY_MATH_METHODS_H_
#define _SPARSEARRAY_MATH_METHODS_H_

#include <Rdefines.h>

SEXP C_SVT_Math(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP op
);

SEXP C_SVT_Math2(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP op,
	SEXP digits
);

#endif  /* _SPARSEARRAY_MATH_METHODS_H_ */

