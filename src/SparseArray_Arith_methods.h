#ifndef _SPARSEARRAY_ARITH_METHODS_H_
#define _SPARSEARRAY_ARITH_METHODS_H_

#include <Rdefines.h>

SEXP C_unary_minus_SVT(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT
);

SEXP C_Arith_SVT1_v2(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP x_na_background,
	SEXP v2,
	SEXP op,
	SEXP ans_type
);

SEXP C_Arith_v1_SVT2(
	SEXP v1,
	SEXP y_dim,
	SEXP y_type,
	SEXP y_SVT,
	SEXP y_na_background,
	SEXP op,
	SEXP ans_type
);

SEXP C_Arith_SVT1_SVT2(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP x_na_background,
	SEXP y_dim,
	SEXP y_type,
	SEXP y_SVT,
	SEXP y_na_background,
	SEXP op,
	SEXP ans_type
);

#endif  /* _SPARSEARRAY_ARITH_METHODS_H_ */

