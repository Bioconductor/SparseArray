#ifndef _SPARSEARRAY_OPS_METHODS_H_
#define _SPARSEARRAY_OPS_METHODS_H_

#include <Rdefines.h>

SEXP C_SVT_Arith(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP y_dim,
	SEXP y_type,
	SEXP y_SVT,
	SEXP op,
	SEXP ans_type
);

SEXP C_SVT_Compare(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP y_dim,
	SEXP y_type,
	SEXP y_SVT,
	SEXP op
);

SEXP C_SVT_Logic(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP y_dim,
	SEXP y_type,
	SEXP y_SVT,
	SEXP op
);

#endif  /* _SPARSEARRAY_OPS_METHODS_H_ */

