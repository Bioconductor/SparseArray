#ifndef _SPARSEARRAY_OPS_METHODS_H_
#define _SPARSEARRAY_OPS_METHODS_H_

#include <Rdefines.h>

SEXP C_Arith_SVT_num(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP y,
	SEXP op,
	SEXP ans_type
);

SEXP C_Arith_SVT_SVT(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP y_dim,
	SEXP y_type,
	SEXP y_SVT,
	SEXP op,
	SEXP ans_type
);

SEXP C_Compare_SVT_SVT(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP y_dim,
	SEXP y_type,
	SEXP y_SVT,
	SEXP op
);

SEXP C_Logic_SVT_SVT(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP y_dim,
	SEXP y_type,
	SEXP y_SVT,
	SEXP op
);

#endif  /* _SPARSEARRAY_OPS_METHODS_H_ */

