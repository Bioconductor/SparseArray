#ifndef _SPARSEMATRIX_MULT_H_
#define _SPARSEMATRIX_MULT_H_

#include <Rdefines.h>

SEXP C_SVT_crossprod1(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP ans_type,
	SEXP ans_dimnames
);

SEXP C_SVT_crossprod2(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP y_dim,
	SEXP y_type,
	SEXP y_SVT,
	SEXP ans_type,
	SEXP ans_dimnames
);

#endif  /* _SPARSEMATRIX_MULT_H_ */

