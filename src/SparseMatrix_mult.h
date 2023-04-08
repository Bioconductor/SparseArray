#ifndef _SPARSEMATRIX_MULT_H_
#define _SPARSEMATRIX_MULT_H_

#include <Rdefines.h>

SEXP C_crossprod2_SVT_mat(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP y,
	SEXP transpose_y,
	SEXP ans_type,
	SEXP ans_dimnames
);

SEXP C_crossprod2_mat_SVT(
	SEXP x,
	SEXP y_dim,
	SEXP y_type,
	SEXP y_SVT,
	SEXP transpose_x,
	SEXP ans_type,
	SEXP ans_dimnames
);

SEXP C_crossprod2_SVT_SVT(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP y_dim,
	SEXP y_type,
	SEXP y_SVT,
	SEXP ans_type,
	SEXP ans_dimnames
);

SEXP C_crossprod1_SVT(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP ans_type,
	SEXP ans_dimnames
);

#endif  /* _SPARSEMATRIX_MULT_H_ */

