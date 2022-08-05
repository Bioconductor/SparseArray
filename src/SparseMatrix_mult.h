#ifndef _SPARSEMATRIX_MULT_H_
#define _SPARSEMATRIX_MULT_H_

#include <Rdefines.h>

SEXP C_SVT_SparseMatrix_crossprod(
	SEXP x_dim,
	SEXP x_SVT,
	SEXP y_dim,
	SEXP y_SVT,
	SEXP ans_type
);

#endif  /* _SPARSEMATRIX_MULT_H_ */

