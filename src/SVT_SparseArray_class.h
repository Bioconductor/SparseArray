#ifndef _SVT_SPARSEARRAY_CLASS_H_
#define _SVT_SPARSEARRAY_CLASS_H_

#include <Rdefines.h>

SEXP C_set_SVT_SparseArray_type(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP new_type
);

SEXP _coerce_SVT(
	SEXP SVT,
	const int *dim,
	int ndim,
	SEXPTYPE from_Rtype,
	SEXPTYPE to_Rtype,
	int *offs_buf
);

R_xlen_t _REC_nzcount_SVT(
	SEXP SVT,
	int ndim
);

SEXP C_nzcount_SVT_SparseArray(
	SEXP x_dim,
	SEXP x_SVT
);

SEXP C_nzwhich_SVT_SparseArray(
	SEXP x_dim,
	SEXP x_SVT,
	SEXP arr_ind
);

SEXP C_from_SVT_SparseArray_to_Rarray(
	SEXP x_dim,
	SEXP x_dimnames,
	SEXP x_type,
	SEXP x_SVT,
	SEXP na_background
);

SEXP C_build_SVT_from_Rarray(
	SEXP x,
	SEXP ans_type
);

SEXP C_from_SVT_SparseMatrix_to_CsparseMatrix(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP as_ngCMatrix
);

SEXP C_build_SVT_from_CSC(
	SEXP dim,
	SEXP indptr,
	SEXP data,
	SEXP indices,
	SEXP indices_are_1based
);

SEXP C_build_SVT_from_CsparseMatrix(
	SEXP x,
	SEXP ans_type
);

SEXP C_from_SVT_SparseArray_to_COO_SparseArray(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT
);

#endif  /* _SVT_SPARSEARRAY_CLASS_H_ */

