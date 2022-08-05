#ifndef _SVT_SPARSEARRAY_CLASS_H_
#define _SVT_SPARSEARRAY_CLASS_H_

#include <Rdefines.h>

SEXPTYPE _get_Rtype_from_Rstring(SEXP type);

SEXP C_get_SVT_SparseArray_nzcount(
	SEXP x_dim,
	SEXP x_SVT
);

SEXP C_set_SVT_SparseArray_type(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP new_type
);

SEXP C_from_SVT_SparseArray_to_Rarray(
	SEXP x_dim,
	SEXP x_dimnames,
	SEXP x_type,
	SEXP x_SVT
);

SEXP C_build_SVT_from_Rarray(
	SEXP x,
	SEXP ans_type
);

SEXP C_from_SVT_SparseMatrix_to_CsparseMatrix(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT
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

SEXP C_transpose_SVT_SparseMatrix(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT
);

#endif  /* _SVT_SPARSEARRAY_CLASS_H_ */

