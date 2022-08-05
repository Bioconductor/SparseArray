#ifndef _SPARSEARRAY_SUBSETTING_H_
#define _SPARSEARRAY_SUBSETTING_H_

#include <Rdefines.h>

SEXP C_drop_SVT_SparseArray_ineffective_dims(
	SEXP x_dim,
	SEXP x_dimnames,
	SEXP x_type,
	SEXP x_SVT
);

SEXP C_subset_SVT_SparseArray(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP index
);

#endif  /* _SPARSEARRAY_SUBSETTING_H_ */

