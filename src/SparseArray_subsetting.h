#ifndef _SPARSEARRAY_SUBSETTING_H_
#define _SPARSEARRAY_SUBSETTING_H_

#include <Rdefines.h>

SEXP C_subset_SVT_SparseArray(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP index
);

#endif  /* _SPARSEARRAY_SUBSETTING_H_ */

