#ifndef _SPARSEARRAY_SUBSETTING_H_
#define _SPARSEARRAY_SUBSETTING_H_

#include <Rdefines.h>

SEXP C_subset_SVT_by_Lindex(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP x_na_background,
	SEXP Lindex
);

SEXP C_subset_SVT_by_Mindex(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP x_na_background,
	SEXP Mindex
);

SEXP C_subset_SVT_by_Nindex(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP Nindex
);

#endif  /* _SPARSEARRAY_SUBSETTING_H_ */

