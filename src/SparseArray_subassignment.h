#ifndef _SPARSEARRAY_SUBASSIGNMENT_H_
#define _SPARSEARRAY_SUBASSIGNMENT_H_

#include <Rdefines.h>

SEXP C_subassign_SVT_by_Lindex(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP x_na_background,
	SEXP Lindex,
	SEXP vals
);

SEXP C_subassign_SVT_by_Mindex(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP x_na_background,
	SEXP Mindex,
	SEXP vals
);

SEXP C_subassign_SVT_with_short_Rvector(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP Nindex,
	SEXP Rvector
);

SEXP C_subassign_SVT_with_Rarray(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP x_na_background,
	SEXP Noffs,
	SEXP Rarray
);

SEXP C_subassign_SVT_with_SVT(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP x_na_background,
	SEXP Noffs,
	SEXP y_dim,
	SEXP y_type,
	SEXP y_SVT,
	SEXP y_na_background
);

#endif  /* _SPARSEARRAY_SUBASSIGNMENT_H_ */

