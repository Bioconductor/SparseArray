#ifndef _SPARSEARRAY_SUBASSIGNMENT_OLD_H_
#define _SPARSEARRAY_SUBASSIGNMENT_OLD_H_

#include <Rdefines.h>

SEXP C_subassign_SVT_by_Lindex_OLD(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP Lindex,
	SEXP vals
);

SEXP C_subassign_SVT_by_Mindex_OLD(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP Mindex,
	SEXP vals
);

#endif  /* _SPARSEARRAY_SUBASSIGNMENT_OLD_H_ */

