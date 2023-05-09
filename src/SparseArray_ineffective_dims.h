#ifndef _SPARSEARRAY_INEFFECTIVE_DIMS_H_
#define	_SPARSEARRAY_INEFFECTIVE_DIMS_H_ 

#include <Rdefines.h>

SEXP C_select_SVT_dims(
	SEXP x_dim,
	SEXP x_dimnames,
	SEXP x_type,
	SEXP x_SVT,
	SEXP dim_selection
);

SEXP C_drop_SVT_ineffective_dims(
	SEXP x_dim,
	SEXP x_dimnames,
	SEXP x_type,
	SEXP x_SVT
);

#endif  /* _SPARSEARRAY_INEFFECTIVE_DIMS_H_ */

