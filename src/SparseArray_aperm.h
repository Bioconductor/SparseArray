#ifndef _SPARSEARRAY_APERM_H_
#define	_SPARSEARRAY_APERM_H_

#include <Rdefines.h>

SEXP C_transpose_2D_SVT(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT
);

SEXP C_aperm0_SVT(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP perm
);

SEXP C_aperm_SVT(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP perm
);

#endif  /* _SPARSEARRAY_APERM_H_ */

