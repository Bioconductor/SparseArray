#ifndef _SPARSEARRAY_MATRIXSTATS_H_
#define	_SPARSEARRAY_MATRIXSTATS_H_

#include <Rdefines.h>

SEXP C_colStats_SVT(
	SEXP x_dim,
	SEXP x_dimnames,
	SEXP x_type,
	SEXP x_SVT,
	SEXP na_background,
	SEXP op,
	SEXP na_rm,
	SEXP center,
	SEXP dims
);

SEXP C_rowStats_SVT(
	SEXP x_dim,
	SEXP x_dimnames,
	SEXP x_type,
	SEXP x_SVT,
	SEXP na_background,
	SEXP op,
	SEXP na_rm,
	SEXP center,
	SEXP dims
);

#endif  /* _SPARSEARRAY_MATRIXSTATS_H_ */

