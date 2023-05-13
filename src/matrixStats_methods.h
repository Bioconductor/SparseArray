#ifndef _MATRIXSTATS_METHODS_H_
#define	_MATRIXSTATS_METHODS_H_

#include <Rdefines.h>

SEXP C_colStats1_SVT(
	SEXP x_dim,
	SEXP x_dimnames,
	SEXP x_type,
	SEXP x_SVT,
	SEXP op,
	SEXP dims
);

SEXP C_colStats2_SVT(
	SEXP x_dim,
	SEXP x_dimnames,
	SEXP x_type,
	SEXP x_SVT,
	SEXP op,
	SEXP na_rm,
	SEXP dims
);

SEXP C_colStats3_SVT(
	SEXP x_dim,
	SEXP x_dimnames,
	SEXP x_type,
	SEXP x_SVT,
	SEXP op,
	SEXP na_rm,
	SEXP center,
	SEXP dims
);

#endif  /* _MATRIXSTATS_METHODS_H_ */

