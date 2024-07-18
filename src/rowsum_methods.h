#ifndef _ROWSUM_METHODS_H_
#define _ROWSUM_METHODS_H_

#include <Rdefines.h>

SEXP C_rowsum_SVT(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP group,
	SEXP ngroup,
	SEXP na_rm
);

SEXP C_rowsum_dgCMatrix(
	SEXP x,
	SEXP group,
	SEXP ngroup,
	SEXP na_rm
);

SEXP C_colsum_SVT(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP group,
	SEXP ngroup,
	SEXP na_rm
);

SEXP C_colsum_dgCMatrix(
	SEXP x,
	SEXP group,
	SEXP ngroup,
	SEXP na_rm
);

#endif  /* _ROWSUM_METHODS_H_ */

