#ifndef _SPARSEARRAY_LOGIC_METHODS_H_
#define _SPARSEARRAY_LOGIC_METHODS_H_

#include <Rdefines.h>

SEXP C_logical_neg_NaSVT(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_NaSVT
);

SEXP C_Logic_NaSVT1_na(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_NaSVT,
	SEXP op
);

SEXP C_Logic_SVT1_SVT2(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP x_na_background,
	SEXP y_dim,
	SEXP y_type,
	SEXP y_SVT,
	SEXP y_na_background,
	SEXP op
);

#endif  /* _SPARSEARRAY_LOGIC_METHODS_H_ */

