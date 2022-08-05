#ifndef _RANDOM_SPARSEARRAY_H_
#define	_RANDOM_SPARSEARRAY_H_

#include <Rdefines.h>

SEXP C_simple_rpois(
	SEXP n,
	SEXP lambda
);

SEXP C_poissonSparseArray(
	SEXP dim,
	SEXP lambda
);

#endif  /* _RANDOM_SPARSEARRAY_H_ */

