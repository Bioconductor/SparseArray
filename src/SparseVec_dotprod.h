#ifndef _SPARSEVEC_DOTPROD_H_
#define _SPARSEVEC_DOTPROD_H_

#include "SparseVec.h"

double _dotprod_doubleSV_doubleSV(
	const SparseVec *sv1,
	const SparseVec *sv2
);

double _dotprod_doubleSV_finite_doubles(
	const SparseVec *sv1,
	const double *x2
);

double _dotprod_doubleSV_doubles(
	const SparseVec *sv1,
	const double *x2
);

double _dotprod_intSV_noNA_ints(
	const SparseVec *sv1,
	const int *x2
);

double _dotprod_intSV_ints(
	const SparseVec *sv1,
	const int *x2
);

double _dotprod_doubles_zero(
	const double *x,
	int x_len
);

double _dotprod_ints_zero(
	const int *x,
	int x_len
);

double _dotprod_doubleSV_zero(const SparseVec *sv);

double _dotprod_intSV_zero(const SparseVec *sv);

#endif  /* _SPARSEVEC_DOTPROD_H_ */

