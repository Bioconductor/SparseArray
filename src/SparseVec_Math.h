#ifndef _SPARSEVEC_MATH_H_
#define _SPARSEVEC_MATH_H_ 

#include <Rdefines.h>

#include "SparseVec.h"

typedef double (*MathFUN)(double);

MathFUN _get_MathFUN(const char *op);

int _Math_doubleSV(
	MathFUN fun,
	const SparseVec *sv,
	double digits,
	double *out_nzvals,
	int *out_nzoffs,
	int *newNaNs
);

#endif  /* _SPARSEVEC_MATH_H_ */

