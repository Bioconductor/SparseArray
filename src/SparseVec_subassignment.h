#ifndef _SPARSEVEC_SUBASSIGNMENT_H_
#define _SPARSEVEC_SUBASSIGNMENT_H_

#include <Rdefines.h>

#include "SparseVec.h"

void _fill_SV_with_vals(
	const void *vals,
	const int *offs,
	int n,
	SparseVec *out_sv
);

void _subassign_SV1_with_v2(
	const SparseVec *sv1,
	const int *offs2,
	const void *vals2,
	int n2,
	SparseVec *out_sv
);

#endif  /* _SPARSEVEC_SUBASSIGNMENT_H_ */

