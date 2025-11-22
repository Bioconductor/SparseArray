#ifndef _SPARSEVEC_SUBSETTING_H_
#define _SPARSEVEC_SUBSETTING_H_

#include <Rdefines.h>

#include "SparseVec.h"

void _subset_SV(
	const SparseVec *sv,
	const int *offs,
	SparseVec *out_sv,
	int *lookup_table
);

#endif  /* _SPARSEVEC_SUBSETTING_H_ */

