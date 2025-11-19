#ifndef _SPARSEVEC_SUBASSIGNMENT_H_
#define _SPARSEVEC_SUBASSIGNMENT_H_

#include <Rdefines.h>

#include "SparseVec.h"

int _subassign_SV_with_Rvector_block(
	const SparseVec *sv,
	const int *offs,
	int n,
	SEXP Rvector,
	R_xlen_t block_offset,
	SparseVec *out_sv
);

int _subassign_SV_with_Rvector_subset(
	const SparseVec *sv,
	const int *offs,
	int n,
	SEXP Rvector,
	const int *selection,
	SparseVec *out_sv
);

int _subassign_SV_with_SV(
	const SparseVec *sv1,
	const int *offs,
	const SparseVec *sv2,
	SparseVec *out_sv
);

int _subassign_full_SV_with_Rvector_block(
	const SparseVec *sv,
	SEXP Rvector,
	R_xlen_t block_offset,
	SparseVec *out_sv
);

#endif  /* _SPARSEVEC_SUBASSIGNMENT_H_ */

