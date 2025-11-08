#ifndef _SPARSEVEC_SUBASSIGNMENT_H_
#define _SPARSEVEC_SUBASSIGNMENT_H_

#include <Rdefines.h>

#include "SparseVec.h"

int _subassign_SV_with_Rsubvec(
	const SparseVec *sv,
	const int *offs,
	int n,
	SEXP Rvector,
	R_xlen_t subvec_offset,
	SparseVec *out_sv
);

int _subassign_full_SV_with_Rsubvec(
	const SparseVec *sv,
	SEXP Rvector,
	R_xlen_t subvec_offset,
	SparseVec *out_sv
);

int _subassign_SV_with_Rvector_selection(
	const SparseVec *sv,
	const int *offs,
	int n,
	SEXP Rvector,
	const int *selection,
	SparseVec *out_sv
);

#endif  /* _SPARSEVEC_SUBASSIGNMENT_H_ */

