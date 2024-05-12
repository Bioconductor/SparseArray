#ifndef _SPARSEVEC_COMPARE_H_
#define _SPARSEVEC_COMPARE_H_

#include <Rdefines.h>

#include "SparseVec.h"

/* Operations from 'Compare' group */
#define	EQ_OPCODE	1  /* "==" */
#define	NE_OPCODE	2  /* "!=" */
#define	LE_OPCODE	3  /* "<=" */
#define	GE_OPCODE	4  /* ">=" */
#define	LT_OPCODE	5  /* "<" */
#define	GT_OPCODE	6  /* ">" */

/* Special value returned by the _Compare_sv1_zero() or _Compare_sv1_scalar()
   functions below to indicate that the result of the 'Commpare' operation
   is a logical sparse vector where all the nonzero values are TRUEs and the
   corresponding offsets are the same as the input ones ('sv1->nzoffs').
   IMPORTANT: If this is the case then the functions don't even write anything
   to output buffers 'out_nzvals' or 'out_nzoffs' so the caller should ignore
   them. */
#define	COMPARE_IS_NOOP -1  /* make sure to use a **negative** int */

static inline int flip_opcode(int opcode)
{
	switch (opcode) {
	    case NE_OPCODE: return opcode;
	    case LT_OPCODE: return GT_OPCODE;
	    case GT_OPCODE: return LT_OPCODE;
	}
	error("SparseArray internal error in flip_opcode():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

int _get_Compare_opcode(SEXP op);

int _Compare_sv1_zero(
	int opcode,
	const SparseVec *sv1,
	int *out_nzvals,
	int *out_nzoffs
);

int _Compare_sv1_scalar(
	int opcode,
	const SparseVec *sv1,
	SEXP scalar,
	int *out_nzvals,
	int *out_nzoffs
);

int _Compare_sv1_sv2(
	int opcode,
	const SparseVec *sv1,
	const SparseVec *sv2,
	int *out_nzvals,
	int *out_nzoffs
);

#endif  /* _SPARSEVEC_COMPARE_H_ */

