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

static inline int flip_Compare_opcode(int opcode)
{
	switch (opcode) {
	    case EQ_OPCODE: case NE_OPCODE: return opcode;
	    case LE_OPCODE: return GE_OPCODE;
	    case GE_OPCODE: return LE_OPCODE;
	    case LT_OPCODE: return GT_OPCODE;
	    case GT_OPCODE: return LT_OPCODE;
	}
	error("SparseArray internal error in flip_Compare_opcode():\n"
	      "    invalid 'Compare' opcode: %d", opcode);
	return 0;  /* will never reach this */
}

int _get_Compare_opcode(SEXP op);

void _Compare_sv1_zero(
	int opcode,
	const SparseVec *sv1,
	SparseVec *out_sv
);

void _Compare_sv1_scalar(
	int opcode,
	const SparseVec *sv1,
	SEXP scalar,
	SparseVec *out_sv
);

void _Compare_sv1_sv2(
	int opcode,
	const SparseVec *sv1,
	const SparseVec *sv2,
	SparseVec *out_sv
);

#endif  /* _SPARSEVEC_COMPARE_H_ */

