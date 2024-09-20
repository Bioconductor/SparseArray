#ifndef _SPARSEVEC_LOGIC_H_
#define _SPARSEVEC_LOGIC_H_

#include <Rdefines.h>

#include "SparseVec.h"

/* Operations from 'Logic' group */
#define	AND_OPCODE	1  /* "&" */
#define	OR_OPCODE	2  /* "|" */

int _get_Logic_opcode(SEXP op);

void _Logic_intSV_na(
	int opcode,
	const SparseVec *sv1,
	SEXPTYPE Rtype2,
	SparseVec *out_sv
);

void _Logic_intSV_intSV(
	int opcode,
	const SparseVec *sv1,
	const SparseVec *sv2,
	SparseVec *out_sv
);

#endif  /* _SPARSEVEC_LOGIC_H_ */

