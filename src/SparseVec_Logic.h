#ifndef _SPARSEVEC_LOGIC_H_
#define _SPARSEVEC_LOGIC_H_

#include <Rdefines.h>

#include "SparseVec.h"

/* Operations from 'Logic' group */
#define	AND_OPCODE	1  /* "&" */
#define	OR_OPCODE	2  /* "|" */

int _get_Logic_opcode(SEXP op);

int _Logic_intSV_intSV(
	int opcode,
	const SparseVec *sv1,
	const SparseVec *sv2,
	int *out_nzvals,
	int *out_nzoffs
);

#endif  /* _SPARSEVEC_LOGIC_H_ */

