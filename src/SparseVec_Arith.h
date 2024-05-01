#ifndef _SPARSEVEC_ARITH_H_
#define _SPARSEVEC_ARITH_H_

#include <Rdefines.h>

#include "SparseVec.h"

/* Operations from 'Arith' group */
#define	ADD_OPCODE	1  /* "+" */
#define	SUB_OPCODE	2  /* "-" */
#define	MULT_OPCODE	3  /* "*" */
#define	DIV_OPCODE	4  /* "/" */
#define	POW_OPCODE	5  /* "^" */
#define	MOD_OPCODE	6  /* "%%" */
#define	IDIV_OPCODE	7  /* "%/%" */

int _get_Arith_opcode(SEXP op);

int _Arith_sv1_scalar(
	int opcode,
	const SparseVec *sv1,
	SEXP scalar,
	SEXPTYPE expected_outRtype,
	int *out_nzoffs,
	void *out_nzvals,
	int *ovflow
);

int _mult_SV_zero(
	const SparseVec *sv,
	SEXPTYPE outRtype,
	int *out_nzoffs,
	void *out_nzvals
);

int _Arith_sv1_sv2(
	int opcode,
	const SparseVec *sv1,
	const SparseVec *sv2,
	SEXPTYPE expected_outRtype,
	int *out_nzoffs,
	void *out_nzvals,
	int *ovflow
);

#endif  /* _SPARSEVEC_ARITH_H_ */

