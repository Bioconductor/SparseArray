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

void _Arith_sv1_scalar(
	int opcode,
	const SparseVec *sv1,
	SEXP scalar,
	SparseVec *out_sv,
	int *ovflow
);

void _Arith_scalar_sv2(
	int opcode,
	SEXP scalar,
	const SparseVec *sv2,
	SparseVec *out_sv,
	int *ovflow
);

void _Arith_sv1_zero(
	int opcode,
	const SparseVec *sv1,
	SEXPTYPE Rtype2,
	SparseVec *out_sv
);

void _Arith_sv1_na(
	int opcode,
	const SparseVec *sv1,
	SEXPTYPE Rtype2,
	SparseVec *out_sv
);

void _Arith_zero_sv2(
	int opcode,
	SEXPTYPE Rtype1,
	const SparseVec *sv2,
	SparseVec *out_sv
);

void _Arith_na_sv2(
	int opcode,
	SEXPTYPE Rtype1,
	const SparseVec *sv2,
	SparseVec *out_sv
);

void _Arith_sv1_sv2(
	int opcode,
	const SparseVec *sv1,
	const SparseVec *sv2,
	SparseVec *out_sv,
	int *ovflow
);

#endif  /* _SPARSEVEC_ARITH_H_ */

