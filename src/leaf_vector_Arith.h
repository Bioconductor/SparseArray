#ifndef _LEAF_VECTOR_ARITH_H_
#define _LEAF_VECTOR_ARITH_H_

#include <Rdefines.h>

/* Operations from 'Arith' group */
#define	ADD_OPCODE	1  /* "+" */
#define	SUB_OPCODE	2  /* "-" */
#define	MULT_OPCODE	3  /* "*" */
#define	DIV_OPCODE	4  /* "/" */
#define	POW_OPCODE	5  /* "^" */
#define	MOD_OPCODE	6  /* "%%" */
#define	IDIV_OPCODE	7  /* "%/%" */

int _get_Arith_opcode(SEXP op);

SEXP _unary_minus_leaf_vector(
	SEXP lv,
	SEXPTYPE ans_Rtype
);

SEXP _Arith_lv1_v2(
	SEXP lv1,
	SEXP v2,
	int opcode,
	SEXPTYPE ans_Rtype,
	int *offs_buf,
	void *vals_buf,
	int *ovflow
);

SEXP _Arith_lv1_lv2(
	SEXP lv1,
	SEXP lv2,
	int opcode,
	SEXPTYPE ans_Rtype,
	int *offs_buf,
	void *vals_buf,
	int *ovflow
);

#endif  /* _LEAF_VECTOR_ARITH_H_ */

