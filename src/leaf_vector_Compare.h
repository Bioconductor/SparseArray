#ifndef _LEAF_VECTOR_COMPARE_H_
#define	_LEAF_VECTOR_COMPARE_H_ 

#include <Rdefines.h>

/* Operations from 'Compare' group */
#define	EQ_OPCODE	1  /* "==" */
#define	NE_OPCODE	2  /* "!=" */
#define	LE_OPCODE	3  /* "<=" */
#define	GE_OPCODE	4  /* ">=" */
#define	LT_OPCODE	5  /* "<" */
#define	GT_OPCODE	6  /* ">" */

int _get_Compare_opcode(SEXP op);

SEXP _Compare_lv1_v2(
	SEXP lv1,
	SEXP v2,
	int opcode,
	int *offs_buf,
	int *vals_buf
);

SEXP _Compare_lv1_lv2(
	SEXP lv1,
	SEXP lv2,
	int opcode,
	int *offs_buf,
	int *vals_buf
);

#endif  /* _LEAF_VECTOR_COMPARE_H_ */

