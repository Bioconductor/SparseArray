#ifndef _LEAF_VECTOR_OPS_H_
#define _LEAF_VECTOR_OPS_H_

#include <Rdefines.h>

/* Operations from 'Arith' group */
#define	ADD_OPCODE	1  /* "+" */
#define	SUB_OPCODE	2  /* "-" */
#define	MULT_OPCODE	3  /* "*" */
#define	DIV_OPCODE	4  /* "/" */
#define	POW_OPCODE	5  /* "^" */
#define	MOD_OPCODE	6  /* "%%" */
#define	IDIV_OPCODE	7  /* "%/%" */

/* Operations from 'Compare' group */
#define	EQ_OPCODE	1  /* "==" */
#define	NE_OPCODE	2  /* "!=" */
#define	LE_OPCODE	3  /* "<=" */
#define	GE_OPCODE	4  /* ">=" */
#define	LT_OPCODE	5  /* "<" */
#define	GT_OPCODE	6  /* ">" */

/* Operations from 'Logic' group */
#define	AND_OPCODE	1  /* "&" */
#define	OR_OPCODE	2  /* "|" */

int _get_Arith_opcode(
	SEXP op,
	SEXPTYPE x_Rtype,
	SEXPTYPE y_Rtype
);

int _get_Compare_opcode(
	SEXP op,
	SEXPTYPE x_Rtype,
	SEXPTYPE y_Rtype
);

int _get_Logic_opcode(
	SEXP op,
	SEXPTYPE x_Rtype,
	SEXPTYPE y_Rtype
);

#endif  /* _LEAF_VECTOR_OPS_H_ */

