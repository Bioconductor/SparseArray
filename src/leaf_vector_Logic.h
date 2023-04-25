#ifndef _LEAF_VECTOR_LOGIC_H_
#define	_LEAF_VECTOR_LOGIC_H_ 

#include <Rdefines.h>

/* Operations from 'Logic' group */
#define	AND_OPCODE	1  /* "&" */
#define	OR_OPCODE	2  /* "|" */

int _get_Logic_opcode(SEXP op);

#endif  /* _LEAF_VECTOR_LOGIC_H_ */

