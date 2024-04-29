#ifndef _SPARSE_VEC_LOGIC_H_
#define _SPARSE_VEC_LOGIC_H_

#include <Rdefines.h>

#include "sparse_vec.h"

/* Operations from 'Logic' group */
#define	AND_OPCODE	1  /* "&" */
#define	OR_OPCODE	2  /* "|" */

int _get_Logic_opcode(SEXP op);

int _sparse_vec_Logic_ints_ints(
	int opcode,
	const struct sparse_vec *sv1,
	const struct sparse_vec *sv2,
	int *out_nzoffs,
	int *out_nzvals
);

#endif  /* _SPARSE_VEC_LOGIC_H_ */

