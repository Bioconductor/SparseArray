#ifndef _LEAF_VECTOR_SUMMARIZATION_H_
#define	_LEAF_VECTOR_SUMMARIZATION_H_

#include <Rdefines.h>
#include "Rvector_summarization.h"

int _summarize_leaf_vector(
	SEXP lv,
	int d,
	const SummarizeOp *summarize_op,
	SummarizeResult *res
);

#endif  /* _LEAF_VECTOR_SUMMARIZATION_H_ */

