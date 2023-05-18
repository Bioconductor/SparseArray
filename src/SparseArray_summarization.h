#ifndef _SPARSEARRAY_SUMMARIZATION_H_
#define _SPARSEARRAY_SUMMARIZATION_H_

#include <Rdefines.h>
#include "Rvector_summarization.h"

SummarizeResult _summarize_SVT(
	SEXP SVT,
	const int *dim,
	int ndim,
	const SummarizeOp *summarize_op
);

SEXP C_summarize_SVT_SparseArray(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP op,
	SEXP na_rm,
	SEXP shift
);

SEXP C_count_SVT_SparseArray_NAs(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT
);

SEXP C_anyNA_SVT_SparseArray(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT
);

#endif  /* _SPARSEARRAY_SUMMARIZATION_H_ */

