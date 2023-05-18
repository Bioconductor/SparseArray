/****************************************************************************
 *                     Summarization of a "leaf vector"                     *
 ****************************************************************************/
#include "leaf_vector_summarization.h"

#include "Rvector_summarization.h"
#include "leaf_vector_utils.h"


int _summarize_leaf_vector(SEXP lv, int d,
		const SummarizeOp *summarize_op, SummarizeResult *res)
{
	int lv_len, bailout;
	SEXP lv_offs, lv_vals;

	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	res->totalcount += d - lv_len;
	res->nzcount += lv_len;  /* assuming 'lv_vals' contains no zeros! */
	bailout = _summarize_Rvector(lv_vals, summarize_op, res);
	if (bailout)
		return 1;
	if (lv_len == d || summarize_op->opcode == SUM_SHIFTED_X2_OPCODE)
		return 0;
	/* TODO: Reconsider the need to call _summarize_one_zero() here.
	   Does it still make sense? Is this the right place for this?
	   Is it correct? Do we need to use this trick in _summarize_SVT()
	   too? (see file SparseArray_summarization.c) */
	return _summarize_one_zero(summarize_op, res);
}

