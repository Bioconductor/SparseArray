/****************************************************************************
 *              Summarization methods for SparseArray objects               *
 ****************************************************************************/
#include "SparseArray_summarization.h"

#include "Rvector_utils.h"
#include "Rvector_summarization.h"
#include "leaf_vector_summarization.h"


/****************************************************************************
 * C_summarize_SVT()
 */

/* Recursive. */
static void REC_summarize_SVT(SEXP SVT, const int *dim, int ndim,
		const SummarizeOp *summarize_op, SummarizeResult *res)
{
	R_xlen_t in_len;
	int along, SVT_len, i;
	SEXP subSVT;

	if (SVT == R_NilValue) {
		in_len = 1;
		for (along = 0; along < ndim; along++)
			in_len *= dim[along];
		res->in_length += in_len;
		return;
	}

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		_summarize_leaf_vector(SVT, dim[0], summarize_op, res);
		return;
	}

	/* 'SVT' is a regular node (list). */
	SVT_len = LENGTH(SVT);
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		REC_summarize_SVT(subSVT, dim, ndim - 1, summarize_op, res);
		if (res->outbuf_status == OUTBUF_IS_SET_WITH_BREAKING_VALUE)
			return;  /* Bail out early. */
	}
	return;
}

static SummarizeOp replace_SummarizeOp_center_with_mean(
		SEXP SVT, const int *dim, int ndim,
		const SummarizeOp *summarize_op)
{
	SummarizeOp tmp_op;
	SummarizeResult res;
	double SVT_mean;

	/* Compute 'mean(SVT)'. */
	tmp_op = *summarize_op;
	tmp_op.opcode = MEAN_OPCODE;
	_init_SummarizeResult(&tmp_op, &res);
	REC_summarize_SVT(SVT, dim, ndim, &tmp_op, &res);
	_postprocess_SummarizeResult(&tmp_op, &res);
	SVT_mean = res.outbuf.one_double[0];

	/* Set 'center' with 'mean(SVT)'. */
	tmp_op = *summarize_op;
	tmp_op.center = SVT_mean;
	return tmp_op;
}

SummarizeResult _summarize_SVT(SEXP SVT, const int *dim, int ndim,
			       const SummarizeOp *summarize_op)
{
	SummarizeOp tmp_op;
	SummarizeResult res;

	if (ISNAN(summarize_op->center) &&
		(summarize_op->opcode == SUM_CENTERED_X2_OPCODE ||
		 summarize_op->opcode == VAR1_OPCODE ||
		 summarize_op->opcode == SD1_OPCODE))
	{
		tmp_op = replace_SummarizeOp_center_with_mean(SVT, dim, ndim,
							      summarize_op);
		summarize_op = &tmp_op;
	}

	_init_SummarizeResult(summarize_op, &res);
	REC_summarize_SVT(SVT, dim, ndim, summarize_op, &res);
	_postprocess_SummarizeResult(summarize_op, &res);
	return res;
}

/* --- .Call ENTRY POINT --- */
SEXP C_summarize_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		SEXP op, SEXP na_rm, SEXP center)
{
	SEXPTYPE x_Rtype;
	int opcode, narm;
	SummarizeOp summarize_op;
	SummarizeResult res;

	x_Rtype = _get_Rtype_from_Rstring(x_type);
	if (x_Rtype == 0)
		error("SparseArray internal error in "
		      "C_summarize_SVT():\n"
		      "    SVT_SparseArray object has invalid type");

	opcode = _get_summarize_opcode(op, x_Rtype);

	if (!(IS_LOGICAL(na_rm) && LENGTH(na_rm) == 1))
		error("'na.rm' must be TRUE or FALSE");
	narm = LOGICAL(na_rm)[0];

	if (!IS_NUMERIC(center) || LENGTH(center) != 1)
		error("SparseArray internal error in "
		      "C_summarize_SVT():\n"
		      "    'center' must be a single number");

	summarize_op = _make_SummarizeOp(opcode, x_Rtype, narm,
					 REAL(center)[0]);
	res = _summarize_SVT(x_SVT, INTEGER(x_dim), LENGTH(x_dim),
			     &summarize_op);
	if (res.warn)
		warning("NAs introduced by coercion of "
			"infinite values to integers");
	return _make_SEXP_from_summarize_result(&summarize_op, &res);
}

