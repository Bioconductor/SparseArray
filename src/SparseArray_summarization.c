/****************************************************************************
 *              Summarization methods for SparseArray objects               *
 ****************************************************************************/
#include "SparseArray_summarization.h"

#include "Rvector_utils.h"
#include "Rvector_summarization.h"
#include "leaf_vector_summarization.h"


/****************************************************************************
 * C_summarize_SVT_SparseArray()
 */

/* Recursive. */
static int REC_summarize_SVT(SEXP SVT, const int *dim, int ndim,
		const SummarizeOp *summarize_op, SummarizeResult *res)
{
	R_xlen_t count;
	int along, SVT_len, i, bailout;
	SEXP subSVT;

	if (SVT == R_NilValue) {
		count = 1;
		for (along = 0; along < ndim; along++)
			count *= dim[along];
		res->totalcount += count;
		return 0;
	}

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		return _summarize_leaf_vector(SVT, dim[0],
					      summarize_op, res);
	}

	/* 'SVT' is a regular node (list). */
	SVT_len = LENGTH(SVT);
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		bailout = REC_summarize_SVT(subSVT, dim, ndim - 1,
						    summarize_op, res);
		if (bailout != 0)
			return bailout;
	}
	return 0;
}

SummarizeResult _summarize_SVT(SEXP SVT, const int *dim, int ndim,
			       const SummarizeOp *summarize_op)
{
	SummarizeResult res;
	int bailout;

	_init_SummarizeResult(summarize_op, &res);
	bailout = REC_summarize_SVT(SVT, dim, ndim, summarize_op, &res);
	if (bailout ||
	    res.nzcount == res.totalcount ||
	    summarize_op->opcode == SUM_SHIFTED_X2_OPCODE)
		return res;
	/* TODO: Reconsider the need to call _summarize_one_zero() here.
	   Does it still make sense? Is this the right place for this?
	   Is it correct? Do we also need to use this trick in
	   _summarize_leaf_vector()? (see file leaf_vector_summarization.c) */
	_summarize_one_zero(summarize_op, &res);
	return res;
}

SEXP C_summarize_SVT_SparseArray(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		SEXP op, SEXP na_rm, SEXP shift)
{
	SEXPTYPE Rtype;
	int opcode, narm;
	SummarizeOp summarize_op;
	SummarizeResult res;

	Rtype = _get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("SparseArray internal error in "
		      "C_summarize_SVT_SparseArray():\n"
		      "    SVT_SparseArray object has invalid type");

	opcode = _get_summarize_opcode(op, Rtype);

	if (!(IS_LOGICAL(na_rm) && LENGTH(na_rm) == 1))
		error("'na.rm' must be TRUE or FALSE");
	narm = LOGICAL(na_rm)[0];

	if (!IS_NUMERIC(shift) || LENGTH(shift) != 1)
		error("SparseArray internal error in "
		      "C_summarize_SVT_SparseArray():\n"
		      "    'shift' must be a single numeric value");

	summarize_op = _make_SummarizeOp(opcode, Rtype, narm, REAL(shift)[0]);
	res = _summarize_SVT(x_SVT, INTEGER(x_dim), LENGTH(x_dim),
			     &summarize_op);
	return _make_SEXP_from_summarize_result(&summarize_op, &res);
}


/****************************************************************************
 * C_count_SVT_SparseArray_NAs()
 */

/* Recursive. */
static R_xlen_t REC_count_SVT_NAs(SEXP SVT, int ndim)
{
	R_xlen_t na_count;
	int SVT_len, i;
	SEXP subSVT;

	if (SVT == R_NilValue)
		return 0;

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		return _count_Rvector_NAs(VECTOR_ELT(SVT, 1));
	}

	/* 'SVT' is a regular node (list). */
	na_count = 0;
	SVT_len = LENGTH(SVT);
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		na_count += REC_count_SVT_NAs(subSVT, ndim - 1);
	}
	return na_count;
}

/* --- .Call ENTRY POINT --- */
SEXP C_count_SVT_SparseArray_NAs(SEXP x_dim, SEXP x_type, SEXP x_SVT)
{
	R_xlen_t na_count;

	na_count = REC_count_SVT_NAs(x_SVT, LENGTH(x_dim));
	if (na_count > INT_MAX)
		return ScalarReal((double) na_count);
	return ScalarInteger((int) na_count);
}


/****************************************************************************
 * C_anyNA_SVT_SparseArray()
 */

/* Recursive. */
static int REC_SVT_has_any_NA(SEXP SVT, int ndim)
{
	int SVT_len, i;
	SEXP subSVT;

	if (SVT == R_NilValue)
		return 0;

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		return _Rvector_has_any_NA(VECTOR_ELT(SVT, 1));
	}

	/* 'SVT' is a regular node (list). */
	SVT_len = LENGTH(SVT);
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		if (REC_SVT_has_any_NA(subSVT, ndim - 1))
			return 1;
	}
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_anyNA_SVT_SparseArray(SEXP x_dim, SEXP x_type, SEXP x_SVT)
{
	return ScalarLogical(REC_SVT_has_any_NA(x_SVT, LENGTH(x_dim)));
}

