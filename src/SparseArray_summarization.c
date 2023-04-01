/****************************************************************************
 *              Summarization methods for SparseArray objects               *
 ****************************************************************************/
#include "SparseArray_summarization.h"

#include "Rvector_summarization.h"
#include "leaf_vector_utils.h"


/****************************************************************************
 * C_summarize_SVT_SparseArray()
 */


/* Recursive. */
static int REC_summarize_SVT(SEXP SVT, const int *dim, int ndim,
		const SummarizeOp *summarize_op,
		void *init, R_xlen_t *na_rm_count, int status,
		int *has_null_leaves)
{
	int SVT_len, i;
	SEXP subSVT;

	if (SVT == R_NilValue) {
		*has_null_leaves = 1;
		return status;
	}

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		return _summarize_leaf_vector(SVT, dim[0],
					summarize_op,
					init, na_rm_count, status);
	}

	/* 'SVT' is a regular node (list). */
	SVT_len = LENGTH(SVT);
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		status = REC_summarize_SVT(subSVT, dim, ndim - 1,
					summarize_op,
					init, na_rm_count, status,
					has_null_leaves);
		if (status == 2)
			break;
	}
	return status;
}

static int summarize_SVT(SEXP SVT, const int *dim, int ndim,
		int opcode, SEXPTYPE Rtype,
		double *init, int na_rm, R_xlen_t *na_rm_count,
		double shift)
{
	SummarizeOp summarize_op;
	int has_null_leaves, status;

	summarize_op = _init_SummarizeOp(opcode, Rtype, na_rm, shift, init);
	*na_rm_count = 0;

	if (SVT == R_NilValue)
		return 0;
	
	status = has_null_leaves = 0;
	status = REC_summarize_SVT(SVT, dim, ndim,
				   &summarize_op,
				   init, na_rm_count, status,
				   &has_null_leaves);
	if (status == 2 || !has_null_leaves || opcode == SUM_SHIFTED_X2_OPCODE)
		return status;
	if (Rtype == INTSXP) {
		int zero = 0;
		status = _apply_summarize_op(&summarize_op,
					     init, &zero, 1,
					     na_rm_count, status);
	} else {
		double zero = 0.0;
		status = _apply_summarize_op(&summarize_op,
					     init, &zero, 1,
					     na_rm_count, status);
	}
	return status;
}

SEXP C_summarize_SVT_SparseArray(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		SEXP op, SEXP na_rm, SEXP shift)
{
	SEXPTYPE Rtype;
	int opcode, narm0, status;
	R_xlen_t na_rm_count;
	double init[2];  /* 'init' will store 1 or 2 ints or doubles */

	Rtype = _get_Rtype_from_Rstring(x_type);
        if (Rtype == 0)
		error("SparseArray internal error in "
		      "C_summarize_SVT_SparseArray():\n"
		      "    SVT_SparseArray object has invalid type");

	opcode = _get_summarize_opcode(op, Rtype);

	if (!(IS_LOGICAL(na_rm) && LENGTH(na_rm) == 1))
		error("'na.rm' must be TRUE or FALSE");
	narm0 = LOGICAL(na_rm)[0];

	if (!IS_NUMERIC(shift) || LENGTH(shift) != 1)
		error("SparseArray internal error in "
		      "C_summarize_SVT_SparseArray():\n"
		      "    'shift' must be a single numeric value");

	status = summarize_SVT(x_SVT, INTEGER(x_dim), LENGTH(x_dim),
			       opcode, Rtype,
			       init, narm0, &na_rm_count,
			       REAL(shift)[0]);

	return _make_SEXP_from_summarize_result(opcode, Rtype,
			       init, narm0, na_rm_count, status);
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

