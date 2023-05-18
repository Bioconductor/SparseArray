/****************************************************************************
 *               matrixStats methods for SparseMatrix objects               *
 ****************************************************************************/
#include "matrixStats_methods.h"

#include "Rvector_utils.h"
#include "Rvector_summarization.h"
#include "SparseArray_summarization.h"

#include <string.h>  /* for memcpy() */


static SEXPTYPE get_ans_Rtype(const SummarizeOp *summarize_op)
{
	SummarizeResult res;

	_init_SummarizeResult(summarize_op, &res);
	return res.out_Rtype;
}

/* Return 'tail(dim(x), n=-d)'. */
static SEXP get_ans_dim(SEXP x_dim, SEXP dims)
{
	int d, x_ndim, ans_ndim;
	SEXP ans_dim;

	if (TYPEOF(dims) != INTSXP || LENGTH(dims) != 1)
		error("'dims' must be a single integer");
	d = INTEGER(dims)[0];
	x_ndim = LENGTH(x_dim);
	if (d == NA_INTEGER || d < 1 || d > x_ndim)
		error("'dims' must be a single integer "
		      ">= 1 and <= length(dim(x))");
	ans_ndim = x_ndim - d;
	ans_dim = PROTECT(NEW_INTEGER(ans_ndim));
	memcpy(INTEGER(ans_dim), INTEGER(x_dim) + d, sizeof(int) * ans_ndim);
	UNPROTECT(1);
	return ans_dim;
}

static SEXP alloc_ans(SEXPTYPE Rtype, SEXP ans_dim, long long int *out_incs)
{
	int ans_ndim, ans_len, along;
	SEXP ans;
	long long int out_inc;

	ans_ndim = LENGTH(ans_dim);
	if (ans_ndim <= 1) {
		ans_len = ans_ndim == 1 ? INTEGER(ans_dim)[0] : 1;
		ans = PROTECT(allocVector(Rtype, ans_len));
	} else {
		ans = PROTECT(allocArray(Rtype, ans_dim));
	}
	out_inc = 1;
	for (along = 0; along < ans_ndim; along++) {
		out_incs[along] = out_inc;
		out_inc *= INTEGER(ans_dim)[along];
	}
	UNPROTECT(1);
	return ans;
}

static inline void *increment_out(const void *out, SEXPTYPE out_Rtype,
				  long long int inc)
{
	switch (out_Rtype) {
	    case LGLSXP: case INTSXP: return ((int    *) out) + inc;
	    case REALSXP:             return ((double *) out) + inc;
	}
	error("SparseArray internal error in increment_out():\n",
	      "    output type \"%s\" is not supported", type2char(out_Rtype));
	return NULL;  /* will never reach this */
}

/* Only the first value in 'res->outbuf' gets copied. */
static inline void copy_result_to_out(const SummarizeResult *res,
		void *out, SEXPTYPE out_Rtype)
{
	if (out_Rtype != res->out_Rtype)
		error("SparseArray internal error in "
		      "copy_result_to_out():\n"
		      "    out_Rtype != res->out_Rtype");
	switch (out_Rtype) {
	    case LGLSXP: case INTSXP:
		*((int *) out) = res->outbuf.one_int[0];
		return;
	    case REALSXP:
		*((double *) out) = res->outbuf.one_double[0];
		return;
	}
	error("SparseArray internal error in copy_result_to_out():\n",
	      "    output type \"%s\" is not supported", type2char(out_Rtype));
	return;  /* will never reach this */
}


/****************************************************************************
 * C_colStats1_SVT()
 */

/* --- .Call ENTRY POINT ---
   Operations using "interface 1": FUN(x, dims) */
SEXP C_colStats1_SVT(SEXP x_dim, SEXP x_dimnames, SEXP x_type, SEXP x_SVT,
		     SEXP op, SEXP dims)
{
	//int opcode;

	error("C_colStats1_SVT() is not ready yet, sorry!");
	//opcode = _get_matrixStats_opcode(op);
	return R_NilValue;
}


/****************************************************************************
 * C_colStats2_SVT()
 */

/* Recursive. */
static void REC_colStats2_SVT(SEXP SVT, const int *dims, int ndim,
		const SummarizeOp *summarize_op,
		void *out, SEXPTYPE out_Rtype,
		const long long int *out_incs, int out_ndim,
		int *warn)
{
	SummarizeResult res;
	int SVT_len, i;
	long long int out_inc;
	SEXP subSVT;

	if (out_ndim == 0) {
		res = _summarize_SVT(SVT, dims, ndim, summarize_op);
		if (res.warn)
			*warn = 1;
		copy_result_to_out(&res, out, out_Rtype);
		return;
	}
	SVT_len = dims[ndim - 1];
	out_inc = out_incs[out_ndim - 1];
	subSVT = R_NilValue;
	for (i = 0; i < SVT_len; i++) {
		if (SVT != R_NilValue)
			subSVT = VECTOR_ELT(SVT, i);
		REC_colStats2_SVT(subSVT, dims, ndim - 1,
				  summarize_op,
				  out, out_Rtype, out_incs, out_ndim - 1,
				  warn);
		out = increment_out(out, out_Rtype, out_inc);
	}
	return;
}


/* --- .Call ENTRY POINT ---
   Operations using "interface 2": FUN(x, na.rm, dims) */
SEXP C_colStats2_SVT(SEXP x_dim, SEXP x_dimnames, SEXP x_type, SEXP x_SVT,
		     SEXP op, SEXP na_rm, SEXP dims)
{
	SEXPTYPE x_Rtype, ans_Rtype;
	int opcode, narm, ans_ndim, warn;
	SummarizeOp summarize_op;
	SEXP ans_dim, ans;
	long long int *out_incs;

	x_Rtype = _get_Rtype_from_Rstring(x_type);
	if (x_Rtype == 0)
		error("SparseArray internal error in "
		      "C_colStats2_SVT():\n"
		      "    SVT_SparseArray object has invalid type");

	opcode = _get_summarize_opcode(op, x_Rtype);

	if (!(IS_LOGICAL(na_rm) && LENGTH(na_rm) == 1))
		error("'na.rm' must be TRUE or FALSE");
	narm = LOGICAL(na_rm)[0];

	summarize_op = _make_SummarizeOp(opcode, x_Rtype, narm, 0.0);
	ans_Rtype = get_ans_Rtype(&summarize_op);

	ans_dim = PROTECT(get_ans_dim(x_dim, dims));
	ans_ndim = LENGTH(ans_dim);
	out_incs = (long long int *) R_alloc(ans_ndim, sizeof(long long int));
	ans = PROTECT(alloc_ans(ans_Rtype, ans_dim, out_incs));
	warn = 0;
	REC_colStats2_SVT(x_SVT, INTEGER(x_dim), LENGTH(x_dim),
			  &summarize_op,
			  DATAPTR(ans), ans_Rtype, out_incs, ans_ndim,
			  &warn);
	if (warn)
		warning("NAs introduced by coercion of "
			"infinite values to integer range");
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * C_colStats3_SVT()
 */

/* --- .Call ENTRY POINT ---
   Operations using "interface 3": FUN(x, na.rm, center, dims) */
SEXP C_colStats3_SVT(SEXP x_dim, SEXP x_dimnames, SEXP x_type, SEXP x_SVT,
		     SEXP op, SEXP na_rm, SEXP center, SEXP dims)
{
	//int opcode;

	error("C_colStats3_SVT() is not ready yet, sorry!");
	//opcode = _get_matrixStats_opcode(op);
	return R_NilValue;
}

