/****************************************************************************
 *               matrixStats methods for SparseMatrix objects               *
 ****************************************************************************/
#include "matrixStats_methods.h"

#include "Rvector_utils.h"
#include "Rvector_summarization.h"
#include "SparseArray_summarization.h"

#include <string.h>  /* for memcpy() */


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

/* Recursive. */
static void REC_colStats2_SVT(SEXP SVT, const int *dims, int ndim,
		const SummarizeOp *summarize_op,
		double *out, const long long int *out_incs, int out_ndim)
{
	SummarizeResult res;
	int SVT_len, i;
	long long int out_inc;
	SEXP subSVT;

	if (out_ndim == 0) {
		res = _summarize_SVT(SVT, dims, ndim, summarize_op);
		*out =  res.outbuf.one_double[0];
		return;
	}
	SVT_len = dims[ndim - 1];
	out_inc = out_incs[out_ndim - 1];
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		REC_colStats2_SVT(subSVT, dims, ndim - 1,
				  summarize_op,
				  out, out_incs, out_ndim - 1);
		out += out_inc;
	}
	return;
}

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

/* --- .Call ENTRY POINT ---
   Operations using "interface 2": FUN(x, na.rm, dims) */
SEXP C_colStats2_SVT(SEXP x_dim, SEXP x_dimnames, SEXP x_type, SEXP x_SVT,
		     SEXP op, SEXP na_rm, SEXP dims)
{
	SEXPTYPE Rtype;
	int opcode, narm, ans_ndim;
	SummarizeOp summarize_op;
	SEXP ans_dim, ans;
	long long int *out_incs;

	Rtype = _get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("SparseArray internal error in "
		      "C_colStats2_SVT():\n"
		      "    SVT_SparseArray object has invalid type");

	opcode = _get_summarize_opcode(op, Rtype);

	if (!(IS_LOGICAL(na_rm) && LENGTH(na_rm) == 1))
		error("'na.rm' must be TRUE or FALSE");
	narm = LOGICAL(na_rm)[0];

	summarize_op = _make_SummarizeOp(opcode, Rtype, narm, 0.0);

	ans_dim = PROTECT(get_ans_dim(x_dim, dims));
	ans_ndim = LENGTH(ans_dim);
	out_incs = (long long int *) R_alloc(ans_ndim, sizeof(long long int));
	ans = PROTECT(alloc_ans(REALSXP, ans_dim, out_incs));
	REC_colStats2_SVT(x_SVT, INTEGER(x_dim), LENGTH(x_dim),
			  &summarize_op,
			  REAL(ans), out_incs, ans_ndim);
	UNPROTECT(2);
	return ans;
}

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


