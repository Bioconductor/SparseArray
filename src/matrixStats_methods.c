/****************************************************************************
 ****************************************************************************
 **									   **
 **									   **
 **              matrixStats methods for SparseMatrix objects              **
 **									   **
 **									   **
 ****************************************************************************
 ****************************************************************************/
#include "matrixStats_methods.h"

#include "thread_control.h"  /* for which_max() */
#include "Rvector_utils.h"
#include "Rvector_summarization.h"
#include "SparseVec.h"
#include "leaf_utils.h"
#include "SparseArray_summarization.h"

#include <string.h>  /* for memcpy() */


static SEXPTYPE compute_ans_Rtype(const SummarizeOp *summarize_op)
{
	SummarizeResult res;
	_init_SummarizeResult(summarize_op, &res);
	return res.out_Rtype;
}

static int check_dims(SEXP dims, int min, int max)
{
	if (!IS_INTEGER(dims) || LENGTH(dims) != 1)
		error("'dims' must be a single integer");
	int d = INTEGER(dims)[0];
	if (d == NA_INTEGER || d < min || d > max)
		error("'dims' must be >= %d and <= %d", min, max);
	return d;
}

/* Returns 'tail(dim(x), n=-dims)'. */
static SEXP compute_colStats_ans_dim(SEXP x_dim, SEXP dims)
{
	int x_ndim = LENGTH(x_dim);
	int d = check_dims(dims, 1, x_ndim);
	int ans_ndim = x_ndim - d;
	SEXP ans_dim = PROTECT(NEW_INTEGER(ans_ndim));
	memcpy(INTEGER(ans_dim), INTEGER(x_dim) + d, sizeof(int) * ans_ndim);
	UNPROTECT(1);
	return ans_dim;
}

/* Returns 'head(dim(x), n=dims)'. */
static SEXP compute_rowStats_ans_dim(SEXP x_dim, SEXP dims)
{
	int x_ndim = LENGTH(x_dim);
	int d = check_dims(dims, 1, x_ndim - 1);
	SEXP ans_dim = PROTECT(NEW_INTEGER(d));
	memcpy(INTEGER(ans_dim), INTEGER(x_dim), sizeof(int) * d);
	UNPROTECT(1);
	return ans_dim;
}

/* Returns 'S4Arrays:::simplify_NULL_dimnames(tail(dimnames(x), n=-dims))'. */
static SEXP compute_colStats_ans_dimnames(SEXP x_dimnames, int dims)
{
	if (x_dimnames == R_NilValue)
		return R_NilValue;
	int x_ndim = LENGTH(x_dimnames);
	int any_retained = 0;
	for (int along = dims; along < x_ndim; along++) {
		if (VECTOR_ELT(x_dimnames, along) != R_NilValue) {
			any_retained = 1;
			break;
		}
	}
	if (!any_retained)
		return R_NilValue;
	SEXP ans_dimnames = PROTECT(NEW_LIST(x_ndim - dims));
	for (int along = dims; along < x_ndim; along++)
		SET_VECTOR_ELT(ans_dimnames, along - dims,
			       VECTOR_ELT(x_dimnames, along));
	UNPROTECT(1);
	return ans_dimnames;
}

/* Returns 'S4Arrays:::simplify_NULL_dimnames(head(dimnames(x), n=dims))'. */
static SEXP compute_rowStats_ans_dimnames(SEXP x_dimnames, int dims)
{
	if (x_dimnames == R_NilValue)
		return R_NilValue;
	int any_retained = 0;
	for (int along = 0; along < dims; along++) {
		if (VECTOR_ELT(x_dimnames, along) != R_NilValue) {
			any_retained = 1;
			break;
		}
	}
	if (!any_retained)
		return R_NilValue;
	SEXP ans_dimnames = PROTECT(NEW_LIST(dims));
	for (int along = 0; along < dims; along++)
		SET_VECTOR_ELT(ans_dimnames, along,
			       VECTOR_ELT(x_dimnames, along));
	UNPROTECT(1);
	return ans_dimnames;
}

/* Returns an (uninitialized) atomic vector or an array. */
static SEXP alloc_ans(SEXPTYPE Rtype, SEXP ans_dim, R_xlen_t *out_incs)
{
	SEXP ans;
	int ans_ndim = LENGTH(ans_dim);
	if (ans_ndim <= 1) {
		int ans_len = ans_ndim == 1 ? INTEGER(ans_dim)[0] : 1;
		ans = PROTECT(allocVector(Rtype, ans_len));
	} else {
		ans = PROTECT(allocArray(Rtype, ans_dim));
	}
	R_xlen_t out_inc = 1;
	for (int along = 0; along < ans_ndim; along++) {
		out_incs[along] = out_inc;
		out_inc *= INTEGER(ans_dim)[along];
	}
	UNPROTECT(1);
	return ans;
}

static void propagate_colStats_dimnames(SEXP ans, SEXP x_dimnames, int dims)
{
	if (x_dimnames == R_NilValue)
		return;
	int ans_ndim = LENGTH(x_dimnames) - dims;
	if (ans_ndim == 0)
		return;
	if (ans_ndim == 1) {
		SEXP ans_names = VECTOR_ELT(x_dimnames, dims);
		if (ans_names == R_NilValue)
			return;
		SET_NAMES(ans, ans_names);
		return;
	}
	SEXP ans_dimnames = compute_colStats_ans_dimnames(x_dimnames, dims);
	if (ans_dimnames == R_NilValue)
		return;
	PROTECT(ans_dimnames);
	SET_DIMNAMES(ans, ans_dimnames);
	UNPROTECT(1);
	return;
}

/* 'dims' is guaranteed to be >= 1. */
static void propagate_rowStats_dimnames(SEXP ans, SEXP x_dimnames, int dims)
{
	if (x_dimnames == R_NilValue)
		return;
	if (dims == 1) {
		SEXP ans_names = VECTOR_ELT(x_dimnames, 0);
		if (ans_names == R_NilValue)
			return;
		SET_NAMES(ans, ans_names);
		return;
	}
	SEXP ans_dimnames = compute_rowStats_ans_dimnames(x_dimnames, dims);
	if (ans_dimnames == R_NilValue)
		return;
	PROTECT(ans_dimnames);
	SET_DIMNAMES(ans, ans_dimnames);
	UNPROTECT(1);
	return;
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
	    case INTSXP: case LGLSXP:
		*((int *) out) = res->outbuf.one_int[0];
		return;
	    case REALSXP:
		*((double *) out) = res->outbuf.one_double[0];
		return;
	}
	error("SparseArray internal error in copy_result_to_out():\n"
	      "    output type \"%s\" is not supported", type2char(out_Rtype));
	return;  /* will never reach this */
}


/****************************************************************************
 * C_colStats_SVT()
 */

/* Recursive. */
static void REC_colStats_SVT(SEXP SVT, const int *dims, int ndim,
		const SummarizeOp *summarize_op,
		void *out, SEXPTYPE out_Rtype,
		const R_xlen_t *out_incs, int out_ndim, int pardim,
		int *warn)
{
	if (out_ndim == 0) {
		SummarizeResult res = _summarize_SVT(SVT, dims, ndim,
						     summarize_op);
		if (res.warn)
			*warn = 1;
		copy_result_to_out(&res, out, out_Rtype);
		return;
	}
	int SVT_len = dims[ndim - 1];
	R_xlen_t out_inc = out_incs[out_ndim - 1];
	/* Parallel execution along the biggest dimension only. */
	#pragma omp parallel for schedule(static) if(out_ndim == pardim)
	for (int i = 0; i < SVT_len; i++) {
		SEXP subSVT = SVT == R_NilValue ? R_NilValue
						: VECTOR_ELT(SVT, i);
		void *subout = shift_dataptr(out_Rtype, out, out_inc * i);
		REC_colStats_SVT(subSVT, dims, ndim - 1,
				 summarize_op,
				 subout, out_Rtype,
				 out_incs, out_ndim - 1, pardim,
				 warn);
	}
	return;
}

/* --- .Call ENTRY POINT --- */
SEXP C_colStats_SVT(SEXP x_dim, SEXP x_dimnames, SEXP x_type, SEXP x_SVT,
		    SEXP op, SEXP na_rm, SEXP center, SEXP dims)
{
	SEXPTYPE x_Rtype = _get_Rtype_from_Rstring(x_type);
	if (x_Rtype == 0)
		error("SparseArray internal error in "
		      "C_colStats_SVT():\n"
		      "    SVT_SparseArray object has invalid type");

	int opcode = _get_summarize_opcode(op, x_Rtype);

	if (!(IS_LOGICAL(na_rm) && LENGTH(na_rm) == 1))
		error("'na.rm' must be TRUE or FALSE");
	int narm = LOGICAL(na_rm)[0];

	if (!IS_NUMERIC(center) || LENGTH(center) != 1)
		error("SparseArray internal error in "
		      "C_colStats_SVT():\n"
		      "    'center' must be a single number");

	SummarizeOp summarize_op = _make_SummarizeOp(opcode, x_Rtype, narm,
						     REAL(center)[0]);
	SEXPTYPE ans_Rtype = compute_ans_Rtype(&summarize_op);

	SEXP ans_dim = PROTECT(compute_colStats_ans_dim(x_dim, dims));
	int ans_ndim = LENGTH(ans_dim);  /* = x_ndim - INTEGER(dims)[0] */
	/* Get 1-based rank of biggest dimension. Parallel execution will
	   be along that dimension. */
	int pardim = which_max(INTEGER(ans_dim), ans_ndim) + 1;

	R_xlen_t *out_incs = NULL;
	if (ans_ndim != 0)
		out_incs = (R_xlen_t *) R_alloc(ans_ndim, sizeof(R_xlen_t));

	SEXP ans = PROTECT(alloc_ans(ans_Rtype, ans_dim, out_incs));
	propagate_colStats_dimnames(ans, x_dimnames, INTEGER(dims)[0]);

	int warn = 0;
	REC_colStats_SVT(x_SVT, INTEGER(x_dim), LENGTH(x_dim),
			 &summarize_op,
			 DATAPTR(ans), ans_Rtype,
			 out_incs, ans_ndim, pardim,
			 &warn);
	if (warn)
		warning("NAs introduced by coercion of "
			"infinite values to integers");

	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * C_rowStats_SVT()
 */

static void add_SV_to_doubles(const SparseVec *sv, double *out, int narm)
{
	int nzcount = get_SV_nzcount(sv);
	if (sv->nzvals == NULL) {
		/* lacunar leaf */
		for (int k = 0; k < nzcount; k++)
			out[sv->nzoffs[k]] += double1;
	} else {
		/* regular leaf */
		for (int k = 0; k < nzcount; k++) {
			double x;
			if (get_SV_Rtype(sv) == REALSXP) {
				const double *nzvals_p = sv->nzvals;
				x = nzvals_p[k];
				/* ISNAN(): True for *both* NA and NaN.
				   See <R_ext/Arith.h> */
				if (narm && ISNAN(x))
					continue;
			} else {
				const int *nzvals_p = sv->nzvals;
				int v = nzvals_p[k];
				if (v == NA_INTEGER) {
					if (narm)
						continue;
					x = NA_REAL;
				} else {
					x = (double) v;
				}
			}
			out[sv->nzoffs[k]] += x;
		}
	}
	return;
}

static void rowStats_leaf(SEXP leaf, int dim0,
		const SummarizeOp *summarize_op,
		void *out, SEXPTYPE out_Rtype, int *warn)
{
	SparseVec sv = leaf2SV(leaf, summarize_op->in_Rtype, dim0);
	switch (summarize_op->opcode) {
	    case SUM_OPCODE:
		add_SV_to_doubles(&sv, out, summarize_op->na_rm);
		return;
	}
	error("SparseArray internal error in rowStats_leaf():\n"
	      "    operation not supported");
	return;
}

/* Recursive. */
static void REC_rowStats_SVT(SEXP SVT, const int *dims, int ndim,
		const SummarizeOp *summarize_op,
		void *out, SEXPTYPE out_Rtype,
		const R_xlen_t *out_incs, int out_ndim,
		int *warn)
{
	if (SVT == R_NilValue)
		return;

	if (ndim == 1) { /* 'out_ndim' also guaranteed to be 1 */
		/* 'SVT' is a leaf (i.e. a 1D SVT). */
		rowStats_leaf(SVT, dims[0], summarize_op, out, out_Rtype, warn);
		return;
	}
	int SVT_len = dims[ndim - 1];
	R_xlen_t out_inc = 0;
	if (ndim <= out_ndim) {
		out_ndim--;
		out_inc = out_incs[out_ndim];
	}
	for (int i = 0; i < SVT_len; i++) {
		SEXP subSVT = VECTOR_ELT(SVT, i);
		void *subout = shift_dataptr(out_Rtype, out, out_inc * i);
		REC_rowStats_SVT(subSVT, dims, ndim - 1,
				 summarize_op,
				 subout, out_Rtype,
				 out_incs, out_ndim,
				 warn);
	}
	return;
}

/* --- .Call ENTRY POINT --- */
SEXP C_rowStats_SVT(SEXP x_dim, SEXP x_dimnames, SEXP x_type, SEXP x_SVT,
		    SEXP op, SEXP na_rm, SEXP center, SEXP dims)
{
	SEXPTYPE x_Rtype = _get_Rtype_from_Rstring(x_type);
	if (x_Rtype == 0)
		error("SparseArray internal error in "
		      "C_rowStats_SVT():\n"
		      "    SVT_SparseArray object has invalid type");

	int opcode = _get_summarize_opcode(op, x_Rtype);

	if (!(IS_LOGICAL(na_rm) && LENGTH(na_rm) == 1))
		error("'na.rm' must be TRUE or FALSE");
	int narm = LOGICAL(na_rm)[0];

	if (!IS_NUMERIC(center) || LENGTH(center) != 1)
		error("SparseArray internal error in "
		      "C_rowStats_SVT():\n"
		      "    'center' must be a single number");

	SummarizeOp summarize_op = _make_SummarizeOp(opcode, x_Rtype, narm,
						     REAL(center)[0]);
	SEXPTYPE ans_Rtype = compute_ans_Rtype(&summarize_op);

	SEXP ans_dim = PROTECT(compute_rowStats_ans_dim(x_dim, dims));
	int ans_ndim = LENGTH(ans_dim);  /* = INTEGER(dims)[0] */

	/* 'ans_ndim' is guaranteed to be >= 1. */
	R_xlen_t *out_incs = (R_xlen_t *) R_alloc(ans_ndim, sizeof(R_xlen_t));

	SEXP ans = PROTECT(alloc_ans(ans_Rtype, ans_dim, out_incs));
	propagate_rowStats_dimnames(ans, x_dimnames, ans_ndim);
	_set_Rvector_elts_to_zero(ans);

	int warn = 0;
	REC_rowStats_SVT(x_SVT, INTEGER(x_dim), LENGTH(x_dim),
			 &summarize_op,
			 DATAPTR(ans), ans_Rtype,
			 out_incs, ans_ndim,
			 &warn);
	if (warn)
		warning("NAs introduced by coercion of "
			"infinite values to integers");

	UNPROTECT(2);
	return ans;
}

