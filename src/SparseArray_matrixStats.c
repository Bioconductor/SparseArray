/****************************************************************************
 ****************************************************************************
 **									   **
 **									   **
 **              matrixStats methods for SparseArray objects               **
 **									   **
 **									   **
 ****************************************************************************
 ****************************************************************************/
#include "SparseArray_matrixStats.h"

#include "argcheck_utils.h"
#include "thread_control.h"  /* for which_max() */
#include "Rvector_summarization.h"
#include "SparseVec.h"
#include "leaf_utils.h"
#include "SparseArray_summarization.h"

#include <string.h>  /* for memcpy() and memset() */


/****************************************************************************
 * Some helper functions used by C_colStats_SVT() and/or C_rowStats_SVT()
 */

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
static SEXP compute_colStats_ans_dim(SEXP x_dim, int dims)
{
	int ans_ndim = LENGTH(x_dim) - dims;
	SEXP ans_dim = PROTECT(NEW_INTEGER(ans_ndim));
	memcpy(INTEGER(ans_dim), INTEGER(x_dim) + dims, sizeof(int) * ans_ndim);
	UNPROTECT(1);
	return ans_dim;
}

/* Returns 'head(dim(x), n=ans_ndim)'. */
static SEXP compute_rowStats_ans_dim(SEXP x_dim, int ans_ndim)
{
	SEXP ans_dim = PROTECT(NEW_INTEGER(ans_ndim));
	memcpy(INTEGER(ans_dim), INTEGER(x_dim), sizeof(int) * ans_ndim);
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



/****************************************************************************
 * C_colStats_SVT()                                                         *
 ****************************************************************************/


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

/* Recursive. */
static void REC_colStats_SVT(SEXP SVT, int na_background,
		const int *dim, int ndim,
		const SummarizeOp *summarize_op,
		void *out, SEXPTYPE out_Rtype,
		const R_xlen_t *out_incs, int out_ndim, int pardim,
		int *warn)
{
	if (out_ndim == 0) {
		SummarizeResult res = _summarize_SVT(SVT, na_background,
						     dim, ndim,
						     summarize_op);
		if (res.warn)
			*warn = 1;
		copy_result_to_out(&res, out, out_Rtype);
		return;
	}
	int SVT_len = dim[ndim - 1];
	R_xlen_t out_inc = out_incs[out_ndim - 1];
	/* Parallel execution along the biggest dimension only. */
	#pragma omp parallel for schedule(static) if(out_ndim == pardim)
	for (int i = 0; i < SVT_len; i++) {
		SEXP subSVT = SVT == R_NilValue ? R_NilValue
						: VECTOR_ELT(SVT, i);
		void *subout = shift_dataptr(out_Rtype, out, out_inc * i);
		REC_colStats_SVT(subSVT, na_background, dim, ndim - 1,
				 summarize_op,
				 subout, out_Rtype,
				 out_incs, out_ndim - 1, pardim,
				 warn);
	}
	return;
}

/* --- .Call ENTRY POINT --- */
SEXP C_colStats_SVT(SEXP x_dim, SEXP x_dimnames, SEXP x_type,
		    SEXP x_SVT, SEXP x_na_background,
		    SEXP op, SEXP na_rm, SEXP center, SEXP dims)
{
	SEXPTYPE x_Rtype = _get_and_check_Rtype_from_Rstring(x_type,
					"C_colStats_SVT", "x_type");
	int x_has_NAbg = _get_and_check_na_background(x_na_background,
					"C_colStats_SVT", "x_na_background");

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

	int d = check_dims(dims, 1, LENGTH(x_dim));
	SEXP ans_dim = PROTECT(compute_colStats_ans_dim(x_dim, d));
	int ans_ndim = LENGTH(ans_dim);  /* = x_ndim - d */
	/* Get 1-based rank of biggest dimension. Parallel execution will
	   be along that dimension. */
	int pardim = which_max(INTEGER(ans_dim), ans_ndim) + 1;

	R_xlen_t *out_incs = NULL;
	if (ans_ndim != 0)
		out_incs = (R_xlen_t *) R_alloc(ans_ndim, sizeof(R_xlen_t));

	SEXP ans = PROTECT(alloc_ans(ans_Rtype, ans_dim, out_incs));
	propagate_colStats_dimnames(ans, x_dimnames, d);

	int warn = 0;
	REC_colStats_SVT(x_SVT, x_has_NAbg, INTEGER(x_dim), LENGTH(x_dim),
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
 * C_rowStats_SVT()                                                         *
 ****************************************************************************/


/****************************************************************************
 * A bunch of low-level helpers to support mid-level functions
 * update_out_for_rowStats() and update_out_for_rowStats_NULL(), both
 * used by C_rowStats_SVT()
 *
 * TODO: Consider moving this stuff to a different file, maybe.
 */

/* If 'narm' is True then 'out_is_not_set' is ignored. In this case '*out'
   is assumed to have been initialized with NA. */
static inline void update_out_with_int_min(int x, int narm,
		int *out, int out_is_not_set)
{
	if (narm) {
		if (x == NA_INTEGER)
			return;
		if (*out == NA_INTEGER) {
			*out = x;
			return;
		}
	} else {
		if (out_is_not_set || x == NA_INTEGER) {
			*out = x;
			return;
		}
		if (*out == NA_INTEGER)
			return;
	}
	if (x < *out)
		*out = x;
	return;
}

/* If 'narm' is True then 'out_is_not_set' is ignored. In this case '*out'
   is assumed to have been initialized with NA. */
static inline void update_out_with_int_max(int x, int narm,
		int *out, int out_is_not_set)
{
	if (narm) {
		if (x == NA_INTEGER)
			return;
		if (*out == NA_INTEGER) {
			*out = x;
			return;
		}
	} else {
		if (out_is_not_set || x == NA_INTEGER) {
			*out = x;
			return;
		}
		if (*out == NA_INTEGER)
			return;
	}
	if (x > *out)
		*out = x;
	return;
}

/* If 'narm' is True then 'out_is_not_set' is ignored. In this case '*out'
   is assumed to have been initialized with NA. */
static inline void update_out_with_double_min(double x, int narm,
		double *out, int out_is_not_set)
{
	if (narm) {
		if (ISNAN(x))  // True for *both* NA and NaN
			return;
		if (R_IsNA(*out)) {
			*out = x;
			return;
		}
	} else {
		if (out_is_not_set || R_IsNA(x)) {
			*out = x;
			return;
		}
		if (ISNAN(*out))  // True for *both* NA and NaN
			return;
		if (R_IsNaN(x)) {
			*out = x;
			return;
		}
	}
	if (x < *out)
		*out = x;
	return;
}

/* If 'narm' is True then 'out_is_not_set' is ignored. In this case '*out'
   is assumed to have been initialized with NA. */
static inline void update_out_with_double_max(double x, int narm,
		double *out, int out_is_not_set)
{
	if (narm) {
		if (ISNAN(x))  // True for *both* NA and NaN
			return;
		if (R_IsNA(*out)) {
			*out = x;
			return;
		}
	} else {
		if (out_is_not_set || R_IsNA(x)) {
			*out = x;
			return;
		}
		if (ISNAN(*out))  // True for *both* NA and NaN
			return;
		if (R_IsNaN(x)) {
			*out = x;
			return;
		}
	}
	if (x > *out)
		*out = x;
	return;
}

static inline void update_out_with_int_sum(int x, int narm, double *out)
{
	double xx;
	if (x == NA_INTEGER) {
		if (narm)
			return;
		xx = NA_REAL;
	} else {
		xx = (double) x;
	}
	*out += xx;
	return;
}

static inline void update_out_with_double_sum(double x, int narm, double *out)
{
	/* ISNAN(): True for *both* NA and NaN. See <R_ext/Arith.h> */
	if (narm && ISNAN(x))
		return;
	*out += x;
	return;
}

static inline int is_na(SEXPTYPE Rtype, const void *x, int i)
{
	switch (Rtype) {
	    case INTSXP: case LGLSXP: {
		const int *x_p = x;
		return x_p[i] == NA_INTEGER;
	    }
	    case REALSXP: {
		const double *x_p = x;
		/* ISNAN(): True for *both* NA and NaN.
		   See <R_ext/Arith.h> */
		return ISNAN(x_p[i]);
	    }
	    case CPLXSXP: {
		const Rcomplex *x_p = x;
		return RCOMPLEX_IS_NA_OR_NaN(x_p + i);
	    }
	    case RAWSXP:
		return 0;
	    case STRSXP:
		return STRING_ELT((SEXP) x, i) == NA_STRING;
	}
	error("SparseArray internal error in is_na():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return 0;  /* will never reach this */
}

static inline void add_nzval_to_out(const SparseVec *sv, int k,
				    int narm, double *out)
{
	SEXPTYPE Rtype = get_SV_Rtype(sv);
	switch (Rtype) {
	    case INTSXP: case LGLSXP: {
		update_out_with_int_sum(get_intSV_nzvals_p(sv)[k],
				narm, out);
		return;
	    }
	    case REALSXP: {
		update_out_with_double_sum(get_doubleSV_nzvals_p(sv)[k],
				narm, out);
		return;
	    }
	}
	error("SparseArray internal error in add_nzval_to_out():\n"
	      "    type \"%s\" is not supported",
	      type2char(Rtype));
}


/****************************************************************************
 * update_out_for_rowStats()
 * update_out_for_rowStats_NULL()
 *
 * The two workhorses behind rowStats_SVT().
 */

static inline void check_out_Rtype(SEXPTYPE out_Rtype, SEXPTYPE expected,
		const char *fun)
{
	if (out_Rtype == expected)
		return;
	error("SparseArray internal error in %s():\n"
	      "    out_Rtype (\"%s\") != expected out_Rtype (\"%s\")",
	      fun, type2char(out_Rtype), type2char(expected));
}

static void update_out_for_rowAnyNAs(const SparseVec *sv, int *out)
{
	if (sv->na_background)
		error("SparseArray internal error in "
		      "update_out_for_rowAnyNAs():\n"
		      "    operation not yet supported on NaArray objects");
	if (sv->nzvals == NULL)  /* lacunar leaf */
		return;
	/* regular leaf */
	int nzcount = get_SV_nzcount(sv);
	SEXPTYPE Rtype = get_SV_Rtype(sv);
	for (int k = 0; k < nzcount; k++) {
		if (is_na(Rtype, sv->nzvals, k))
			out[sv->nzoffs[k]] = 1;
	}
	return;
}

static void update_out_for_rowCountNAs(const SparseVec *sv, double *out)
{
	if (!sv->na_background && sv->nzvals == NULL)  /* lacunar leaf */
		return;
	int nzcount = get_SV_nzcount(sv);
	SEXPTYPE Rtype = get_SV_Rtype(sv);
	for (int k = 0; k < nzcount; k++) {
		if (sv->na_background) {
			if (sv->nzvals == NULL || !is_na(Rtype, sv->nzvals, k))
				out[sv->nzoffs[k]]--;
		} else {
			/* regular leaf */
			if (is_na(Rtype, sv->nzvals, k))
				out[sv->nzoffs[k]]++;
		}
	}
	return;
}

static void update_out_for_int_rowMins(const SparseVec *sv, int narm,
		int *out, R_xlen_t *nzcvg)
{
	const int *nzvals_p = get_intSV_nzvals_p(sv);
	int nzcount = get_SV_nzcount(sv);
	int x = int1;
	for (int k = 0; k < nzcount; k++) {
		if (nzvals_p != NULL)
			x = nzvals_p[k];
		int out_is_not_set = nzcvg[sv->nzoffs[k]]++ == 0;
		update_out_with_int_min(x, narm, out + sv->nzoffs[k],
					out_is_not_set);
	}
	return;
}

static void update_out_for_int_rowMaxs(const SparseVec *sv, int narm,
		int *out, R_xlen_t *nzcvg)
{
	const int *nzvals_p = get_intSV_nzvals_p(sv);
	int nzcount = get_SV_nzcount(sv);
	int x = int1;
	for (int k = 0; k < nzcount; k++) {
		if (nzvals_p != NULL)
			x = nzvals_p[k];
		int out_is_not_set = nzcvg[sv->nzoffs[k]]++ == 0;
		update_out_with_int_max(x, narm, out + sv->nzoffs[k],
					out_is_not_set);
	}
	return;
}

static void update_out_for_double_rowMins(const SparseVec *sv, int narm,
		double *out, R_xlen_t *nzcvg)
{
	const double *nzvals_p = get_doubleSV_nzvals_p(sv);
	int nzcount = get_SV_nzcount(sv);
	double x = double1;
	for (int k = 0; k < nzcount; k++) {
		if (nzvals_p != NULL)
			x = nzvals_p[k];
		int out_is_not_set = nzcvg[sv->nzoffs[k]]++ == 0;
		update_out_with_double_min(x, narm, out + sv->nzoffs[k],
					   out_is_not_set);
	}
	return;
}

static void update_out_for_double_rowMaxs(const SparseVec *sv, int narm,
		double *out, R_xlen_t *nzcvg)
{
	const double *nzvals_p = get_doubleSV_nzvals_p(sv);
	int nzcount = get_SV_nzcount(sv);
	double x = double1;
	for (int k = 0; k < nzcount; k++) {
		if (nzvals_p != NULL)
			x = nzvals_p[k];
		int out_is_not_set = nzcvg[sv->nzoffs[k]]++ == 0;
		update_out_with_double_max(x, narm, out + sv->nzoffs[k],
					   out_is_not_set);
	}
	return;
}

static void update_out_for_rowSums(const SparseVec *sv, int narm, double *out)
{
	int nzcount = get_SV_nzcount(sv);
	if (!sv->na_background || narm) {
		if (sv->nzvals == NULL) {
			/* lacunar leaf */
			for (int k = 0; k < nzcount; k++)
				out[sv->nzoffs[k]] += double1;
			return;
		}
		/* regular leaf */
		for (int k = 0; k < nzcount; k++)
			add_nzval_to_out(sv, k, narm, out + sv->nzoffs[k]);
		return;
	}
	/* The rowSums(<NaArray>, na.rm=FALSE) case.
	   Note that using 'na.rm=FALSE' on an NaArray object is atypical.
	   The implementation below is slow because we walk on _all_ the
	   values of SparseVec 'sv' i.e. on its non-NA and (implicit) NA
	   values. */
	int k = 0;
	for (int i = 0; i < sv->len; i++) {
		if (k < get_SV_nzcount(sv) && sv->nzoffs[k] == i) {
			if (sv->nzvals == NULL) {
				/* lacunar leaf */
				out[i] += double1;
			} else {
				add_nzval_to_out(sv, k, 0, out + i);
			}
			k++;
		} else {
			out[i] = doubleNA;
		}
	}
	return;
}

static void update_out_for_rowCenteredX2Sum(const SparseVec *sv,
		int narm, const double *center, double *out)
{
	if (sv->na_background)
		error("SparseArray internal error in "
		      "update_out_for_rowCenteredX2Sum():\n"
		      "    operation not yet supported on NaArray objects");
	int nzcount = get_SV_nzcount(sv);
	if (sv->nzvals == NULL) {
		/* lacunar leaf */
		for (int k = 0; k < nzcount; k++) {
			int i = sv->nzoffs[k];
			double x = double1;
			if (center != NULL)
				x -= 2 * center[i];
			out[i] += x;
		}
		return;
	}
	/* regular leaf */
	SEXPTYPE Rtype = get_SV_Rtype(sv);
	for (int k = 0; k < nzcount; k++) {
		int i = sv->nzoffs[k];
		double c = center == NULL ? 0.0 : center[i];
		double x;
		switch (Rtype) {
		    case INTSXP: case LGLSXP: {
			const int *nzvals_p = sv->nzvals;
			int v = nzvals_p[k];
			if (v == intNA) {
				if (narm) {
					out[i] -= c * c;
					continue;
				}
				x = doubleNA;
			} else {
				x = (double) v;
			}
			break;
		    }
		    case REALSXP: {
			const double *nzvals_p = sv->nzvals;
			x = nzvals_p[k];
			/* ISNAN(): True for *both* NA and NaN.
			   See <R_ext/Arith.h> */
			if (narm && ISNAN(x)) {
				out[i] -= c * c;
				continue;
			}
			break;
		    }
		    default:
			error("SparseArray internal error in "
			      "update_out_for_rowCenteredX2Sum():\n"
			      "    type \"%s\" is not supported",
			      type2char(Rtype));
		}
		out[i] += x * (x - 2 * c);
	}
	return;
}

static void update_out_for_rowStats(const SparseVec *sv,
		const SummarizeOp *summarize_op, const double *center,
		void *out, SEXPTYPE out_Rtype, R_xlen_t *nzcvg)
{
	switch (summarize_op->opcode) {
	    case ANYNA_OPCODE:
		check_out_Rtype(out_Rtype, LGLSXP, "update_out_for_rowStats");
		update_out_for_rowAnyNAs(sv, (int *) out);
		return;
	    case COUNTNAS_OPCODE:
		check_out_Rtype(out_Rtype, REALSXP, "update_out_for_rowStats");
		update_out_for_rowCountNAs(sv, (double *) out);
		return;
	    case MIN_OPCODE: case MAX_OPCODE:
		if (sv->Rtype == INTSXP || sv->Rtype == LGLSXP) {
			check_out_Rtype(out_Rtype, INTSXP,
					"update_out_for_rowStats");
			if (summarize_op->opcode == MIN_OPCODE) {
			    update_out_for_int_rowMins(sv,
						summarize_op->na_rm,
						(int *) out, nzcvg);
			} else {
			    update_out_for_int_rowMaxs(sv,
						summarize_op->na_rm,
						(int *) out, nzcvg);
			}
			return;
		}
		if (sv->Rtype == REALSXP) {
			check_out_Rtype(out_Rtype, REALSXP,
					"update_out_for_rowStats");
			if (summarize_op->opcode == MIN_OPCODE) {
			    update_out_for_double_rowMins(sv,
						summarize_op->na_rm,
						(double *) out, nzcvg);
			} else {
			    update_out_for_double_rowMaxs(sv,
						summarize_op->na_rm,
						(double *) out, nzcvg);
			}
			return;
		}
		break;
	    case SUM_OPCODE:
		check_out_Rtype(out_Rtype, REALSXP, "update_out_for_rowStats");
		update_out_for_rowSums(sv, summarize_op->na_rm,
				(double *) out);
		return;
	    case CENTERED_X2_SUM_OPCODE:
		check_out_Rtype(out_Rtype, REALSXP, "update_out_for_rowStats");
		update_out_for_rowCenteredX2Sum(sv, summarize_op->na_rm,
				center, (double *) out);
		return;
	}
	error("SparseArray internal error in update_out_for_rowStats():\n"
	      "    operation not supported");
}

static void update_out_for_rowStats_NULL(int na_background,
		const SummarizeOp *summarize_op, const double *center,
		void *out, int out_len, SEXPTYPE out_Rtype)
{
	if (!na_background)
		return;
	if (summarize_op->opcode == COUNTNAS_OPCODE || summarize_op->na_rm)
		return;
	_set_elts_to_NA(out_Rtype, out, 0, out_len);
	return;
}


/****************************************************************************
 * rowStats_SVT()
 */

/* Recursive. 'strata_counter' not used for anything at the moment. */
static void REC_rowStats_SVT(SEXP SVT, int na_background,
		const int *dim, int ndim,
		const SummarizeOp *summarize_op, const double *center,
		void *out, SEXPTYPE out_Rtype,
		const R_xlen_t *out_incs, int out_ndim,
		R_xlen_t *nzcvg, R_xlen_t *strata_counter)
{
	if (SVT == R_NilValue) {
		R_xlen_t out_len = 1;
		for (int along = 0; along < ndim && along < out_ndim; along++)
			out_len *= dim[along];
		update_out_for_rowStats_NULL(na_background,
			summarize_op, center,
			out, out_len, out_Rtype);
		if (ndim >= out_ndim) {
			R_xlen_t inc = 1;
			for (int along = out_ndim; along < ndim; along++)
				inc *= dim[along];
			*strata_counter += inc;
		}
		return;
	}

	if (ndim == 1) {
		/* 'SVT' is a leaf (i.e. a 1D SVT). */
		SparseVec sv = leaf2SV(SVT, summarize_op->in_Rtype,
				       dim[0], na_background);
		update_out_for_rowStats(&sv, summarize_op, center,
					out, out_Rtype, nzcvg);
		if (out_ndim == 1)
			(*strata_counter)++;
		return;
	}

	int SVT_len = dim[ndim - 1];
	R_xlen_t out_inc = 0;
	if (ndim <= out_ndim)
		out_inc = out_incs[ndim - 1];
	for (int i = 0; i < SVT_len; i++) {
		SEXP subSVT = VECTOR_ELT(SVT, i);
		const double *subcenter = NULL;
		if (center != NULL)
			subcenter = center + out_inc * i;
		void *subout = shift_dataptr(out_Rtype, out, out_inc * i);
		R_xlen_t *subnzcvg = NULL;
		if (nzcvg != NULL)
			subnzcvg = nzcvg + out_inc * i;
		REC_rowStats_SVT(subSVT, na_background, dim, ndim - 1,
				 summarize_op, subcenter,
				 subout, out_Rtype, out_incs, out_ndim,
				 subnzcvg, strata_counter);
	}
	if (ndim == out_ndim)
		(*strata_counter)++;
	return;
}

static void rowStats_SVT(SEXP SVT, int na_background,
		const int *dim, int ndim,
		const SummarizeOp *summarize_op, const double *center,
		void *out, SEXPTYPE out_Rtype,
		const R_xlen_t *out_incs, int out_ndim,
		R_xlen_t *nzcvg, R_xlen_t nstrata)
{
	R_xlen_t strata_counter = 0;
	REC_rowStats_SVT(SVT, na_background, dim, ndim,
			 summarize_op, center,
			 out, out_Rtype, out_incs, out_ndim,
			 nzcvg, &strata_counter);
	//printf("strata_counter = %ld\n", strata_counter);
	if (strata_counter != nstrata)
		error("SparseArray internal error in "
		      "rowStats_SVT():\n"
		      "    strata_counter != nstrata");
	return;
}


/****************************************************************************
 * The SVT_row*() functions
 */

static void SVT_rowCountNAs(SEXP SVT, SEXPTYPE Rtype, int na_background,
		const int *dim, int ndim,
		double *out, R_xlen_t out_len,
		const R_xlen_t *out_incs, int out_ndim,
		R_xlen_t nstrata)
{
	/* Initialization. */
	double v = na_background ? (double) nstrata : 0.0;
	for (R_xlen_t i = 0; i < out_len; i++)
		out[i] = v;

	/* Recursive tree traversal. */
	if (nstrata != 0) {
		SummarizeOp summarize_op =
			_make_SummarizeOp(COUNTNAS_OPCODE,
					  Rtype, 0, NA_REAL);
		rowStats_SVT(SVT, na_background, dim, ndim,
			     &summarize_op, NULL,
			     out, REALSXP, out_incs, out_ndim,
			     NULL, nstrata);
	}
	return;
}

static void SVT_rowAnyNAs(SEXP SVT, SEXPTYPE Rtype, int na_background,
		const int *dim, int ndim,
		int *out, R_xlen_t out_len,
		const R_xlen_t *out_incs, int out_ndim,
		R_xlen_t nstrata)
{
	if (!na_background) {
		/* Initialization. */
		_set_elts_to_zero(LGLSXP, out, 0, out_len);

		/* Recursive tree traversal. */
		if (nstrata != 0) {
			SummarizeOp summarize_op =
				_make_SummarizeOp(ANYNA_OPCODE,
						  Rtype, 0, NA_REAL);
			rowStats_SVT(SVT, 0, dim, ndim,
				     &summarize_op, NULL,
				     out, LGLSXP, out_incs, out_ndim,
				     NULL, nstrata);
		}
		return;
	}

	double *count_NAs = (double *) R_alloc(out_len, sizeof(double));
	SVT_rowCountNAs(SVT, Rtype, 1,
			dim, ndim,
			count_NAs, out_len,
			out_incs, out_ndim,
			nstrata);
	for (R_xlen_t i = 0; i < out_len; i++)
		out[i] = count_NAs[i] != 0.0;
	return;
}

static void postprocess_int_rowMinsMaxs(int na_background,
		int opcode, int narm,
		int *out, R_xlen_t out_len,
		const R_xlen_t *nzcvg, R_xlen_t nstrata)
{
	int warn = 0;
	for (R_xlen_t i = 0; i < out_len; i++) {
		R_xlen_t nzcvg_i = nzcvg[i];
		if (nzcvg_i < nstrata) {
			int background = na_background ? intNA : int0;
			if (opcode == MIN_OPCODE)
				update_out_with_int_min(background, narm,
							out + i, nzcvg_i == 0);
			else
				update_out_with_int_max(background, narm,
							out + i, nzcvg_i == 0);
		}
		if (narm && out[i] == NA_INTEGER)
			warn = 1;
	}
	if (warn)
		warning("NAs introduced by coercion of "
			"infinite values to integers");
	return;
}

static void postprocess_double_rowMinsMaxs(int na_background,
		int opcode, int narm,
		double *out, R_xlen_t out_len,
		const R_xlen_t *nzcvg, R_xlen_t nstrata)
{
	double inf = opcode == MIN_OPCODE ? R_PosInf : R_NegInf;
	for (R_xlen_t i = 0; i < out_len; i++) {
		R_xlen_t nzcvg_i = nzcvg[i];
		if (nzcvg_i < nstrata) {
			double background = na_background ? doubleNA : double0;
			if (opcode == MIN_OPCODE)
				update_out_with_double_min(background, narm,
							out + i, nzcvg_i == 0);
			else
				update_out_with_double_max(background, narm,
							out + i, nzcvg_i == 0);
		}
		if (narm && R_IsNA(out[i]))
			out[i] = inf;
	}
	return;
}

static void SVT_rowMinsMaxs(SEXP SVT, SEXPTYPE Rtype, int na_background,
		const int *dim, int ndim,
		int opcode, int narm,
		void *out, R_xlen_t out_len, SEXPTYPE out_Rtype,
		const R_xlen_t *out_incs, int out_ndim,
		R_xlen_t nstrata)
{
	/* Initialization. */

	if (nstrata == 0) {
		if (out_Rtype == REALSXP) {
			double inf = opcode == MIN_OPCODE ? R_PosInf : R_NegInf;
			double *out_p = (double *) out;
			for (R_xlen_t i = 0; i < out_len; i++)
				out_p[i] = inf;
		} else {
			_set_elts_to_NA(out_Rtype, out, 0, out_len);
			warning("NAs introduced by coercion of "
				"infinite values to integers");
		}
		return;
	}
	if (narm)
		_set_elts_to_NA(out_Rtype, out, 0, out_len);

	/* Recursive tree traversal. */

	SummarizeOp summarize_op =
		_make_SummarizeOp(opcode, Rtype, narm, NA_REAL);
	/* 'nzcvg' is an array "parallel" to 'out' (i.e. same dimensions).
	   We will use it to compute the coverage of 'out' as we
	   traverse 'SVT'. More precisely, each 'nzcvg[i]' will be used
	   to keep track of the number of times that 'out[i]' is visited
	   as we walk on all the leaves in 'SVT'. Note that 'nzcvg[i]' is
	   guaranteed to be >= 0 and <= nstrata. */
	R_xlen_t *nzcvg = (R_xlen_t *) R_alloc(out_len, sizeof(R_xlen_t));
	memset(nzcvg, 0, sizeof(R_xlen_t) * out_len);
	rowStats_SVT(SVT, na_background, dim, ndim,
		     &summarize_op, NULL,
		     out, out_Rtype, out_incs, out_ndim,
		     nzcvg, nstrata);

	/* Postprocessing. */

	if (out_Rtype == INTSXP) {
		postprocess_int_rowMinsMaxs(na_background,
				opcode, narm,
				(int *) out, out_len,
				nzcvg, nstrata);
	} else {
		postprocess_double_rowMinsMaxs(na_background,
				opcode, narm,
				(double *) out, out_len,
				nzcvg, nstrata);
	}
	return;
}

static void SVT_rowSums(SEXP SVT, SEXPTYPE Rtype, int na_background,
		const int *dim, int ndim,
		int narm,
		double *out, R_xlen_t out_len,
		const R_xlen_t *out_incs, int out_ndim,
		R_xlen_t nstrata)
{
	/* Initialization. */
	_set_elts_to_zero(REALSXP, out, 0, out_len);

	/* Recursive tree traversal. */
	if (nstrata != 0) {
		SummarizeOp summarize_op =
			_make_SummarizeOp(SUM_OPCODE,
					  Rtype, narm, NA_REAL);
		rowStats_SVT(SVT, na_background, dim, ndim,
			     &summarize_op, NULL,
			     out, REALSXP, out_incs, out_ndim,
			     NULL, nstrata);
	}
	return;
}

static void SVT_rowCenteredX2Sum(SEXP SVT, SEXPTYPE Rtype, int na_background,
		const int *dim, int ndim,
		int narm, const double *center,
		double *out, R_xlen_t out_len,
		const R_xlen_t *out_incs, int out_ndim,
		R_xlen_t nstrata)
{
	/* Initialization. */
	if (center == NULL) {
		_set_elts_to_zero(REALSXP, out, 0, out_len);
	} else {
		for (R_xlen_t i = 0; i < out_len; i++) {
			double c = center[i];
			out[i] = c * c * nstrata;
		}
	}

	/* Recursive tree traversal. */
	if (nstrata != 0) {
		SummarizeOp summarize_op =
			_make_SummarizeOp(CENTERED_X2_SUM_OPCODE,
					  Rtype, narm, NA_REAL);
		rowStats_SVT(SVT, na_background, dim, ndim,
			     &summarize_op, center,
			     out, REALSXP, out_incs, out_ndim,
			     NULL, nstrata);
	}
	return;
}


/****************************************************************************
 * C_rowStats_SVT()
 */

static const double *check_rowStats_center(SEXP center,
		SEXP x_dim, int ans_ndim)
{
	if (center == R_NilValue)
		return NULL;
	if (!IS_NUMERIC(center))
		error("SparseArray internal error in check_rowStats_center():\n"
		      "    'center' must be NULL or a numeric array");
	/* We only check the length of 'center', not its actual dimensions. */
	R_xlen_t ans_len = 1;
	for (int along = 0; along < ans_ndim; along++)
		ans_len *= INTEGER(x_dim)[along];
	if (LENGTH(center) != ans_len)
		error("SparseArray internal error in check_rowStats_center():\n"
		      "    unexpected 'center' length");
	return REAL(center);
}

/* Compute 'nstrata', the number of "strata" in an array-like object 'x'.
   The "strata" are the generalized **columns** in 'x'. Their number and
   geometry are determined by the dimensions of 'x' and the 'dims' argument:
     o dimensions of each strata = head(dim(x), n=dims)
     o number of strata = prod(tail(dim(x), n=-dims))
   Note that:
     o The dimensions of a strata are the dimensions of 'ans', where 'ans'
       is the output of a row*() function.
     o The number of strata is the number of array elements in each
       generalized row.
     o length(ans) * nstrata = length(x)
     o 'x' can be seen as a pile of strata stacked up on top of each others.
       An row*() function summarizes vertically i.e. **across** strata, and
       not **within** each strata. */
static R_xlen_t compute_nstrata(const int *x_dim, int x_ndim, int dims)
{
	R_xlen_t nstrata = 1;
	for (int along = dims; along < x_ndim; along++)
		nstrata *= x_dim[along];
	//printf("nstrata = %ld\n", nstrata);
	return nstrata;
}

/* --- .Call ENTRY POINT --- */
SEXP C_rowStats_SVT(SEXP x_dim, SEXP x_dimnames, SEXP x_type,
		    SEXP x_SVT, SEXP x_na_background,
		    SEXP op, SEXP na_rm, SEXP center, SEXP dims)
{
	SEXPTYPE x_Rtype = _get_and_check_Rtype_from_Rstring(x_type,
					"C_rowStats_SVT", "x_type");
	int x_has_NAbg = _get_and_check_na_background(x_na_background,
					"C_colStats_SVT", "x_na_background");

	int opcode = _get_summarize_opcode(op, x_Rtype);

	if (!(IS_LOGICAL(na_rm) && LENGTH(na_rm) == 1))
		error("'na.rm' must be TRUE or FALSE");
	int narm = LOGICAL(na_rm)[0];

	int ans_ndim = check_dims(dims, 1, LENGTH(x_dim) - 1);
	const double *center_p = check_rowStats_center(center, x_dim, ans_ndim);

	/* 4th argument (center) will be ignored. */
	SummarizeOp summarize_op = _make_SummarizeOp(opcode, x_Rtype, narm,
						     NA_REAL);
	SEXPTYPE ans_Rtype = compute_ans_Rtype(&summarize_op);

	SEXP ans_dim = PROTECT(compute_rowStats_ans_dim(x_dim, ans_ndim));

	/* 'ans_ndim' is guaranteed to be >= 1. */
	R_xlen_t *out_incs = (R_xlen_t *) R_alloc(ans_ndim, sizeof(R_xlen_t));

	SEXP ans = PROTECT(alloc_ans(ans_Rtype, ans_dim, out_incs));
	propagate_rowStats_dimnames(ans, x_dimnames, ans_ndim);

	if (LENGTH(ans) == 0) {
		UNPROTECT(2);
		return ans;
	}

	R_xlen_t nstrata = compute_nstrata(INTEGER(x_dim),
					   LENGTH(x_dim), ans_ndim);
	switch (opcode) {
	    case COUNTNAS_OPCODE:
		check_out_Rtype(ans_Rtype, REALSXP, "C_rowStats_SVT");
		SVT_rowCountNAs(x_SVT, x_Rtype, x_has_NAbg,
			INTEGER(x_dim), LENGTH(x_dim),
			REAL(ans), LENGTH(ans),
			out_incs, ans_ndim, nstrata);
		break;
	    case ANYNA_OPCODE:
		check_out_Rtype(ans_Rtype, LGLSXP, "C_rowStats_SVT");
		SVT_rowAnyNAs(x_SVT, x_Rtype, x_has_NAbg,
			INTEGER(x_dim), LENGTH(x_dim),
			LOGICAL(ans), LENGTH(ans),
			out_incs, ans_ndim, nstrata);
		break;
	    case MIN_OPCODE: case MAX_OPCODE:
		SVT_rowMinsMaxs(x_SVT, x_Rtype, x_has_NAbg,
			INTEGER(x_dim), LENGTH(x_dim),
			opcode, narm,
			DATAPTR(ans), LENGTH(ans), ans_Rtype,
			out_incs, ans_ndim, nstrata);
		break;
	    case SUM_OPCODE:
		check_out_Rtype(ans_Rtype, REALSXP, "C_rowStats_SVT");
		SVT_rowSums(x_SVT, x_Rtype, x_has_NAbg,
			INTEGER(x_dim), LENGTH(x_dim),
			narm,
			REAL(ans), LENGTH(ans),
			out_incs, ans_ndim, nstrata);
		break;
	    case CENTERED_X2_SUM_OPCODE:
		check_out_Rtype(ans_Rtype, REALSXP, "C_rowStats_SVT");
		SVT_rowCenteredX2Sum(x_SVT, x_Rtype, x_has_NAbg,
			INTEGER(x_dim), LENGTH(x_dim),
			narm, center_p,
			REAL(ans), LENGTH(ans),
			out_incs, ans_ndim, nstrata);

		break;
	    default:
		error("SparseArray internal error in C_rowStats_SVT():\n"
		      "    operation not supported");
	}

	UNPROTECT(2);
	return ans;
}

