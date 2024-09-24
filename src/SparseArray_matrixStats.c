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

static void init_rowStats_ans(SEXP ans, const SummarizeOp *summarize_op,
		const double *center_p, R_xlen_t nstrata)
{
	if (summarize_op->opcode == MIN_OPCODE ||
	    summarize_op->opcode == MAX_OPCODE)
	{
		if (nstrata != 0) {
			_set_Rvector_elts_to_zero(ans);
			return;
		}
		if (TYPEOF(ans) == REALSXP) {
			double init = summarize_op->opcode == MIN_OPCODE ?
						R_PosInf : R_NegInf;
			double *ans_p = REAL(ans);
			for (int i = 0; i < LENGTH(ans); i++)
				ans_p[i] = init;
			return;
		}
		_set_Rvector_elts_to_NA(ans);
		warning("NAs introduced by coercion of "
			"infinite values to integers");
		return;
	}

	if (summarize_op->opcode != CENTERED_X2_SUM_OPCODE ||
	    center_p == NULL)
	{
		_set_Rvector_elts_to_zero(ans);
		return;
	}
	double *ans_p = REAL(ans);
	for (int i = 0; i < LENGTH(ans); i++) {
		double c = center_p[i];
		ans_p[i] = nstrata * c * c;
	}
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
 * A bunch of low-level helpers to support the various row*_SV() functions
 * used by C_rowStats_SVT()
 *
 * TODO: Consider moving this stuff to a different file, maybe.
 */

static inline void update_out_with_int_min(int x, int narm, int *out)
{
	if (narm) {
		/* '*out' should never be NA when 'narm' is TRUE. */
		if (x == NA_INTEGER)
			return;
	} else {
		if (*out == NA_INTEGER)
			return;
		if (x == NA_INTEGER) {
			*out = NA_INTEGER;
			return;
		}
	}
	if (x < *out)
		*out = x;
	return;
}

static inline void update_out_with_int_max(int x, int narm, int *out)
{
	if (narm) {
		/* '*out' should never be NA when 'narm' is TRUE. */
		if (x == NA_INTEGER)
			return;
	} else {
		if (*out == NA_INTEGER)
			return;
		if (x == NA_INTEGER) {
			*out = NA_INTEGER;
			return;
		}
	}
	if (x > *out)
		*out = x;
	return;
}

static inline void update_out_with_double_min(double x, int narm, double *out)
{
	if (narm) {
		/* '*out' should never be NA or NaN when 'narm' is TRUE. */
		if (ISNAN(x))  // True for *both* NA and NaN
			return;
	} else {
		if (R_IsNA(*out))
			return;
		if (R_IsNA(x)) {
			*out = NA_REAL;
			return;
		}
		if (R_IsNaN(*out))
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

static inline void update_out_with_double_max(double x, int narm, double *out)
{
	if (narm) {
		/* '*out' should never be NA or NaN when 'narm' is TRUE. */
		if (ISNAN(x))  // True for *both* NA and NaN
			return;
	} else {
		if (R_IsNA(*out))
			return;
		if (R_IsNA(x)) {
			*out = NA_REAL;
			return;
		}
		if (R_IsNaN(*out))
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
 * C_rowStats_SVT()
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
	if (sv->na_background)
		error("SparseArray internal error in "
		      "update_out_for_rowCountNAs():\n"
		      "    operation not yet supported on NaArray objects");
	if (sv->nzvals == NULL)  /* lacunar leaf */
		return;
	/* regular leaf */
	int nzcount = get_SV_nzcount(sv);
	SEXPTYPE Rtype = get_SV_Rtype(sv);
	for (int k = 0; k < nzcount; k++)
		out[sv->nzoffs[k]] += is_na(Rtype, sv->nzvals, k);
	return;
}

static void update_out_for_int_rowMins(const SparseVec *sv, int narm, int *out)
{
	int nzcount = get_SV_nzcount(sv);
	if (!sv->na_background || narm) {
		if (sv->nzvals == NULL) {
			/* lacunar leaf */
			return;
		}
		/* regular leaf */
		for (int k = 0; k < nzcount; k++)
			update_out_with_int_min(get_intSV_nzvals_p(sv)[k],
						narm, out + sv->nzoffs[k]);
		return;
	}
	/* The rowMins(<NaArray>, na.rm=FALSE) case.
	   Note that using 'na.rm=FALSE' on an NaArray object is atypical.
	   The implementation below is slow because we walk on _all_ the
	   values of SparseVec 'sv' i.e. on its non-NA and (implicit) NA
	   values. */
	int k = 0;
	for (int i = 0; i < sv->len; i++) {
		if (k < get_SV_nzcount(sv) && sv->nzoffs[k] == i) {
			update_out_with_int_min(get_intSV_nzval(sv, k),
						0, out + i);
			k++;
		} else {
			out[i] = intNA;
		}
	}
	return;
}

static void update_out_for_int_rowMaxs(const SparseVec *sv, int narm, int *out)
{
	int nzcount = get_SV_nzcount(sv);
	if (!sv->na_background || narm) {
		if (sv->nzvals == NULL) {
			/* lacunar leaf */
			for (int k = 0; k < nzcount; k++)
				update_out_with_int_max(int1,
						narm, out + sv->nzoffs[k]);
			return;
		}
		/* regular leaf */
		for (int k = 0; k < nzcount; k++)
			update_out_with_int_max(get_intSV_nzvals_p(sv)[k],
						narm, out + sv->nzoffs[k]);
		return;
	}
	/* The rowMaxs(<NaArray>, na.rm=FALSE) case.
	   Note that using 'na.rm=FALSE' on an NaArray object is atypical.
	   The implementation below is slow because we walk on _all_ the
	   values of SparseVec 'sv' i.e. on its non-NA and (implicit) NA
	   values. */
	int k = 0;
	for (int i = 0; i < sv->len; i++) {
		if (k < get_SV_nzcount(sv) && sv->nzoffs[k] == i) {
			update_out_with_int_max(get_intSV_nzval(sv, k),
						0, out + i);
			k++;
		} else {
			out[i] = intNA;
		}
	}
	return;
}


/* Same as 'update_out_for_int_rowMins()' above but operates on an implicit
   vector of zeros instead of 'sv'. Said otherwise, it computes:
     update_out_for_int_rowMins(<vector-of-zeros>, 0, out) */
static void update_out_for_int_rowMins0(int *out, R_xlen_t out_len)
{
	for (int i = 0; i < out_len; i++)
		update_out_with_int_min(0, 0, out + i);
	return;
}

/* Same as 'update_out_for_int_rowMaxs()' above but operates on an implicit
   vector of zeros instead of 'sv'. Said otherwise, it computes:
     update_out_for_int_rowMaxs(<vector-of-zeros>, 0, out) */
static void update_out_for_int_rowMaxs0(int *out, R_xlen_t out_len)
{
	for (int i = 0; i < out_len; i++)
		update_out_with_int_max(0, 0, out + i);
	return;
}

static void update_out_for_double_rowMins(const SparseVec *sv, int narm,
		double *out)
{
	int nzcount = get_SV_nzcount(sv);
	if (!sv->na_background || narm) {
		if (sv->nzvals == NULL) {
			/* lacunar leaf */
			return;
		}
		/* regular leaf */
		for (int k = 0; k < nzcount; k++)
			update_out_with_double_min(get_doubleSV_nzvals_p(sv)[k],
						   narm, out + sv->nzoffs[k]);
		return;
	}
	/* The rowMins(<NaArray>, na.rm=FALSE) case.
	   Note that using 'na.rm=FALSE' on an NaArray object is atypical.
	   The implementation below is slow because we walk on _all_ the
	   values of SparseVec 'sv' i.e. on its non-NA and (implicit) NA
	   values. */
	int k = 0;
	for (int i = 0; i < sv->len; i++) {
		if (k < get_SV_nzcount(sv) && sv->nzoffs[k] == i) {
			update_out_with_double_min(get_doubleSV_nzval(sv, k),
						   0, out + i);
			k++;
		} else {
			out[i] = doubleNA;
		}
	}
	return;
}

static void update_out_for_double_rowMaxs(const SparseVec *sv, int narm,
		double *out)
{
	int nzcount = get_SV_nzcount(sv);
	if (!sv->na_background || narm) {
		if (sv->nzvals == NULL) {
			/* lacunar leaf */
			for (int k = 0; k < nzcount; k++)
				update_out_with_double_max(double1,
						   narm, out + sv->nzoffs[k]);
			return;
		}
		/* regular leaf */
		for (int k = 0; k < nzcount; k++)
			update_out_with_double_max(get_doubleSV_nzvals_p(sv)[k],
						   narm, out + sv->nzoffs[k]);
		return;
	}
	/* The rowMaxs(<NaArray>, na.rm=FALSE) case.
	   Note that using 'na.rm=FALSE' on an NaArray object is atypical.
	   The implementation below is slow because we walk on _all_ the
	   values of SparseVec 'sv' i.e. on its non-NA and (implicit) NA
	   values. */
	int k = 0;
	for (int i = 0; i < sv->len; i++) {
		if (k < get_SV_nzcount(sv) && sv->nzoffs[k] == i) {
			update_out_with_double_max(get_doubleSV_nzval(sv, k),
						   0, out + i);
			k++;
		} else {
			out[i] = doubleNA;
		}
	}
	return;
}

static void update_out_for_double_rowMins0(double *out, R_xlen_t out_len)
{
	for (int i = 0; i < out_len; i++)
		update_out_with_double_min(0.0, 0, out + i);
	return;
}

static void update_out_for_double_rowMaxs0(double *out, R_xlen_t out_len)
{
	for (int i = 0; i < out_len; i++)
		update_out_with_double_max(0.0, 0, out + i);
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
		void *out, SEXPTYPE out_Rtype, R_xlen_t strata_counter)
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
					summarize_op->na_rm, (int *) out);
			} else {
				update_out_for_int_rowMaxs(sv,
					summarize_op->na_rm, (int *) out);
			}
			return;
		}
		if (sv->Rtype == REALSXP) {
			check_out_Rtype(out_Rtype, REALSXP,
					"update_out_for_rowStats");
			if (summarize_op->opcode == MIN_OPCODE) {
				update_out_for_double_rowMins(sv,
					summarize_op->na_rm, (double *) out);
			} else {
				update_out_for_double_rowMaxs(sv,
					summarize_op->na_rm, (double *) out);
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

static void update_out_for_rowStats0(
		const SummarizeOp *summarize_op, const double *center,
		void *out, R_xlen_t out_len, SEXPTYPE out_Rtype,
		R_xlen_t strata_counter)
{
	SEXPTYPE in_Rtype = summarize_op->in_Rtype;
	switch (summarize_op->opcode) {
	    case MIN_OPCODE: case MAX_OPCODE:
		if (in_Rtype == INTSXP || in_Rtype == LGLSXP) {
			check_out_Rtype(out_Rtype, INTSXP,
					"update_out_for_rowStats0");
			if (summarize_op->opcode == MIN_OPCODE) {
				update_out_for_int_rowMins0((int *) out,
							    out_len);
			} else {
				update_out_for_int_rowMaxs0((int *) out,
							    out_len);
			}
			return;
		}
		if (in_Rtype == REALSXP) {
			check_out_Rtype(out_Rtype, REALSXP,
					"update_out_for_rowStats0");
			if (summarize_op->opcode == MIN_OPCODE) {
				update_out_for_double_rowMins0((double *) out,
							       out_len);
			} else {
				update_out_for_double_rowMaxs0((double *) out,
							       out_len);
			}
			return;
		}
		break;
	    default:
		return;
	}
	error("SparseArray internal error in update_out_for_rowStats0():\n"
	      "    unsupported 'in_Rtype': \"%s\"", type2char(in_Rtype));
}

static void update_out_for_rowStats_NULL(int na_background,
		const SummarizeOp *summarize_op, const double *center,
		void *out, int out_len, SEXPTYPE out_Rtype,
		R_xlen_t strata_counter)
{
	if (na_background) {
		if (!summarize_op->na_rm)
			_set_elts_to_NA(out_Rtype, out, 0, out_len);
		return;
	}
	update_out_for_rowStats0(summarize_op, center,
				 out, out_len, out_Rtype,
				 strata_counter);
	return;
}

/* Recursive. */
static void REC_rowStats_SVT(SEXP SVT, int na_background,
		const int *dim, int ndim,
		const SummarizeOp *summarize_op, const double *center,
		void *out, SEXPTYPE out_Rtype,
		const R_xlen_t *out_incs, int out_ndim,
		R_xlen_t *strata_counter)
{
	if (SVT == R_NilValue) {
		R_xlen_t out_len = 1;
		for (int along = 0; along < ndim && along < out_ndim; along++)
			out_len *= dim[along];
		update_out_for_rowStats_NULL(na_background,
			summarize_op, center,
			out, out_len, out_Rtype,
			*strata_counter);
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
					out, out_Rtype, *strata_counter);
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
		const double *subcenter = center + out_inc * i;
		void *subout = shift_dataptr(out_Rtype, out, out_inc * i);
		REC_rowStats_SVT(subSVT, na_background, dim, ndim - 1,
				 summarize_op, subcenter,
				 subout, out_Rtype, out_incs, out_ndim,
				 strata_counter);
	}
	if (ndim == out_ndim)
		(*strata_counter)++;
	return;
}

static R_xlen_t rowStats_SVT(SEXP SVT, int na_background,
		const int *dim, int ndim,
		const SummarizeOp *summarize_op, const double *center,
		void *out, SEXPTYPE out_Rtype,
		const R_xlen_t *out_incs, int out_ndim)
{
	R_xlen_t strata_counter = 0;
	REC_rowStats_SVT(SVT, na_background, dim, ndim,
			 summarize_op, center,
			 out, out_Rtype, out_incs, out_ndim,
			 &strata_counter);
	return strata_counter;
}

static void postprocess_rowStats_SVT_out(SEXP SVT, int na_background,
		const int *dim, int ndim,
		const SummarizeOp *summarize_op, const double *center,
		void *out, R_xlen_t out_len, SEXPTYPE out_Rtype,
		const R_xlen_t *out_incs, int out_ndim, R_xlen_t nstrata)
{
	if ((summarize_op->opcode != MIN_OPCODE &&
	     summarize_op->opcode != MAX_OPCODE) || !summarize_op->na_rm)
	{
		return;
	}

	double *count_NAs = (double *) R_alloc(out_len, sizeof(double));
	_set_elts_to_zero(REALSXP, count_NAs, 0, out_len);

	SummarizeOp summarize_op2 =
		_make_SummarizeOp(COUNTNAS_OPCODE, summarize_op->in_Rtype,
				  0, NA_REAL);
	rowStats_SVT(SVT, na_background, dim, ndim,
		     &summarize_op2, NULL,
		     count_NAs, REALSXP, out_incs, out_ndim);
	if (out_Rtype == INTSXP || out_Rtype == LGLSXP) {
		int warn = 0;
		int *out_p = (int *) out;
		for (int i = 0; i < out_len; i++) {
			if (count_NAs[i] == nstrata) {
				out_p[i] = NA_INTEGER;
				warn = 1;
			}
		}
		if (warn)
			warning("NAs introduced by coercion of "
				"infinite values to integers");
	} else {
		double v = summarize_op->opcode == MIN_OPCODE ?
						   R_PosInf : R_NegInf;
		double *out_p = (double *) out;
		for (int i = 0; i < out_len; i++) {
			if (count_NAs[i] == nstrata)
				out_p[i] = v;
		}
	}
	return;
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

	R_xlen_t nstrata = 1;
	for (int along = ans_ndim; along < LENGTH(x_dim); along++)
		nstrata *= INTEGER(x_dim)[along];
	//printf("nstrata = %ld\n", nstrata);

	init_rowStats_ans(ans, &summarize_op, center_p, nstrata);

	if (nstrata != 0) {
		R_xlen_t strata_counter =
			rowStats_SVT(x_SVT, x_has_NAbg,
				     INTEGER(x_dim), LENGTH(x_dim),
				     &summarize_op, center_p,
				     DATAPTR(ans), ans_Rtype,
				     out_incs, ans_ndim);
		//printf("strata_counter = %ld\n", strata_counter);
		if (strata_counter != nstrata)
			error("SparseArray internal error in "
			      "C_rowStats_SVT():\n"
			      "    strata_counter != nstrata");
		postprocess_rowStats_SVT_out(x_SVT, x_has_NAbg,
				     INTEGER(x_dim), LENGTH(x_dim),
				     &summarize_op, center_p,
				     DATAPTR(ans), LENGTH(ans), ans_Rtype,
				     out_incs, ans_ndim, nstrata);
	}

	UNPROTECT(2);
	return ans;
}

