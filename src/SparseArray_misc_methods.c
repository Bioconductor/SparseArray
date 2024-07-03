/*****************************************************************************
 *              Miscellaneous operations on SparseArray objects              *
 *****************************************************************************/
#include "SparseArray_misc_methods.h"

#include "leaf_utils.h"

#include <string.h>  /* for memcpy(), strcmp() */


static int collect_int_na_nzoffs(const int *nzvals,
		const int *nzoffs, int nzcount, int *out_nzoffs)
{
	int out_nzcount = 0;
	for (int k = 0; k < nzcount; k++) {
		if (nzvals[k] == NA_INTEGER) {
			out_nzoffs[out_nzcount] = nzoffs[k];
			out_nzcount++;
		}
	}
	return out_nzcount;
}

static int collect_double_na_nzoffs(const double *nzvals,
		const int *nzoffs, int nzcount, int *out_nzoffs)
{
	int out_nzcount = 0;
	for (int k = 0; k < nzcount; k++) {
		double x = nzvals[k];
		if (ISNAN(x)) {  /* True for *both* NA and NaN. */
			out_nzoffs[out_nzcount] = nzoffs[k];
			out_nzcount++;
		}
	}
	return out_nzcount;
}

static int collect_double_nan_nzoffs(const double *nzvals,
		const int *nzoffs, int nzcount, int *out_nzoffs)
{
	int out_nzcount = 0;
	for (int k = 0; k < nzcount; k++) {
		double x = nzvals[k];
		if (R_IsNaN(x)) {  /* True for special NaN, *not* for NA. */
			out_nzoffs[out_nzcount] = nzoffs[k];
			out_nzcount++;
		}
	}
	return out_nzcount;
}

#define	IS_INFINITE(x) ((x) == R_PosInf || (x) == R_NegInf)

static int collect_double_infinite_nzoffs(const double *nzvals,
		const int *nzoffs, int nzcount, int *out_nzoffs)
{
	int out_nzcount = 0;
	for (int k = 0; k < nzcount; k++) {
		double x = nzvals[k];
		if (IS_INFINITE(x)) {
			out_nzoffs[out_nzcount] = nzoffs[k];
			out_nzcount++;
		}
	}
	return out_nzcount;
}

static int collect_Rcomplex_na_nzoffs(const Rcomplex *nzvals,
		const int *nzoffs, int nzcount, int *out_nzoffs)
{
	int out_nzcount = 0;
	for (int k = 0; k < nzcount; k++) {
		const Rcomplex *z = nzvals + k;
		if (ISNAN(z->r) || ISNAN(z->i)) {
			out_nzoffs[out_nzcount] = nzoffs[k];
			out_nzcount++;
		}
	}
	return out_nzcount;
}

static int collect_Rcomplex_infinite_nzoffs(const Rcomplex *nzvals,
		const int *nzoffs, int nzcount, int *out_nzoffs)
{
	int out_nzcount = 0;
	for (int k = 0; k < nzcount; k++) {
		const Rcomplex *z = nzvals + k;
		if (IS_INFINITE(z->r) || IS_INFINITE(z->i)) {
			out_nzoffs[out_nzcount] = nzoffs[k];
			out_nzcount++;
		}
	}
	return out_nzcount;
}

static int collect_Rcomplex_nan_nzoffs(const Rcomplex *nzvals,
		const int *nzoffs, int nzcount, int *out_nzoffs)
{
	int out_nzcount = 0;
	for (int k = 0; k < nzcount; k++) {
		const Rcomplex *z = nzvals + k;
		if (R_IsNaN(z->r) || R_IsNaN(z->i)) {
			out_nzoffs[out_nzcount] = nzoffs[k];
			out_nzcount++;
		}
	}
	return out_nzcount;
}

static int collect_character_na_nzoffs(SEXP nzvals,
		const int *nzoffs, int nzcount, int *out_nzoffs)
{
	int out_nzcount = 0;
	for (int k = 0; k < nzcount; k++) {
		if (STRING_ELT(nzvals, k) == NA_STRING) {
			out_nzoffs[out_nzcount] = nzoffs[k];
			out_nzcount++;
		}
	}
	return out_nzcount;
}


/*****************************************************************************
 * get_CollectnzoffsFUN()
 */

typedef int (*CollectnzoffsFUN)(SEXP nzvals,
		const int *nzoffs, int nzcount, int *out_nzoffs);

static int collect_na_nzoffs(SEXP nzvals,
		const int *nzoffs, int nzcount, int *out_nzoffs)
{
	SEXPTYPE Rtype = TYPEOF(nzvals);
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		return collect_int_na_nzoffs(INTEGER(nzvals),
					nzoffs, nzcount, out_nzoffs);
	    case REALSXP:
		return collect_double_na_nzoffs(REAL(nzvals),
					nzoffs, nzcount, out_nzoffs);
	    case CPLXSXP:
		return collect_Rcomplex_na_nzoffs(COMPLEX(nzvals),
					nzoffs, nzcount, out_nzoffs);
	    case STRSXP:
		return collect_character_na_nzoffs(nzvals,
					nzoffs, nzcount, out_nzoffs);
	}
	error("SparseArray internal error in collect_na_nzoffs():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return 0;  /* will never reach this */
}

static int collect_nan_nzoffs(SEXP nzvals,
		const int *nzoffs, int nzcount, int *out_nzoffs)
{
	SEXPTYPE Rtype = TYPEOF(nzvals);
	switch (Rtype) {
	    case REALSXP:
		return collect_double_nan_nzoffs(REAL(nzvals),
					nzoffs, nzcount, out_nzoffs);
	    case CPLXSXP:
		return collect_Rcomplex_nan_nzoffs(COMPLEX(nzvals),
					nzoffs, nzcount, out_nzoffs);
	}
	error("SparseArray internal error in collect_nan_nzoffs():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return 0;  /* will never reach this */
}

static int collect_infinite_nzoffs(SEXP nzvals,
		const int *nzoffs, int nzcount, int *out_nzoffs)
{
	SEXPTYPE Rtype = TYPEOF(nzvals);
	switch (Rtype) {
	    case REALSXP:
		return collect_double_infinite_nzoffs(REAL(nzvals),
					nzoffs, nzcount, out_nzoffs);
	    case CPLXSXP:
		return collect_Rcomplex_infinite_nzoffs(COMPLEX(nzvals),
					nzoffs, nzcount, out_nzoffs);
	}
	error("SparseArray internal error in collect_infinite_nzoffs():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return 0;  /* will never reach this */
}

static CollectnzoffsFUN get_CollectnzoffsFUN(const char *isFUN)
{
	if (strcmp(isFUN, "is.na") == 0)
		return collect_na_nzoffs;
	if (strcmp(isFUN, "is.nan") == 0)
		return collect_nan_nzoffs;
	if (strcmp(isFUN, "is.infinite") == 0)
		return collect_infinite_nzoffs;
	error("SparseArray internal error in get_CollectnzoffsFUN():\n"
	      "    unsupported function: \"%s\"", isFUN);
	return NULL;  /* will never reach this */
}


/*****************************************************************************
 * C_SVT_apply_isFUN()
 */

static SEXP leaf_apply_isFUN(SEXP leaf, CollectnzoffsFUN fun, int *nzoffs_buf)
{
	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
	if (nzvals == R_NilValue)  /* lacunar leaf */
		return R_NilValue;
	/* regular leaf */
	int ans_nzcount = fun(nzvals, INTEGER(nzoffs), nzcount, nzoffs_buf);
	if (ans_nzcount == 0)
		return R_NilValue;
	SEXP ans_nzoffs = PROTECT(NEW_INTEGER(ans_nzcount));
	memcpy(INTEGER(ans_nzoffs), nzoffs_buf, sizeof(int) * ans_nzcount);
	SEXP ans = _make_lacunar_leaf(ans_nzoffs);
	UNPROTECT(1);
	return ans;
}

/* Recursive tree traversal. */
static SEXP REC_SVT_apply_isFUN(SEXP SVT, const int *dim, int ndim,
		CollectnzoffsFUN fun, int *nzoffs_buf)
{
	if (SVT == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* 'SVT' is a leaf (i.e. 1D SVT). */
		return leaf_apply_isFUN(SVT, fun, nzoffs_buf);
	}

	/* 'SVT' is a list. */
	int ans_len = dim[ndim - 1];  /* same as 'LENGTH(SVT)' */
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	int is_empty = 1;
	for (int i = 0; i < ans_len; i++) {
		SEXP subSVT = VECTOR_ELT(SVT, i);
		SEXP ans_elt = REC_SVT_apply_isFUN(subSVT, dim, ndim - 1,
						   fun, nzoffs_buf);
		if (ans_elt != R_NilValue) {
			PROTECT(ans_elt);
			SET_VECTOR_ELT(ans, i, ans_elt);
			UNPROTECT(1);
			is_empty = 0;
		}
	}
	UNPROTECT(1);
	return is_empty ? R_NilValue : ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_SVT_apply_isFUN(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP isFUN)
{
	SEXPTYPE x_Rtype = _get_Rtype_from_Rstring(x_type);
	if (x_Rtype == 0)
		error("SparseArray internal error in C_SVT_apply_isFUN():\n"
		      "    invalid 'x_type' value");

	if (!IS_CHARACTER(isFUN) || LENGTH(isFUN) != 1)
		error("SparseArray internal error in C_SVT_apply_isFUN():\n"
		      "    'isFUN' must be a single string");
	isFUN = STRING_ELT(isFUN, 0);
	if (isFUN == NA_STRING)
		error("SparseArray internal error in C_SVT_apply_isFUN():\n"
		      "    'isFUN' cannot be NA");

	CollectnzoffsFUN fun = get_CollectnzoffsFUN(CHAR(isFUN));

	if (x_Rtype == VECSXP)
		error("%s() is not supported yet on SVT_SparseArray "
		      "objects of type \"list\"", CHAR(isFUN));

	if (x_Rtype == RAWSXP || (fun != collect_na_nzoffs &&
				  x_Rtype != REALSXP &&
				  x_Rtype != CPLXSXP))
		return R_NilValue;

	int *nzoffs_buf = (int *)
		R_alloc(INTEGER(x_dim)[0], sizeof(int));
	return REC_SVT_apply_isFUN(x_SVT, INTEGER(x_dim), LENGTH(x_dim),
				   fun, nzoffs_buf);
}

