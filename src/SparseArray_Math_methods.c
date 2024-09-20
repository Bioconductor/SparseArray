/*****************************************************************************
 *                   Math methods for SparseArray objects                    *
 *****************************************************************************/
#include "SparseArray_Math_methods.h"

#include "argcheck_utils.h"
#include "SparseVec.h"
#include "SparseVec_Math.h"
#include "leaf_utils.h"


static SEXP Math_leaf(MathFUN fun, SEXP leaf, double digits,
		SparseVec *buf_sv, int *newNaNs)
{
	const SparseVec sv = leaf2SV(leaf, REALSXP, buf_sv->len, 0);
	_Math_doubleSV(fun, &sv, digits, buf_sv, newNaNs);
	if (buf_sv->len == PROPAGATE_NZOFFS)
		return _make_leaf_with_single_shared_nzval(
					      buf_sv->Rtype, buf_sv->nzvals,
					      get_leaf_nzoffs(leaf));
	return SV2leaf(buf_sv);
}


/*****************************************************************************
 * Recursive tree traversal
 */

static SEXP REC_Math_SVT(MathFUN fun, SEXP SVT, double digits,
			 const int *dim, int ndim,
			 SparseVec *buf_sv, int *newNaNs)
{
	if (SVT == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* 'SVT' is a leaf (i.e. 1D SVT). */
		return Math_leaf(fun, SVT, digits, buf_sv, newNaNs);
	}

	/* 'SVT' is a list. */
	int ans_len = dim[ndim - 1];  /* same as 'LENGTH(SVT)' */
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	int is_empty = 1;
	for (int i = 0; i < ans_len; i++) {
		SEXP subSVT = VECTOR_ELT(SVT, i);
		SEXP ans_elt = REC_Math_SVT(fun, subSVT, digits,
					    dim, ndim - 1,
					    buf_sv, newNaNs);
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


/*****************************************************************************
 * C_Math_SVT()
 */

/* --- .Call ENTRY POINT --- */
SEXP C_Math_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP x_na_background,
		SEXP op, SEXP digits)
{
	/* Returned value ignored for now. */
	_get_and_check_Rtype_from_Rstring(x_type, "C_Math_SVT", "x_type");
	int x_has_NAbg = _get_and_check_na_background(x_na_background,
				"C_Math_SVT", "x_na_background");

	if (!IS_CHARACTER(op) || LENGTH(op) != 1)
		error("SparseArray internal error in C_Math_SVT():\n"
		      "    'op' must be a single string");
	op = STRING_ELT(op, 0);
	if (op == NA_STRING)
		error("SparseArray internal error in C_Math_SVT():\n"
		      "    'op' cannot be NA");

	MathFUN fun = _get_MathFUN(CHAR(op));
	double digits0 = REAL(digits)[0];

	int dim0 = INTEGER(x_dim)[0];
	SparseVec buf_sv = alloc_SparseVec(REALSXP, dim0, x_has_NAbg);

	int newNaNs = 0;
	SEXP ans = REC_Math_SVT(fun, x_SVT, digits0,
				INTEGER(x_dim), LENGTH(x_dim),
				&buf_sv, &newNaNs);
	if (newNaNs) {
		PROTECT(ans);
		warning("NaNs produced");
		UNPROTECT(1);
	}
	return ans;
}

