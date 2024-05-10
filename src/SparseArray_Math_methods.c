/****************************************************************************
 *                   Math methods for SparseArray objects                   *
 ****************************************************************************/
#include "SparseArray_Math_methods.h"

#include "Rvector_utils.h"
#include "SparseVec.h"
#include "SparseVec_Math.h"
#include "leaf_utils.h"


static SEXP Math_leaf(MathFUN fun, SEXP leaf, double digits, int dim0,
		double *nzvals_buf, int *nzoffs_buf, int *newNaNs)
{
	const SparseVec sv = leaf2SV(leaf, REALSXP, dim0);
	int buf_len = _Math_doubleSV(fun, &sv, digits,
				     nzvals_buf, nzoffs_buf, newNaNs);
	return _make_leaf_from_two_arrays(REALSXP,
					  nzvals_buf, nzoffs_buf, buf_len);
}


/****************************************************************************
 * Recursive tree traversal
 */

static SEXP REC_Math_SVT(MathFUN fun, SEXP SVT, double digits,
			 const int *dim, int ndim,
			 double *nzvals_buf, int *nzoffs_buf, int *newNaNs)
{
	if (SVT == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* 'SVT' is a leaf. */
		return Math_leaf(fun, SVT, digits, dim[0],
				 nzvals_buf, nzoffs_buf, newNaNs);
	}

	/* 'SVT' is a list. */
	int ans_len = dim[ndim - 1];
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	int is_empty = 1;
	for (int i = 0; i < ans_len; i++) {
		SEXP subSVT = VECTOR_ELT(SVT, i);
		SEXP ans_elt = REC_Math_SVT(fun, subSVT, digits,
					    dim, ndim - 1,
					    nzvals_buf, nzoffs_buf, newNaNs);
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


/****************************************************************************
 * C_Math_SVT()
 */

/* --- .Call ENTRY POINT --- */
SEXP C_Math_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP op, SEXP digits)
{
	SEXPTYPE x_Rtype = _get_Rtype_from_Rstring(x_type);
	if (x_Rtype == 0)
		error("SparseArray internal error in C_Math_SVT():\n"
		      "    invalid 'x_type' value");

	if (!IS_CHARACTER(op) || LENGTH(op) != 1)
		error("SparseArray internal error in C_Math_SVT():\n"
		      "    'op' must be a single string");
	op = STRING_ELT(op, 0);
	if (op == NA_STRING)
		error("SparseArray internal error in C_Math_SVT():\n"
		      "    'op' cannot be NA");

	MathFUN fun = _get_MathFUN(CHAR(op));
	double digits0 = REAL(digits)[0];

	double *nzvals_buf = (double *)
		R_alloc(INTEGER(x_dim)[0], sizeof(double));
	int *nzoffs_buf = (int *)
		R_alloc(INTEGER(x_dim)[0], sizeof(int));
	int newNaNs = 0;
	SEXP ans = REC_Math_SVT(fun, x_SVT, digits0,
				INTEGER(x_dim), LENGTH(x_dim),
				nzvals_buf, nzoffs_buf, &newNaNs);
	if (newNaNs) {
		PROTECT(ans);
		warning("NaNs produced");
		UNPROTECT(1);
	}
	return ans;
}

