/****************************************************************************
 *        crossprod(), tcrossprod(), and %*% of SparseArray objects         *
 ****************************************************************************/
#include "SparseMatrix_mult.h"

#include "Rvector_utils.h"
#include "leaf_vector_utils.h"

#include <string.h>  /* for memset() */

static int lv_is_finite(SEXP lv)
{
	int lv_len, k;
	SEXP lv_offs, lv_vals;
	const double *vals_p;
	double v;

	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	vals_p = REAL(lv_vals);
	for (k = 0; k < lv_len; k++) {
		v = *vals_p;
		/* ISNAN(): True for *both* NA and NaN. See <R_ext/Arith.h> */
		if (ISNAN(v) || v == R_PosInf || v == R_NegInf)
			return 0;
		vals_p++;
	}
	return 1;
}

static void expand_double_lv(SEXP lv, double *x, int x_len)
{
	int lv_len;
	SEXP lv_offs, lv_vals;

	memset(x, 0, sizeof(double) * x_len);
	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	_copy_doubles_to_offsets(REAL(lv_vals), INTEGER(lv_offs), lv_len, x);
	return;
}

static void compute_dotprods_with_finite_col(SEXP SVT, int ncol, int j,
		const double *col, double *out)
{
	double *out1, *out2;
	int i;
	SEXP subSVT1;

	out1 = out + 1;
	out2 = out + ncol;
	for (i = j + 1; i < ncol; i++, out1++, out2 += ncol) {
		subSVT1 = VECTOR_ELT(SVT, i);
		if (subSVT1 == R_NilValue)
			continue;
		*out1 = *out2 =
			_dotprod_leaf_vector_and_finite_col(subSVT1, col);
	}
	return;
}

static void compute_dotprods2_with_finite_col(SEXP SVT, int nrow,
		const double *col, double *out)
{
	int i;
	SEXP subSVT1;

	for (i = 0; i < nrow; i++, out++) {
		subSVT1 = VECTOR_ELT(SVT, i);
		if (subSVT1 == R_NilValue)
			continue;
		*out = _dotprod_leaf_vector_and_finite_col(subSVT1, col);
	}
	return;
}

static void compute_dotprods_with_right_lv(SEXP SVT, int ncol, int j,
		SEXP lv2, double *out)
{
	double *out1, *out2, dp;
	int i;
	SEXP subSVT1;

	out1 = out + 1;
	out2 = out + ncol;
	for (i = j + 1; i < ncol; i++, out1++, out2 += ncol) {
		subSVT1 = VECTOR_ELT(SVT, i);
		if (subSVT1 == R_NilValue) {
			dp = _dotprod0_leaf_vector(lv2);
		} else {
			dp = _dotprod_leaf_vectors(subSVT1, lv2);
		}
		*out1 = *out2 = dp;
	}
	return;
}

static void compute_dotprods2_with_right_lv(SEXP SVT, int nrow,
		SEXP lv2, double *out)
{
	int i;
	SEXP subSVT1;

	for (i = 0; i < nrow; i++, out++) {
		subSVT1 = VECTOR_ELT(SVT, i);
		if (subSVT1 == R_NilValue) {
			*out = _dotprod0_leaf_vector(lv2);
		} else {
			*out = _dotprod_leaf_vectors(subSVT1, lv2);
		}
	}
	return;
}

static void crossprod1_double(SEXP x_SVT, int x_nrow, int x_ncol, double *out)
{
	int j;
	SEXP subSVT2;
	double *col;

	if (x_SVT == R_NilValue)
		return;
	col = (double *) R_alloc(x_nrow, sizeof(double));
	for (j = 0; j < x_ncol; j++, out += x_ncol + 1) {
		subSVT2 = VECTOR_ELT(x_SVT, j);
		if (subSVT2 == R_NilValue) {
			memset(col, 0, sizeof(double) * x_nrow);
			compute_dotprods_with_finite_col(x_SVT, x_ncol, j,
							 col, out);
		} else if (lv_is_finite(subSVT2)) {
			expand_double_lv(subSVT2, col, x_nrow);
			*out =
			  _dotprod_leaf_vector_and_finite_col(subSVT2, col);
			compute_dotprods_with_finite_col(x_SVT, x_ncol, j,
							 col, out);
		} else {
			*out = _dotprod_leaf_vectors(subSVT2, subSVT2);
			compute_dotprods_with_right_lv(x_SVT, x_ncol, j,
						       subSVT2, out);
		}
	}
	return;
}

static void crossprod2_double(SEXP x_SVT, int x_nrow, int x_ncol,
			      SEXP y_SVT, int y_ncol, double *out)
{
	int i, j;
	SEXP subSVT2;
	double *col;

	if (x_SVT == R_NilValue) {
		if (y_SVT == R_NilValue)
			return;
		for (j = 0; j < y_ncol; j++, out += x_ncol) {
			subSVT2 = VECTOR_ELT(y_SVT, j);
			if (subSVT2 == R_NilValue)
				continue;
			for (i = 0; i < x_ncol; i++)
				out[i] = _dotprod0_leaf_vector(subSVT2);
		}
	}
	col = (double *) R_alloc(x_nrow, sizeof(double));
	subSVT2 = R_NilValue;
	for (j = 0; j < y_ncol; j++, out += x_ncol) {
		if (y_SVT != R_NilValue)
			subSVT2 = VECTOR_ELT(y_SVT, j);
		if (subSVT2 == R_NilValue) {
			memset(col, 0, sizeof(double) * x_nrow);
			compute_dotprods2_with_finite_col(x_SVT, x_ncol,
							  col, out);
		} else if (lv_is_finite(subSVT2)) {
			expand_double_lv(subSVT2, col, x_nrow);
			compute_dotprods2_with_finite_col(x_SVT, x_ncol,
							  col, out);
		} else {
			compute_dotprods2_with_right_lv(x_SVT, x_ncol,
						        subSVT2, out);
		}
	}
	return;
}

/****************************************************************************
 * C_SVT_crossprod1() and C_SVT_crossprod2()
 */

/* Like allocMatrix() but with initialization of the matrix elements.
   Also set the dimnames. */
static SEXP new_Rmatrix(SEXPTYPE Rtype, int nrow, int ncol, SEXP dimnames)
{
	SEXP ans;
	size_t Rtype_size;

	ans = PROTECT(allocMatrix(Rtype, nrow, ncol));
	SET_DIMNAMES(ans, dimnames);
	/* allocMatrix() is just a thin wrapper around allocVector() and
	   the latter does NOT initialize the vector elements, except for
	   a list or a character vector. */
	if (Rtype != STRSXP && Rtype != VECSXP) {
		Rtype_size = _get_Rtype_size(Rtype);
		if (Rtype_size == 0) {
			UNPROTECT(1);
			error("SparseArray internal error in new_Rmatrix():\n"
			      "    type \"%s\" is not supported",
			      type2char(Rtype));
		}
		memset(DATAPTR(ans), 0, Rtype_size * XLENGTH(ans));
	}
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_SVT_crossprod1(SEXP x_dim, SEXP x_SVT, SEXP ans_type, SEXP ans_dimnames)
{
	SEXPTYPE ans_Rtype;
	int x_nrow, x_ncol;
	SEXP ans;

	/* Check 'ans_type'. */
	ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (ans_Rtype == 0)
		error("SparseArray internal error in "
		      "C_SVT_crossprod1():\n"
		      "    invalid 'ans_type' value");
	if (ans_Rtype != REALSXP)
		error("SparseArray internal error in "
		      "C_SVT_crossprod1():\n"
		      "    type \"%s\" is not supported yet",
		      type2char(ans_Rtype));

	/* Check 'x_dim'. */
	if (LENGTH(x_dim) != 2)
		error("'x' must have 2 dimensions");
	x_nrow = INTEGER(x_dim)[0];
	x_ncol = INTEGER(x_dim)[1];

	/* Allocate 'ans' and fill with zeros. */
	ans = PROTECT(new_Rmatrix(ans_Rtype, x_ncol, x_ncol, ans_dimnames));

	crossprod1_double(x_SVT, x_nrow, x_ncol, REAL(ans));

	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_SVT_crossprod2(SEXP x_dim, SEXP x_SVT, SEXP y_dim, SEXP y_SVT,
		      SEXP ans_type, SEXP ans_dimnames)
{
	SEXPTYPE ans_Rtype;
	int x_nrow, x_ncol, y_ncol;
	SEXP ans;

	/* Check 'ans_type'. */
	ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (ans_Rtype == 0)
		error("SparseArray internal error in "
		      "C_SVT_crossprod2():\n"
		      "    invalid 'ans_type' value");
	if (ans_Rtype != REALSXP)
		error("SparseArray internal error in "
		      "C_SVT_crossprod2():\n"
		      "    type \"%s\" is not supported yet",
		      type2char(ans_Rtype));

	/* Check 'x_dim' and 'y_dim'. */
	if (LENGTH(x_dim) != 2 || LENGTH(y_dim) != 2)
		error("arguments must have 2 dimensions");
	x_nrow = INTEGER(x_dim)[0];
	if (x_nrow != INTEGER(y_dim)[0])
		error("non-conformable arguments");
	x_ncol = INTEGER(x_dim)[1];
	y_ncol = INTEGER(y_dim)[1];

	/* Allocate 'ans' and fill with zeros. */
	ans = PROTECT(new_Rmatrix(ans_Rtype, x_ncol, y_ncol, ans_dimnames));

	crossprod2_double(x_SVT, x_nrow, x_ncol, y_SVT, y_ncol, REAL(ans));

	UNPROTECT(1);
	return ans;
}

