/****************************************************************************
 *        crossprod(), tcrossprod(), and %*% of SparseArray objects         *
 ****************************************************************************/
#include "SparseMatrix_mult.h"

#include "Rvector_utils.h"
#include "leaf_vector_utils.h"
#include "SVT_SparseArray_class.h"

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

static void compute_sym_dotprods_with_finite_col(SEXP SVT, const double *col,
		int j, double *out, int out_nrow)
{
	double *out1, *out2;
	int i;
	SEXP subSVT;

	out1 = out + 1;
	out2 = out + out_nrow;
	for (i = j + 1; i < out_nrow; i++, out1++, out2 += out_nrow) {
		subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue)
			continue;
		*out1 = *out2 =
			_dotprod_leaf_vector_and_finite_col(subSVT, col);
	}
	return;
}

static void compute_sym_dotprods_with_right_lv(SEXP SVT, SEXP lv2,
		int j, double *out, int out_nrow)
{
	double *out1, *out2, dp;
	int i;
	SEXP subSVT;

	out1 = out + 1;
	out2 = out + out_nrow;
	for (i = j + 1; i < out_nrow; i++, out1++, out2 += out_nrow) {
		subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue) {
			dp = _dotprod0_leaf_vector(lv2);
		} else {
			dp = _dotprod_leaf_vectors(subSVT, lv2);
		}
		*out1 = *out2 = dp;
	}
	return;
}

static void compute_dotprods2_with_finite_Rcol(SEXP SVT, const double *col,
		double *out, int out_nrow)
{
	int i;
	SEXP subSVT;

	for (i = 0; i < out_nrow; i++, out++) {
		subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue)
			continue;
		*out = _dotprod_leaf_vector_and_finite_col(subSVT, col);
	}
	return;
}

static void compute_dotprods2_with_finite_Lcol(const double *col, SEXP SVT,
		double *out, int out_nrow, int out_ncol)
{
	int j;
	SEXP subSVT;

	for (j = 0; j < out_ncol; j++, out += out_nrow) {
		subSVT = VECTOR_ELT(SVT, j);
		if (subSVT == R_NilValue)
			continue;
		*out = _dotprod_leaf_vector_and_finite_col(subSVT, col);
	}
	return;
}

static void compute_dotprods2_with_Rlv(SEXP SVT, SEXP lv2,
		double *out, int out_nrow)
{
	int i;
	SEXP subSVT;

	for (i = 0; i < out_nrow; i++, out++) {
		subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue) {
			*out = _dotprod0_leaf_vector(lv2);
		} else {
			*out = _dotprod_leaf_vectors(subSVT, lv2);
		}
	}
	return;
}

static void compute_dotprods2_with_Llv(SEXP lv1, SEXP SVT,
		double *out, int out_nrow, int out_ncol)
{
	int j;
	SEXP subSVT;

	for (j = 0; j < out_ncol; j++, out += out_nrow) {
		subSVT = VECTOR_ELT(SVT, j);
		if (subSVT == R_NilValue) {
			*out = _dotprod0_leaf_vector(lv1);
		} else {
			*out = _dotprod_leaf_vectors(lv1, subSVT);
		}
	}
	return;
}

/* Cross-product between a fictional matrix of zeros and 'SVT2'. */
static void crossprod_Lzeros_double(SEXP SVT2,
				    double *out, int out_nrow, int out_ncol)
{
	int i, j;
	SEXP subSVT2;

	if (SVT2 == R_NilValue)
		return;
	for (j = 0; j < out_ncol; j++, out += out_nrow) {
		subSVT2 = VECTOR_ELT(SVT2, j);
		if (subSVT2 == R_NilValue)
			continue;
		for (i = 0; i < out_nrow; i++)
			out[i] = _dotprod0_leaf_vector(subSVT2);
	}
	return;
}

/* Cross-product between 'SVT1' and a fictional matrix of zeros. */
static void crossprod_Rzeros_double(SEXP SVT1,
				    double *out, int out_nrow, int out_ncol)
{
	int i, j;
	SEXP subSVT1;

	if (SVT1 == R_NilValue)
		return;
	for (i = 0; i < out_nrow; i++, out++) {
		subSVT1 = VECTOR_ELT(SVT1, i);
		if (subSVT1 == R_NilValue)
			continue;
		for (j = 0; j < out_ncol; j++)
			out[j * out_nrow] = _dotprod0_leaf_vector(subSVT1);
	}
	return;
}

static void crossprod1_double(SEXP SVT, int in_nrow, double *out, int out_ncol)
{
	int j;
	SEXP subSVT;
	double *col;

	if (SVT == R_NilValue)
		return;
	col = (double *) R_alloc(in_nrow, sizeof(double));
	for (j = 0; j < out_ncol; j++, out += out_ncol + 1) {
		subSVT = VECTOR_ELT(SVT, j);
		if (subSVT == R_NilValue) {
			memset(col, 0, sizeof(double) * in_nrow);
			compute_sym_dotprods_with_finite_col(SVT, col,
							j, out, out_ncol);
		} else if (lv_is_finite(subSVT)) {
			expand_double_lv(subSVT, col, in_nrow);
			*out =
			  _dotprod_leaf_vector_and_finite_col(subSVT, col);
			compute_sym_dotprods_with_finite_col(SVT, col,
							j, out, out_ncol);
		} else {
			*out = _dotprod_leaf_vectors(subSVT, subSVT);
			compute_sym_dotprods_with_right_lv(SVT, subSVT,
							j, out, out_ncol);
		}
	}
	return;
}

/* Fills 'out' column by column and performs "right preprocessing". That is,
   for each column in 'out', the corresponding leaf vector in 'SVT2' is
   turned into a dense representation, but only if it contains no infinite
   value (i.e. no NA, NaN, Inf, or -Inf). */
static void crossprod2_Rpp_double(SEXP SVT1, SEXP SVT2, int in_nrow,
				  double *out, int out_nrow, int out_ncol)
{
	int j;
	SEXP subSVT2;
	double *col;

	if (SVT1 == R_NilValue) {
		crossprod_Lzeros_double(SVT2, out, out_nrow, out_ncol);
		return;
	}
	col = (double *) R_alloc(in_nrow, sizeof(double));
	subSVT2 = R_NilValue;
	for (j = 0; j < out_ncol; j++, out += out_nrow) {
		if (SVT2 != R_NilValue)
			subSVT2 = VECTOR_ELT(SVT2, j);
		if (subSVT2 == R_NilValue) {
			memset(col, 0, sizeof(double) * in_nrow);
			compute_dotprods2_with_finite_Rcol(SVT1, col,
						out, out_nrow);
		} else if (lv_is_finite(subSVT2)) {
			/* Preprocess 'subSVT2'. */
			expand_double_lv(subSVT2, col, in_nrow);
			compute_dotprods2_with_finite_Rcol(SVT1, col,
						out, out_nrow);
		} else {
			compute_dotprods2_with_Rlv(SVT1, subSVT2,
						out, out_nrow);
		}
	}
	return;
}

/* Fills 'out' row by row and performs "left preprocessing". That is,
   for each row in 'out', the corresponding leaf vector in 'SVT1' is
   turned into a dense representation, but only if it contains no infinite
   value (i.e. no NA, NaN, Inf, or -Inf). */
static void crossprod2_Lpp_double(SEXP SVT1, SEXP SVT2, int in_nrow,
				  double *out, int out_nrow, int out_ncol)
{
	int i;
	SEXP subSVT1;
	double *col;

	if (SVT2 == R_NilValue) {
		crossprod_Rzeros_double(SVT1, out, out_nrow, out_ncol);
		return;
	}
	col = (double *) R_alloc(in_nrow, sizeof(double));
	subSVT1 = R_NilValue;
	for (i = 0; i < out_nrow; i++, out++) {
		if (SVT1 != R_NilValue)
			subSVT1 = VECTOR_ELT(SVT1, i);
		if (subSVT1 == R_NilValue) {
			memset(col, 0, sizeof(double) * in_nrow);
			compute_dotprods2_with_finite_Lcol(col, SVT2,
						out, out_nrow, out_ncol);
		} else if (lv_is_finite(subSVT1)) {
			/* Preprocess 'subSVT1'. */
			expand_double_lv(subSVT1, col, in_nrow);
			compute_dotprods2_with_finite_Lcol(col, SVT2,
						out, out_nrow, out_ncol);
		} else {
			compute_dotprods2_with_Llv(subSVT1, SVT2,
						out, out_nrow, out_ncol);
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

	crossprod1_double(x_SVT, x_nrow, REAL(ans), x_ncol);

	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_SVT_crossprod2(SEXP x_dim, SEXP x_SVT, SEXP y_dim, SEXP y_SVT,
		      SEXP ans_type, SEXP ans_dimnames)
{
	SEXPTYPE ans_Rtype;
	int in_nrow, x_ncol, y_ncol;
	SEXP ans;
	R_xlen_t Rpp_nops, Lpp_nops;

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
	in_nrow = INTEGER(x_dim)[0];
	if (in_nrow != INTEGER(y_dim)[0])
		error("non-conformable arguments");
	x_ncol = INTEGER(x_dim)[1];
	y_ncol = INTEGER(y_dim)[1];

	/* Allocate 'ans' and fill with zeros. */
	ans = PROTECT(new_Rmatrix(ans_Rtype, x_ncol, y_ncol, ans_dimnames));

	/* Calculate the nb of ops if doing right preprocessing vs
	   left preprocessing. */
	Rpp_nops = _REC_get_SVT_nzcount(x_SVT, LENGTH(x_dim)) * y_ncol;
	Lpp_nops = _REC_get_SVT_nzcount(y_SVT, LENGTH(y_dim)) * x_ncol;
	if (Rpp_nops <= Lpp_nops) {
		/* Will preprocess cols in 'y_SVT' (right preprocessing). */
		crossprod2_Rpp_double(x_SVT, y_SVT, in_nrow,
				      REAL(ans), x_ncol, y_ncol);
	} else {
		/* Will preprocess cols in 'x_SVT' (left preprocessing). */
		crossprod2_Lpp_double(x_SVT, y_SVT, in_nrow,
				      REAL(ans), x_ncol, y_ncol);
	}

	UNPROTECT(1);
	return ans;
}

