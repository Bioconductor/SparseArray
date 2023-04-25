/****************************************************************************
 *        crossprod(), tcrossprod(), and %*% of SparseArray objects         *
 ****************************************************************************/
#include "SparseMatrix_mult.h"

#include "Rvector_utils.h"
#include "leaf_vector_utils.h"
#include "SVT_SparseArray_class.h"

#include <string.h>  /* for memset() */

/* TODO: Maybe move this to Rvector_summarization.c */
static int is_finite_doubles(const double *x, int x_len)
{
	for (int i = 0; i < x_len; i++, x++) {
		/* ISNAN(): True for *both* NA and NaN. See <R_ext/Arith.h> */
		if (ISNAN(*x) || *x == R_PosInf || *x == R_NegInf)
			return 0;
	}
	return 1;
}

static int is_finite_lv(SEXP lv)
{
	int lv_len;
	SEXP lv_offs, lv_vals;

	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	return is_finite_doubles(REAL(lv_vals), lv_len);
}

static int noNA_ints(const int *x, int x_len)
{
	for (int i = 0; i < x_len; i++, x++) {
		if (*x == NA_INTEGER)
			return 0;
	}
	return 1;
}

static int noNA_lv(SEXP lv)
{
	int lv_len;
	SEXP lv_offs, lv_vals;

	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	return noNA_ints(INTEGER(lv_vals), lv_len);
}

static void sym_fill_with_NAs(double *out, int out_nrow, int j)
{
	double *out1, *out2;
	int i;

	*out = NA_REAL;
	out1 = out + 1;
	out2 = out + out_nrow;
	for (i = j + 1; i < out_nrow; i++, out1++, out2 += out_nrow)
		*out1 = *out2 = NA_REAL;
	return;
}

static void fill_col_with_NAs(double *out, int out_nrow)
{
	int i;

	for (i = 0; i < out_nrow; i++, out++)
		*out = NA_REAL;
	return;
}

static void fill_row_with_NAs(double *out, int out_nrow, int out_ncol)
{
	int j;

	for (j = 0; j < out_ncol; j++, out += out_nrow)
		*out = NA_REAL;
	return;
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

static void expand_int_lv(SEXP lv, int *x, int x_len)
{
	int lv_len;
	SEXP lv_offs, lv_vals;

	memset(x, 0, sizeof(int) * x_len);
	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	_copy_ints_to_offsets(INTEGER(lv_vals), INTEGER(lv_offs), lv_len, x);
	return;
}

static void compute_sym_dotprods_with_finite_col(SEXP SVT, const double *col,
		double *out, int out_nrow, int j)
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

static void compute_sym_dotprods_with_noNA_int_col(SEXP SVT, const int *col,
		double *out, int out_nrow, int j)
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
			_dotprod_leaf_vector_and_noNA_int_col(subSVT, col);
	}
	return;
}

static void compute_sym_dotprods_with_lv(SEXP SVT, SEXP lv2,
		double *out, int out_nrow, int j)
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

static void compute_dotprods2_with_double_Rcol(SEXP SVT, const double *col,
		int in_nrow,
		double *out, int out_nrow)
{
	int i;
	SEXP subSVT;

	if (is_finite_doubles(col, in_nrow)) {
		compute_dotprods2_with_finite_Rcol(SVT, col, out, out_nrow);
		return;
	}
	for (i = 0; i < out_nrow; i++, out++) {
		subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue) {
			*out = _dotprod0_double_col(col, in_nrow);
			continue;
		}
		*out = _dotprod_leaf_vector_and_double_col(subSVT,
							   col, in_nrow);
	}
	return;
}

static void compute_dotprods2_with_noNA_int_Rcol(SEXP SVT, const int *col,
		double *out, int out_nrow)
{
	int i;
	SEXP subSVT;

	for (i = 0; i < out_nrow; i++, out++) {
		subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue)
			continue;
		*out = _dotprod_leaf_vector_and_noNA_int_col(subSVT, col);
	}
	return;
}

static void compute_dotprods2_with_int_Rcol(SEXP SVT, const int *col,
		int in_nrow,
		double *out, int out_nrow)
{
	int i;
	SEXP subSVT;

	if (noNA_ints(col, in_nrow)) {
		compute_dotprods2_with_noNA_int_Rcol(SVT, col, out, out_nrow);
		return;
	}
	for (i = 0; i < out_nrow; i++, out++) {
		subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue) {
			*out = _dotprod0_int_col(col, in_nrow);
			continue;
		}
		*out = _dotprod_leaf_vector_and_int_col(subSVT,
							col, in_nrow);
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

static void compute_dotprods2_with_double_Lcol(const double *col, SEXP SVT,
		int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	int j;
	SEXP subSVT;

	if (is_finite_doubles(col, in_nrow)) {
		compute_dotprods2_with_finite_Lcol(col, SVT,
						   out, out_nrow, out_ncol);
		return;
	}
	for (j = 0; j < out_ncol; j++, out += out_nrow) {
		subSVT = VECTOR_ELT(SVT, j);
		if (subSVT == R_NilValue) {
			*out = _dotprod0_double_col(col, in_nrow);
			continue;
		}
		*out = _dotprod_leaf_vector_and_double_col(subSVT,
							   col, in_nrow);
	}
	return;
}

static void compute_dotprods2_with_noNA_int_Lcol(const int *col, SEXP SVT,
		double *out, int out_nrow, int out_ncol)
{
	int j;
	SEXP subSVT;

	for (j = 0; j < out_ncol; j++, out += out_nrow) {
		subSVT = VECTOR_ELT(SVT, j);
		if (subSVT == R_NilValue)
			continue;
		*out = _dotprod_leaf_vector_and_noNA_int_col(subSVT, col);
	}
	return;
}

static void compute_dotprods2_with_int_Lcol(const int *col, SEXP SVT,
		int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	int j;
	SEXP subSVT;

	if (noNA_ints(col, in_nrow)) {
		compute_dotprods2_with_noNA_int_Lcol(col, SVT,
					out, out_nrow, out_ncol);
		return;
	}
	for (j = 0; j < out_ncol; j++, out += out_nrow) {
		subSVT = VECTOR_ELT(SVT, j);
		if (subSVT == R_NilValue) {
			*out = _dotprod0_int_col(col, in_nrow);
			continue;
		}
		*out = _dotprod_leaf_vector_and_int_col(subSVT,
							col, in_nrow);
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


/****************************************************************************
 * Workhorses behind C_crossprod2_SVT_mat() and C_crossprod2_mat_SVT()
 */

/* 'out_ncol' is guaranteed to be the same as 'ncol(mat2)'.
   Computes the dot product of 'lv' with each column in 'mat2'. */
static void crossprod2_lv_mat_double(SEXP lv1, const double *mat2,
		int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	/* Compute the 'ncol(mat2)' dot products by walking on 'lv1'
	   only once, i.e. each value in 'lv1' is visited only once. */
/*
	int lv1_len, k1, j2;
	SEXP lv1_offs, lv1_vals;
	const int *offs1_p;
	const double *vals1_p, *v2_p;
	double v1;

	lv1_len = _split_leaf_vector(lv1, &lv1_offs, &lv1_vals);
	offs1_p = INTEGER(lv1_offs);
	vals1_p = REAL(lv1_vals);
	for (k1 = 0; k1 < lv1_len; k1++) {
		v1 = *vals1_p;
		v2_p = mat2 + *offs1_p;
		for (j2 = 0; j2 < out_ncol; j2++, v2_p += in_nrow)
			out[j2 * out_nrow] += v1 * *v2_p;
		offs1_p++;
		vals1_p++;
	}
*/
	/* Walk on the columns of 'mat2', and, for each column, compute its
	   dot product with 'lv1'. This means that each value in 'lv1' is
	   visited 'ncol(mat2)' times. For some reason, this seems to be
	   significantly faster than the above. */
	int j2;

	for (j2 = 0; j2 < out_ncol; j2++, out += out_nrow, mat2 += in_nrow)
		*out = _dotprod_leaf_vector_and_finite_col(lv1, mat2);
	return;
}

/* 'out_nrow' is guaranteed to be the same as 'ncol(mat1)'.
   Computes the dot product of 'lv' with each column in 'mat1'. */
static void crossprod2_mat_lv_double(const double *mat1, SEXP lv2,
		int in_nrow,
		double *out, int out_nrow)
{
	/* Compute the 'ncol(mat1)' dot products by walking on 'lv2'
	   only once, i.e. each value in 'lv2' is visited only once. */
/*
	int lv2_len, k2, j1;
	SEXP lv2_offs, lv2_vals;
	const int *offs2_p;
	const double *vals2_p, *v1_p;
	double v2;

	lv2_len = _split_leaf_vector(lv2, &lv2_offs, &lv2_vals);
	offs2_p = INTEGER(lv2_offs);
	vals2_p = REAL(lv2_vals);
	for (k2 = 0; k2 < lv2_len; k2++) {
		v2 = *vals2_p;
		v1_p = mat1 + *offs2_p;
		for (j1 = 0; j1 < out_nrow; j1++, v1_p += in_nrow)
			out[j1] += *v1_p * v2;
		offs2_p++;
		vals2_p++;
	}
*/
	/* Walk on the columns of 'mat1', and, for each column, compute its
	   dot product with 'lv2'. This means that each value in 'lv2' is
	   visited 'ncol(mat1)' times. For some reason, this seems to be
	   significantly faster than the above. */
	int j1;

	for (j1 = 0; j1 < out_nrow; j1++, out++, mat1 += in_nrow)
		*out = _dotprod_leaf_vector_and_finite_col(lv2, mat1);
	return;
}

/* Cross-product between 'SVT1' (on the left) and a dense matrix
   of doubles (on the right).
   'SVT1' must contain leaf vectors of type "double". */
static void crossprod2_SVT_mat_double(SEXP SVT1, const double *mat2,
		int tr_mat2, int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	//int i;
	//SEXP subSVT1;

	if (SVT1 == R_NilValue)
		return;
	/* Outer loop: walk on SVT1.
	   Inner loop: walk on the columns of 'mat2'.
	   Does NOT support 'tr_mat2' != 0. */
/*
	for (i = 0; i < out_nrow; i++, out++) {
		subSVT1 = VECTOR_ELT(SVT1, i);
		if (subSVT1 == R_NilValue)
			continue;
		crossprod2_lv_mat_double(subSVT1, REAL(mat2), in_nrow,
					 out, out_nrow, out_ncol);
	}
*/
	/* Outer loop: walk on the columns of 'mat2' (or on its rows
	   if 'tr_mat2' is true).
	   Inner loop: walk on SVT1.
	   For some mysterious reason, this is faster than the above. */
	double *colbuf;
	int i, j;

	if (tr_mat2) {
		colbuf = (double *) R_alloc(in_nrow, sizeof(double));
		for (j = 0; j < out_ncol; j++, out += out_nrow) {
			/* Copy 'mat2[j, ]' to colbuf. */
			for (i = 0; i < in_nrow; i++)
				colbuf[i] = mat2[i * out_ncol];
			compute_dotprods2_with_double_Rcol(SVT1, colbuf,
						in_nrow,
						out, out_nrow);
			mat2++;
		}
	} else {
		for (j = 0; j < out_ncol; j++, out += out_nrow) {
			compute_dotprods2_with_double_Rcol(SVT1, mat2,
						in_nrow,
						out, out_nrow);
			mat2 += in_nrow;
		}
	}
	return;
}

/* Cross-product between a dense matrix of doubles (on the left)
   and 'SVT2' (on the right).
   'SVT2' must contain leaf vectors of type "double". */
static void crossprod2_mat_SVT_double(const double *mat1, SEXP SVT2,
		int tr_mat1, int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	//int j;
	//SEXP subSVT2;

	if (SVT2 == R_NilValue)
		return;
	/* Outer loop: walk on SVT2.
	   Inner loop: walk on the columns of 'mat1'.
	   Does NOT support 'tr_mat1' != 0. */
/*
	for (j = 0; j < out_ncol; j++, out += out_nrow) {
		subSVT2 = VECTOR_ELT(SVT2, j);
		if (subSVT2 == R_NilValue)
			continue;
		crossprod2_mat_lv_double(mat1, subSVT2, in_nrow,
					 out, out_nrow);
	}
*/
	/* Outer loop: walk on the columns of 'mat1' (or on its rows
	   if 'tr_mat1' is true).
	   Inner loop: walk on SVT2.
	   For some mysterious reason, this is faster than the above. */
	double *colbuf;
	int i, j;

	if (tr_mat1) {
		colbuf = (double *) R_alloc(in_nrow, sizeof(double));
		for (i = 0; i < out_nrow; i++, out++) {
			/* Copy 'mat1[i, ]' to colbuf. */
			for (j = 0; j < in_nrow; j++)
				colbuf[j] = mat1[j * out_nrow];
			compute_dotprods2_with_double_Lcol(colbuf, SVT2,
						in_nrow,
						out, out_nrow, out_ncol);
			mat1++;
		}
	} else {
		for (i = 0; i < out_nrow; i++, out++) {
			compute_dotprods2_with_double_Lcol(mat1, SVT2,
						in_nrow,
						out, out_nrow, out_ncol);
			mat1 += in_nrow;
		}
	}
	return;
}

/* Cross-product between 'SVT1' (on the left) and a dense matrix
   of ints (on the right).
   'SVT1' must contain leaf vectors of type "integer". */
static void crossprod2_SVT_mat_int(SEXP SVT1, const int *mat2,
		int tr_mat2, int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	int *colbuf, i, j;

	if (SVT1 == R_NilValue)
		return;
	/* Outer loop: walk on the columns of 'mat2' (or on its rows
	   if 'tr_mat2' is true).
	   Inner loop: walk on SVT1. */
	if (tr_mat2) {
		colbuf = (int *) R_alloc(in_nrow, sizeof(int));
		for (j = 0; j < out_ncol; j++, out += out_nrow) {
			/* Copy 'mat2[j, ]' to colbuf. */
			for (i = 0; i < in_nrow; i++)
				colbuf[i] = mat2[i * out_ncol];
			compute_dotprods2_with_int_Rcol(SVT1, colbuf,
						in_nrow,
						out, out_nrow);
			mat2++;
		}
	} else {
		for (j = 0; j < out_ncol; j++, out += out_nrow) {
			compute_dotprods2_with_int_Rcol(SVT1, mat2,
						in_nrow,
						out, out_nrow);
			mat2 += in_nrow;
		}
	}
	return;
}

/* Cross-product between a dense matrix of ints (on the left)
   and 'SVT2' (on the right).
   'SVT2' must contain leaf vectors of type "integer". */
static void crossprod2_mat_SVT_int(const int *mat1, SEXP SVT2,
		int tr_mat1, int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	int *colbuf, i, j;

	if (SVT2 == R_NilValue)
		return;
	/* Outer loop: walk on the columns of 'mat1' (or on its rows
	   if 'tr_mat1' is true).
	   Inner loop: walk on SVT2. */
	if (tr_mat1) {
		colbuf = (int *) R_alloc(in_nrow, sizeof(int));
		for (i = 0; i < out_nrow; i++, out++) {
			/* Copy 'mat1[i, ]' to colbuf. */
			for (j = 0; j < in_nrow; j++)
				colbuf[j] = mat1[j * out_nrow];
			compute_dotprods2_with_int_Lcol(colbuf, SVT2,
						in_nrow,
						out, out_nrow, out_ncol);
			mat1++;
		}
	} else {
		for (i = 0; i < out_nrow; i++, out++) {
			compute_dotprods2_with_int_Lcol(mat1, SVT2,
						in_nrow,
						out, out_nrow, out_ncol);
			mat1 += in_nrow;
		}
	}
	return;
}


/****************************************************************************
 * Workhorses behind C_crossprod2_SVT_SVT()
 */

/* Cross-product between a fictive matrix of zeros (on the left)
   and 'SVT2' (on the right).
   Note that 'SVT2' must contain leaf vectors of type "double", because
   that's the only type supported by _dotprod0_leaf_vector() for now. */
static void crossprod2_mat0_SVT_double(SEXP SVT2,
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

/* Cross-product between 'SVT1' (on the left) and a fictive matrix
   of zeros (on the right).
   Note that 'SVT1' must contain leaf vectors of type "double", because
   that's the only type supported by _dotprod0_leaf_vector() for now. */
static void crossprod2_SVT_mat0_double(SEXP SVT1,
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

/* Cross-product between a fictive matrix of zeros (on the left)
   and 'SVT2' (on the right).
   'SVT2' must contain leaf vectors of type "integer". */
static void crossprod2_mat0_SVT_int(SEXP SVT2,
		double *out, int out_nrow, int out_ncol)
{
	int j;
	SEXP subSVT2;

	if (SVT2 == R_NilValue)
		return;
	for (j = 0; j < out_ncol; j++, out += out_nrow) {
		subSVT2 = VECTOR_ELT(SVT2, j);
		if (subSVT2 == R_NilValue || noNA_lv(subSVT2))
			continue;
		fill_col_with_NAs(out, out_nrow);
	}
	return;
}

/* Cross-product between 'SVT1' (on the left) and a fictive matrix
   of zeros (on the right).
   'SVT1' must contain leaf vectors of type "integer". */
static void crossprod2_SVT_mat0_int(SEXP SVT1,
		double *out, int out_nrow, int out_ncol)
{
	int i;
	SEXP subSVT1;

	if (SVT1 == R_NilValue)
		return;
	for (i = 0; i < out_nrow; i++, out++) {
		subSVT1 = VECTOR_ELT(SVT1, i);
		if (subSVT1 == R_NilValue || noNA_lv(subSVT1))
			continue;
		fill_row_with_NAs(out, out_nrow, out_ncol);
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
	double *colbuf;

	if (SVT1 == R_NilValue) {
		crossprod2_mat0_SVT_double(SVT2, out, out_nrow, out_ncol);
		return;
	}
	colbuf = (double *) R_alloc(in_nrow, sizeof(double));
	subSVT2 = R_NilValue;
	for (j = 0; j < out_ncol; j++, out += out_nrow) {
		if (SVT2 != R_NilValue)
			subSVT2 = VECTOR_ELT(SVT2, j);
		if (subSVT2 == R_NilValue) {
			memset(colbuf, 0, sizeof(double) * in_nrow);
			compute_dotprods2_with_finite_Rcol(SVT1, colbuf,
						out, out_nrow);
		} else if (is_finite_lv(subSVT2)) {
			/* Preprocess 'subSVT2'. */
			expand_double_lv(subSVT2, colbuf, in_nrow);
			compute_dotprods2_with_finite_Rcol(SVT1, colbuf,
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
	double *colbuf;

	if (SVT2 == R_NilValue) {
		crossprod2_SVT_mat0_double(SVT1, out, out_nrow, out_ncol);
		return;
	}
	colbuf = (double *) R_alloc(in_nrow, sizeof(double));
	subSVT1 = R_NilValue;
	for (i = 0; i < out_nrow; i++, out++) {
		if (SVT1 != R_NilValue)
			subSVT1 = VECTOR_ELT(SVT1, i);
		if (subSVT1 == R_NilValue) {
			memset(colbuf, 0, sizeof(double) * in_nrow);
			compute_dotprods2_with_finite_Lcol(colbuf, SVT2,
						out, out_nrow, out_ncol);
		} else if (is_finite_lv(subSVT1)) {
			/* Preprocess 'subSVT1'. */
			expand_double_lv(subSVT1, colbuf, in_nrow);
			compute_dotprods2_with_finite_Lcol(colbuf, SVT2,
						out, out_nrow, out_ncol);
		} else {
			compute_dotprods2_with_Llv(subSVT1, SVT2,
						out, out_nrow, out_ncol);
		}
	}
	return;
}

/* crossprod2_Rpp_int() and crossprod2_Lpp_int() below are the workhorses
   behind binary crossprod() when the input objects are SVT_SparseMatrix
   objects of type "integer". They work **natively** on the object. However,
   it turns out that going thru crossprod2_Rpp_double() or
   crossprod2_Lpp_double() by doing:

       crossprod(`type<-`(x, "double"), `type<-`(y, "double"))

   is about as fast, if not slightly faster. So only benefit of
   crossprod2_[R|L]pp_int() over crossprod2_[R|L]pp_double() is that they
   keep memory footprint as small as possible. */

static void crossprod2_Rpp_int(SEXP SVT1, SEXP SVT2, int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	int j;
	SEXP subSVT2;
	int *colbuf;

	if (SVT1 == R_NilValue) {
		crossprod2_mat0_SVT_int(SVT2, out, out_nrow, out_ncol);
		return;
	}
	colbuf = (int *) R_alloc(in_nrow, sizeof(int));
	subSVT2 = R_NilValue;
	for (j = 0; j < out_ncol; j++, out += out_nrow) {
		if (SVT2 != R_NilValue)
			subSVT2 = VECTOR_ELT(SVT2, j);
		if (subSVT2 == R_NilValue) {
			memset(colbuf, 0, sizeof(int) * in_nrow);
			compute_dotprods2_with_noNA_int_Rcol(SVT1, colbuf,
						out, out_nrow);
		} else if (noNA_lv(subSVT2)) {
			/* Preprocess 'subSVT2'. */
			expand_int_lv(subSVT2, colbuf, in_nrow);
			compute_dotprods2_with_noNA_int_Rcol(SVT1, colbuf,
						out, out_nrow);
		} else {
			fill_col_with_NAs(out, out_nrow);
		}
	}
	return;
}

static void crossprod2_Lpp_int(SEXP SVT1, SEXP SVT2, int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	int i;
	SEXP subSVT1;
	int *colbuf;

	if (SVT2 == R_NilValue) {
		crossprod2_SVT_mat0_int(SVT1, out, out_nrow, out_ncol);
		return;
	}
	colbuf = (int *) R_alloc(in_nrow, sizeof(int));
	subSVT1 = R_NilValue;
	for (i = 0; i < out_nrow; i++, out++) {
		if (SVT1 != R_NilValue)
			subSVT1 = VECTOR_ELT(SVT1, i);
		if (subSVT1 == R_NilValue) {
			memset(colbuf, 0, sizeof(int) * in_nrow);
			compute_dotprods2_with_noNA_int_Lcol(colbuf, SVT2,
						out, out_nrow, out_ncol);
		} else if (noNA_lv(subSVT1)) {
			/* Preprocess 'subSVT1'. */
			expand_int_lv(subSVT1, colbuf, in_nrow);
			compute_dotprods2_with_noNA_int_Lcol(colbuf, SVT2,
						out, out_nrow, out_ncol);
		} else {
			fill_row_with_NAs(out, out_nrow, out_ncol);
		}
	}
	return;
}


/****************************************************************************
 * Workhorses behind C_crossprod1_SVT()
 */

static void crossprod1_double(SEXP SVT, int in_nrow, double *out, int out_ncol)
{
	int j;
	SEXP subSVT;
	double *colbuf;

	if (SVT == R_NilValue)
		return;
	colbuf = (double *) R_alloc(in_nrow, sizeof(double));
	for (j = 0; j < out_ncol; j++, out += out_ncol + 1) {
		subSVT = VECTOR_ELT(SVT, j);
		if (subSVT == R_NilValue) {
			memset(colbuf, 0, sizeof(double) * in_nrow);
			compute_sym_dotprods_with_finite_col(SVT, colbuf,
						out, out_ncol, j);
		} else if (is_finite_lv(subSVT)) {
			expand_double_lv(subSVT, colbuf, in_nrow);
			*out =
			  _dotprod_leaf_vector_and_finite_col(subSVT, colbuf);
			compute_sym_dotprods_with_finite_col(SVT, colbuf,
						out, out_ncol, j);
		} else {
			*out = _dotprod_leaf_vectors(subSVT, subSVT);
			compute_sym_dotprods_with_lv(SVT, subSVT,
						out, out_ncol, j);
		}
	}
	return;
}

/* crossprod1_int() is the workhorse behind unary crossprod() when the
   input object is an SVT_SparseMatrix object of type "integer". It
   works **natively** on the object. However, it turns out that going thru
   crossprod1_double() by doing:

       crossprod(`type<-`(x, "double"))

   is about as fast, if not slightly faster. So only benefit of
   crossprod1_int() over crossprod1_double() is that it keeps memory
   footprint as small as possible. */
static void crossprod1_int(SEXP SVT, int in_nrow, double *out, int out_ncol)
{
	int j;
	SEXP subSVT;
	int *colbuf;

	if (SVT == R_NilValue)
		return;
	colbuf = (int *) R_alloc(in_nrow, sizeof(int));
	for (j = 0; j < out_ncol; j++, out += out_ncol + 1) {
		subSVT = VECTOR_ELT(SVT, j);
		if (subSVT == R_NilValue) {
			memset(colbuf, 0, sizeof(int) * in_nrow);
			compute_sym_dotprods_with_noNA_int_col(SVT, colbuf,
						out, out_ncol, j);
		} else if (noNA_lv(subSVT)) {
			expand_int_lv(subSVT, colbuf, in_nrow);
			*out =
			  _dotprod_leaf_vector_and_noNA_int_col(subSVT, colbuf);
			compute_sym_dotprods_with_noNA_int_col(SVT, colbuf,
						out, out_ncol, j);
		} else {
			sym_fill_with_NAs(out, out_ncol, j);
		}
	}
	return;
}


/****************************************************************************
 * C_crossprod2_SVT_mat(), C_crossprod2_mat_SVT(), C_crossprod2_SVT_SVT(),
 * and C_crossprod1_SVT()
 */

static SEXPTYPE get_and_check_input_Rtype(SEXP type, const char *what)
{
	SEXPTYPE Rtype;

	Rtype = _get_Rtype_from_Rstring(type);
	if (Rtype == 0)
		error("SparseArray internal error in "
		      "get_and_check_input_Rtype():\n"
		      "    invalid '%s' value", what);
	if (Rtype != REALSXP && Rtype != INTSXP)
		error("SparseArray internal error in "
		      "get_and_check_input_Rtype():\n"
		      "    input type \"%s\" is not supported yet",
		      type2char(Rtype));
	return Rtype;
}

/* --- .Call ENTRY POINT --- */
SEXP C_crossprod2_SVT_mat(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP y,
			  SEXP transpose_y,
			  SEXP ans_type, SEXP ans_dimnames)
{
	int tr_y, x_nrow, x_ncol, y_nrow, y_ncol, ans_ncol;
	SEXP y_dim, ans;
	SEXPTYPE x_Rtype, ans_Rtype;

	tr_y = LOGICAL(transpose_y)[0];

	/* Check 'x_dim' and 'y_dim'. */
	y_dim = GET_DIM(y);
	if (LENGTH(x_dim) != 2 || LENGTH(y_dim) != 2)
		error("input objects must have 2 dimensions");
	x_nrow = INTEGER(x_dim)[0];
	x_ncol = INTEGER(x_dim)[1];
	y_nrow = INTEGER(y_dim)[0];
	y_ncol = INTEGER(y_dim)[1];
	if (x_nrow != (tr_y ? y_ncol : y_nrow))
		error("input objects are non-conformable");

	/* Check 'x_type' and 'TYPEOF(y)'. */
	x_Rtype = get_and_check_input_Rtype(x_type, "x_type");
	if (x_Rtype != TYPEOF(y))
		error("SparseArray internal error in "
		      "C_crossprod2_SVT_mat():\n"
		      "    'x_Rtype != TYPEOF(y)' not supported yet");

	/* Check 'ans_type'. */
	ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (ans_Rtype == 0)
		error("SparseArray internal error in "
		      "C_crossprod2_SVT_mat():\n"
		      "    invalid 'ans_type' value");
	if (ans_Rtype != REALSXP)
		error("SparseArray internal error in "
		      "C_crossprod2_SVT_mat():\n"
		      "    output type \"%s\" is not supported yet",
		      type2char(ans_Rtype));

	/* Allocate 'ans' and fill with zeros. */
	ans_ncol = tr_y ? y_nrow : y_ncol;
	ans = PROTECT(_new_Rmatrix0(ans_Rtype, x_ncol, ans_ncol, ans_dimnames));

	if (x_Rtype == REALSXP) {
		crossprod2_SVT_mat_double(x_SVT, REAL(y), tr_y, x_nrow,
				REAL(ans), x_ncol, ans_ncol);
	} else {
		crossprod2_SVT_mat_int(x_SVT, INTEGER(y), tr_y, x_nrow,
				REAL(ans), x_ncol, ans_ncol);
	}

	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_crossprod2_mat_SVT(SEXP x, SEXP y_dim, SEXP y_type, SEXP y_SVT,
			  SEXP transpose_x,
			  SEXP ans_type, SEXP ans_dimnames)
{
	int tr_x, x_nrow, x_ncol, y_nrow, y_ncol, ans_nrow;
	SEXP x_dim, ans;
	SEXPTYPE y_Rtype, ans_Rtype;

	tr_x = LOGICAL(transpose_x)[0];

	/* Check 'x_dim' and 'y_dim'. */
	x_dim = GET_DIM(x);
	if (LENGTH(x_dim) != 2 || LENGTH(y_dim) != 2)
		error("input objects must have 2 dimensions");
	x_nrow = INTEGER(x_dim)[0];
	x_ncol = INTEGER(x_dim)[1];
	y_nrow = INTEGER(y_dim)[0];
	y_ncol = INTEGER(y_dim)[1];
	if ((tr_x ? x_ncol : x_nrow) != y_nrow)
		error("input objects are non-conformable");

	/* Check 'TYPEOF(x)' and 'y_type'. */
	y_Rtype = get_and_check_input_Rtype(y_type, "y_type");
	if (TYPEOF(x) != y_Rtype)
		error("input objects must have the same type() for now");

	/* Check 'ans_type'. */
	ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (ans_Rtype == 0)
		error("SparseArray internal error in "
		      "C_crossprod2_mat_SVT():\n"
		      "    invalid 'ans_type' value");
	if (ans_Rtype != REALSXP)
		error("SparseArray internal error in "
		      "C_crossprod2_mat_SVT():\n"
		      "    output type \"%s\" is not supported yet",
		      type2char(ans_Rtype));

	/* Allocate 'ans' and fill with zeros. */
	ans_nrow = tr_x ? x_nrow : x_ncol;
	ans = PROTECT(_new_Rmatrix0(ans_Rtype, ans_nrow, y_ncol, ans_dimnames));

	if (y_Rtype == REALSXP) {
		crossprod2_mat_SVT_double(REAL(x), y_SVT, tr_x, y_nrow,
				REAL(ans), ans_nrow, y_ncol);
	} else {
		crossprod2_mat_SVT_int(INTEGER(x), y_SVT, tr_x, y_nrow,
				REAL(ans), ans_nrow, y_ncol);
	}

	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_crossprod2_SVT_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT,
			  SEXP y_dim, SEXP y_type, SEXP y_SVT,
			  SEXP ans_type, SEXP ans_dimnames)
{
	int in_nrow, x_ncol, y_ncol;
	SEXPTYPE x_Rtype, y_Rtype, ans_Rtype;
	SEXP ans;
	R_xlen_t Rpp_nops, Lpp_nops;

	/* Check 'x_dim' and 'y_dim'. */
	if (LENGTH(x_dim) != 2 || LENGTH(y_dim) != 2)
		error("input objects must have 2 dimensions");
	in_nrow = INTEGER(x_dim)[0];
	if (in_nrow != INTEGER(y_dim)[0])
		error("input SVT_SparseMatrix objects "
		      "are non-conformable");
	x_ncol = INTEGER(x_dim)[1];
	y_ncol = INTEGER(y_dim)[1];

	/* Check 'x_type' and 'y_type'. */
	x_Rtype = get_and_check_input_Rtype(x_type, "x_type");
	y_Rtype = get_and_check_input_Rtype(y_type, "y_type");
	if (x_Rtype != y_Rtype)
		error("input SVT_SparseMatrix objects "
		      "must have the same type() for now");

	/* Check 'ans_type'. */
	ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (ans_Rtype == 0)
		error("SparseArray internal error in "
		      "C_crossprod2_SVT_SVT():\n"
		      "    invalid 'ans_type' value");
	if (ans_Rtype != REALSXP)
		error("SparseArray internal error in "
		      "C_crossprod2_SVT_SVT():\n"
		      "    output type \"%s\" is not supported yet",
		      type2char(ans_Rtype));

	/* Allocate 'ans' and fill with zeros. */
	ans = PROTECT(_new_Rmatrix0(ans_Rtype, x_ncol, y_ncol, ans_dimnames));

	/* Calculate the nb of ops if doing right preprocessing vs
	   doing left preprocessing. */
	Rpp_nops = _REC_get_SVT_nzcount(x_SVT, LENGTH(x_dim)) * y_ncol;
	Lpp_nops = _REC_get_SVT_nzcount(y_SVT, LENGTH(y_dim)) * x_ncol;
	if (Rpp_nops <= Lpp_nops) {
		/* Will preprocess cols in 'y_SVT' (right preprocessing). */
		if (x_Rtype == REALSXP) {
			crossprod2_Rpp_double(x_SVT, y_SVT, in_nrow,
					REAL(ans), x_ncol, y_ncol);
		} else {
			crossprod2_Rpp_int(x_SVT, y_SVT, in_nrow,
					REAL(ans), x_ncol, y_ncol);
		}
	} else {
		/* Will preprocess cols in 'x_SVT' (left preprocessing). */
		if (x_Rtype == REALSXP) {
			crossprod2_Lpp_double(x_SVT, y_SVT, in_nrow,
					REAL(ans), x_ncol, y_ncol);
		} else {
			crossprod2_Lpp_int(x_SVT, y_SVT, in_nrow,
					REAL(ans), x_ncol, y_ncol);
		}
	}

	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_crossprod1_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		      SEXP ans_type, SEXP ans_dimnames)
{
	int x_nrow, x_ncol;
	SEXPTYPE x_Rtype, ans_Rtype;
	SEXP ans;

	/* Check 'x_dim'. */
	if (LENGTH(x_dim) != 2)
		error("'x' must have 2 dimensions");
	x_nrow = INTEGER(x_dim)[0];
	x_ncol = INTEGER(x_dim)[1];

	/* Check 'x_type'. */
	x_Rtype = get_and_check_input_Rtype(x_type, "x_type");

	/* Check 'ans_type'. */
	ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (ans_Rtype == 0)
		error("SparseArray internal error in "
		      "C_crossprod1_SVT():\n"
		      "    invalid 'ans_type' value");
	if (ans_Rtype != REALSXP)
		error("SparseArray internal error in "
		      "C_crossprod1_SVT():\n"
		      "    output type \"%s\" is not supported yet",
		      type2char(ans_Rtype));

	/* Allocate 'ans' and fill with zeros. */
	ans = PROTECT(_new_Rmatrix0(ans_Rtype, x_ncol, x_ncol, ans_dimnames));

	if (x_Rtype == REALSXP) {
		crossprod1_double(x_SVT, x_nrow, REAL(ans), x_ncol);
	} else {
		crossprod1_int(x_SVT, x_nrow, REAL(ans), x_ncol);
	}

	UNPROTECT(1);
	return ans;
}

