/****************************************************************************
 *        crossprod(), tcrossprod(), and %*% of SparseMatrix objects        *
 ****************************************************************************/
#include "SparseMatrix_mult.h"

#include "Rvector_utils.h"
#include "SparseVec.h"
#include "SparseVec_dotprod.h"
#include "leaf_utils.h"             /* for leaf2SV() */
#include "SVT_SparseArray_class.h"  /* for _REC_nzcount_SVT() */

#include <string.h>  /* for memset() */


/* TODO: Maybe move this to Rvector_summarization.c */
static int has_no_NaN_or_Inf(const double *x, int x_len)
{
	for (int i = 0; i < x_len; i++, x++)
		if (!R_FINITE(*x)) return 0;
	return 1;
}

/* TODO: Maybe move this to Rvector_summarization.c */
static int has_no_NA(const int *x, int x_len)
{
	for (int i = 0; i < x_len; i++, x++)
		if (*x == NA_INTEGER) return 0;
	return 1;
}

static int doubleSV_has_no_NaN_or_Inf(const SparseVec *sv)
{
	return has_no_NaN_or_Inf(get_doubleSV_nzvals(sv),
				 get_SV_nzcount(sv));
}

static int intSV_has_no_NA(const SparseVec *sv)
{
	return has_no_NA(get_intSV_nzvals(sv),
			 get_SV_nzcount(sv));
}

static void fill_col(double *out, int out_nrow, double v)
{
	int i;

	for (i = 0; i < out_nrow; i++, out++)
		*out = v;
	return;
}

static void fill_row(double *out, int out_nrow, int out_ncol, double v)
{
	int j;

	for (j = 0; j < out_ncol; j++, out += out_nrow)
		*out = v;
	return;
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

static void expand_doubleSV(const SparseVec *sv, double *out)
{
	memset(out, 0, sizeof(double) * sv->len);
	_copy_doubles_to_offsets(get_doubleSV_nzvals(sv),
				 sv->nzoffs, get_SV_nzcount(sv), out);
	return;
}

static void expand_intSV(const SparseVec *sv, int *out)
{
	memset(out, 0, sizeof(int) * sv->len);
	_copy_ints_to_offsets(get_intSV_nzvals(sv),
			      sv->nzoffs, get_SV_nzcount(sv), out);
	return;
}

static double dotprod_leaf_finite_doubles(SEXP leaf1,
					  const double *x2, int x2_len)
{
	if (leaf1 == R_NilValue)
		return 0.0;
	const SparseVec sv1 = leaf2SV(leaf1, x2_len);
	return _dotprod_doubleSV_finite_doubles(&sv1, x2);
}

static double dotprod_leaf_noNA_ints(SEXP leaf1,
				     const int *x2, int x2_len)
{
	if (leaf1 == R_NilValue)
		return 0.0;
	const SparseVec sv1 = leaf2SV(leaf1, x2_len);
	return _dotprod_intSV_noNA_ints(&sv1, x2);
}

static double dotprod_leaf_doubles(SEXP leaf1, const double *x2, int x2_len)
{
	if (leaf1 == R_NilValue)
		return _dotprod_doubles_zero(x2, x2_len);
	const SparseVec sv1 = leaf2SV(leaf1, x2_len);
	return _dotprod_doubleSV_doubles(&sv1, x2);
}

static double dotprod_leaf_ints(SEXP leaf1, const int *x2, int x2_len)
{
	if (leaf1 == R_NilValue)
		return _dotprod_ints_zero(x2, x2_len);
	const SparseVec sv1 = leaf2SV(leaf1, x2_len);
	return _dotprod_intSV_ints(&sv1, x2);
}

static double dotprod_leaf_doubleSV(SEXP leaf1, const SparseVec *sv2)
{
	if (leaf1 == R_NilValue)
		return _dotprod_doubleSV_zero(sv2);
	const SparseVec sv1 = leaf2SV(leaf1, sv2->len);
	return _dotprod_doubleSV_doubleSV(&sv1, sv2);
}


/****************************************************************************
 * Core multithreaded routines
 *
 * Note that all the functions in this section take for granted that they get
 * passed a 'SVT' (Sparse Vector Tree) that is **not** R_NilValue so that
 * better be true.
 */

static void compute_dotprods2_with_finite_Lcol(const double *Lcol, int Lcol_len,
		SEXP SVT, double *out, int out_nrow, int out_ncol)
{
	#pragma omp parallel for schedule(static)
	for (int j = 0; j < out_ncol; j++) {
		SEXP leaf = VECTOR_ELT(SVT, j);
		double dp = dotprod_leaf_finite_doubles(leaf, Lcol, Lcol_len);
		out[j * out_nrow] = dp;
	}
	return;
}

static void compute_dotprods2_with_finite_Rcol(SEXP SVT,
		const double *Rcol, int Rcol_len, double *out, int out_nrow)
{
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < out_nrow; i++) {
		SEXP leaf = VECTOR_ELT(SVT, i);
		out[i] = dotprod_leaf_finite_doubles(leaf, Rcol, Rcol_len);
	}
	return;
}

static void compute_dotprods2_with_noNA_int_Lcol(const int *Lcol, int Lcol_len,
		SEXP SVT, double *out, int out_nrow, int out_ncol)
{
	#pragma omp parallel for schedule(static)
	for (int j = 0; j < out_ncol; j++) {
		SEXP leaf = VECTOR_ELT(SVT, j);
		double dp = dotprod_leaf_noNA_ints(leaf, Lcol, Lcol_len);
		out[j * out_nrow] = dp;
	}
	return;
}

static void compute_dotprods2_with_noNA_int_Rcol(SEXP SVT,
		const int *Rcol, int Rcol_len, double *out, int out_nrow)
{
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < out_nrow; i++) {
		SEXP leaf = VECTOR_ELT(SVT, i);
		out[i] = dotprod_leaf_noNA_ints(leaf, Rcol, Rcol_len);
	}
	return;
}

static void compute_dotprods2_with_double_Lcol(const double *Lcol, int Lcol_len,
		SEXP SVT, double *out, int out_nrow, int out_ncol)
{
	if (has_no_NaN_or_Inf(Lcol, Lcol_len)) {
		compute_dotprods2_with_finite_Lcol(Lcol, Lcol_len, SVT,
						   out, out_nrow, out_ncol);
		return;
	}
	#pragma omp parallel for schedule(static)
	for (int j = 0; j < out_ncol; j++) {
		SEXP leaf = VECTOR_ELT(SVT, j);
		out[j * out_nrow] = dotprod_leaf_doubles(leaf, Lcol, Lcol_len);
	}
	return;
}

static void compute_dotprods2_with_double_Rcol(SEXP SVT,
		const double *Rcol, int Rcol_len, double *out, int out_nrow)
{
	if (has_no_NaN_or_Inf(Rcol, Rcol_len)) {
		compute_dotprods2_with_finite_Rcol(SVT, Rcol, Rcol_len,
						   out, out_nrow);
		return;
	}
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < out_nrow; i++) {
		SEXP leaf = VECTOR_ELT(SVT, i);
		out[i] = dotprod_leaf_doubles(leaf, Rcol, Rcol_len);
	}
	return;
}

static void compute_dotprods2_with_int_Lcol(const int *Lcol, int Lcol_len,
		SEXP SVT, double *out, int out_nrow, int out_ncol)
{
	if (has_no_NA(Lcol, Lcol_len)) {
		compute_dotprods2_with_noNA_int_Lcol(Lcol, Lcol_len, SVT,
						     out, out_nrow, out_ncol);
		return;
	}
	#pragma omp parallel for schedule(static)
	for (int j = 0; j < out_ncol; j++) {
		SEXP leaf = VECTOR_ELT(SVT, j);
		out[j * out_nrow] = dotprod_leaf_ints(leaf, Lcol, Lcol_len);
	}
	return;
}

static void compute_dotprods2_with_int_Rcol(SEXP SVT,
		const int *Rcol, int Rcol_len, double *out, int out_nrow)
{
	if (has_no_NA(Rcol, Rcol_len)) {
		compute_dotprods2_with_noNA_int_Rcol(SVT, Rcol, Rcol_len,
						     out, out_nrow);
		return;
	}
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < out_nrow; i++) {
		SEXP leaf = VECTOR_ELT(SVT, i);
		out[i] = dotprod_leaf_ints(leaf, Rcol, Rcol_len);
	}
	return;
}

static void compute_dotprods2_with_Lsv(const SparseVec *sv1, SEXP SVT,
		double *out, int out_nrow, int out_ncol)
{
	#pragma omp parallel for schedule(static)
	for (int j = 0; j < out_ncol; j++) {
		SEXP leaf = VECTOR_ELT(SVT, j);
		out[j * out_nrow] = dotprod_leaf_doubleSV(leaf, sv1);
	}
	return;
}

static void compute_dotprods2_with_Rsv(SEXP SVT, const SparseVec *sv2,
		double *out, int out_nrow)
{
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < out_nrow; i++) {
		SEXP leaf = VECTOR_ELT(SVT, i);
		out[i] = dotprod_leaf_doubleSV(leaf, sv2);
	}
	return;
}

static void compute_sym_dotprods_with_finite_col(SEXP SVT, int j,
		const double *col, int col_len, double *out, int out_nrow)
{
	#pragma omp parallel for schedule(static)
	for (int k = out_nrow - 1 - j; k >= 1; k--) {
		SEXP leaf = VECTOR_ELT(SVT, j + k);
		out[k] = out[k * out_nrow] =
			dotprod_leaf_finite_doubles(leaf, col, col_len);
	}
	return;
}

static void compute_sym_dotprods_with_noNA_int_col(SEXP SVT, int j,
		const int *col, int col_len, double *out, int out_nrow)
{
	#pragma omp parallel for schedule(static)
	for (int k = out_nrow - 1 - j; k >= 1; k--) {
		SEXP leaf = VECTOR_ELT(SVT, j + k);
		out[k] = out[k * out_nrow] =
			dotprod_leaf_noNA_ints(leaf, col, col_len);
	}
	return;
}

static void compute_sym_dotprods_with_doubleSV(SEXP SVT, int j,
		const SparseVec *sv2, double *out, int out_nrow)
{
	#pragma omp parallel for schedule(static)
	for (int k = out_nrow - 1 - j; k >= 1; k--) {
		SEXP leaf = VECTOR_ELT(SVT, j + k);
		out[k] = out[k * out_nrow] = dotprod_leaf_doubleSV(leaf, sv2);
	}
	return;
}


/****************************************************************************
 * Workhorses behind C_crossprod2_SVT_mat() and C_crossprod2_mat_SVT()
 */

/* 'out_ncol' is guaranteed to be the same as 'ncol(mat2)'.
   Computes the dot product of 'sv1' with each column in 'mat2'.
   CURRENTLY UNUSED */
static void crossprod2_doubleSV_doublemat(
		const SparseVec *sv1, const double *mat2,
		int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	/* Compute the 'ncol(mat2)' dot products by walking on 'sv1'
	   only once, i.e. each value in 'sv1' is visited only once. */
/*
	int lv1_len, k1, j2;
	SEXP lv1_offs, lv1_vals;
	const int *offs1_p;
	const double *vals1_p, *v2_p;
	double v1;

	lv1_len = unzip_leaf(lv1, &lv1_offs, &lv1_vals);
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
	   dot product with 'sv1'. This means that each value in 'sv1' is
	   visited 'ncol(mat2)' times. For some reason, this seems to be
	   significantly faster than the above. */
	int j2;

	for (j2 = 0; j2 < out_ncol; j2++, out += out_nrow, mat2 += in_nrow)
		*out = _dotprod_doubleSV_finite_doubles(sv1, mat2);
	return;
}

/* 'out_nrow' is guaranteed to be the same as 'ncol(mat1)'.
   Computes the dot product of each column in 'mat1' with 'sv2'.
   CURRENTLY UNUSED */
static void crossprod2_doublemat_doubleSV(
		const double *mat1, const SparseVec *sv2,
		int in_nrow,
		double *out, int out_nrow)
{
	/* Compute the 'ncol(mat1)' dot products by walking on 'sv2'
	   only once, i.e. each value in 'sv2' is visited only once. */
/*
	int lv2_len, k2, j1;
	SEXP lv2_offs, lv2_vals;
	const int *offs2_p;
	const double *vals2_p, *v1_p;
	double v2;

	lv2_len = unzip_leaf(lv2, &lv2_offs, &lv2_vals);
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
	   dot product with 'sv2'. This means that each value in 'sv2' is
	   visited 'ncol(mat1)' times. For some reason, this seems to be
	   significantly faster than the above. */
	int j1;

	for (j1 = 0; j1 < out_nrow; j1++, out++, mat1 += in_nrow)
		*out = _dotprod_doubleSV_finite_doubles(sv2, mat1);
	return;
}

/* Cross-product between 'SVT1' (on the left) and a dense matrix
   of doubles (on the right).
   'SVT1' must contain leaf vectors of type "double". */
static void crossprod2_SVT_mat_double(SEXP SVT1, const double *mat2,
		int tr_mat2, int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	if (SVT1 == R_NilValue)
		return;
	/* Outer loop: walk on SVT1.
	   Inner loop: walk on the columns of 'mat2'.
	   Does NOT support 'tr_mat2' != 0. */
/*
	for (int i = 0; i < out_nrow; i++, out++) {
		SEXP leaf = VECTOR_ELT(SVT1, i);
		if (leaf == R_NilValue)
			continue;
		const SparseVec sv1 = leaf2SV(leaf, in_nrow);
		crossprod2_doubleSV_doublemat(&sv1, REAL(mat2), in_nrow,
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
			/* Copy 'mat2[j, ]' to 'colbuf'. */
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
	if (SVT2 == R_NilValue)
		return;
	/* Outer loop: walk on SVT2.
	   Inner loop: walk on the columns of 'mat1'.
	   Does NOT support 'tr_mat1' != 0. */
/*
	for (int j = 0; j < out_ncol; j++, out += out_nrow) {
		SEXP leaf = VECTOR_ELT(SVT2, j);
		if (leaf == R_NilValue)
			continue;
		const SparseVec sv2 = leaf2SV(leaf, in_nrow);
		crossprod2_doublemat_doubleSV(mat1, &sv2, in_nrow,
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
			/* Copy 'mat1[i, ]' to 'colbuf'. */
			for (j = 0; j < in_nrow; j++)
				colbuf[j] = mat1[j * out_nrow];
			compute_dotprods2_with_double_Lcol(colbuf, in_nrow,
						SVT2, out, out_nrow, out_ncol);
			mat1++;
		}
	} else {
		for (i = 0; i < out_nrow; i++, out++) {
			compute_dotprods2_with_double_Lcol(mat1, in_nrow,
						SVT2, out, out_nrow, out_ncol);
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
			/* Copy 'mat2[j, ]' to 'colbuf'. */
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
			/* Copy 'mat1[i, ]' to 'colbuf'. */
			for (j = 0; j < in_nrow; j++)
				colbuf[j] = mat1[j * out_nrow];
			compute_dotprods2_with_int_Lcol(colbuf, in_nrow,
						SVT2, out, out_nrow, out_ncol);
			mat1++;
		}
	} else {
		for (i = 0; i < out_nrow; i++, out++) {
			compute_dotprods2_with_int_Lcol(mat1, in_nrow,
						SVT2, out, out_nrow, out_ncol);
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
   Note that 'SVT2' must be a leaf of type "double" because at the time
   we implemented this _dotprod_intSV_zero() was not available. */
static void crossprod2_mat0_SVT_double(SEXP SVT2, int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	if (SVT2 == R_NilValue)
		return;
	for (int j = 0; j < out_ncol; j++, out += out_nrow) {
		SEXP leaf = VECTOR_ELT(SVT2, j);
		if (leaf == R_NilValue)
			continue;
		const SparseVec sv = leaf2SV(leaf, in_nrow);
		fill_col(out, out_nrow, _dotprod_doubleSV_zero(&sv));
	}
	return;
}

/* Cross-product between 'SVT1' (on the left) and a fictive matrix
   of zeros (on the right).
   Note that 'SVT1' must be a leaf of type "double" because at the time
   we implemented this _dotprod_intSV_zero() was not available. */
static void crossprod2_SVT_mat0_double(SEXP SVT1, int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	if (SVT1 == R_NilValue)
		return;
	for (int i = 0; i < out_nrow; i++, out++) {
		SEXP leaf = VECTOR_ELT(SVT1, i);
		if (leaf == R_NilValue)
			continue;
		const SparseVec sv = leaf2SV(leaf, in_nrow);
		fill_row(out, out_nrow, out_ncol, _dotprod_doubleSV_zero(&sv));
	}
	return;
}

/* Cross-product between a fictive matrix of zeros (on the left)
   and 'SVT2' (on the right).
   'SVT2' must contain leaf vectors of type "integer". */
static void crossprod2_mat0_SVT_int(SEXP SVT2, int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	if (SVT2 == R_NilValue)
		return;
	for (int j = 0; j < out_ncol; j++, out += out_nrow) {
		SEXP leaf = VECTOR_ELT(SVT2, j);
		if (leaf == R_NilValue)
			continue;
		const SparseVec sv = leaf2SV(leaf, in_nrow);
		if (intSV_has_no_NA(&sv))
			continue;
		fill_col(out, out_nrow, NA_REAL);
	}
	return;
}

/* Cross-product between 'SVT1' (on the left) and a fictive matrix
   of zeros (on the right).
   'SVT1' must contain leaf vectors of type "integer". */
static void crossprod2_SVT_mat0_int(SEXP SVT1, int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	if (SVT1 == R_NilValue)
		return;
	for (int i = 0; i < out_nrow; i++, out++) {
		SEXP leaf = VECTOR_ELT(SVT1, i);
		if (leaf == R_NilValue)
			continue;
		const SparseVec sv = leaf2SV(leaf, in_nrow);
		if (intSV_has_no_NA(&sv))
			continue;
		fill_row(out, out_nrow, out_ncol, NA_REAL);
	}
	return;
}

/* Preprocesses 'leaf1' by turning it into a dense vector, but only if it
   contains no infinite values (i.e. no NA, NaN, Inf, or -Inf). */
static void compute_dotprods2_with_left_double_leaf(SEXP leaf1, SEXP SVT2,
		double *densebuf, int dense_len,
		double *out, int out_nrow, int out_ncol)
{
	if (leaf1 == R_NilValue) {
		memset(densebuf, 0, sizeof(double) * dense_len);
		compute_dotprods2_with_finite_Lcol(densebuf, dense_len, SVT2,
						   out, out_nrow, out_ncol);
		return;
	}
	const SparseVec sv1 = leaf2SV(leaf1, dense_len);
	if (doubleSV_has_no_NaN_or_Inf(&sv1)) {
		/* Turn 'sv1' into dense vector. */
		expand_doubleSV(&sv1, densebuf);
		compute_dotprods2_with_finite_Lcol(densebuf, dense_len, SVT2,
						   out, out_nrow, out_ncol);
		return;
	}
	compute_dotprods2_with_Lsv(&sv1, SVT2, out, out_nrow, out_ncol);
	return;
}

/* Preprocesses 'leaf2' by turning it into a dense vector, but only if it
   contains no infinite values (i.e. no NA, NaN, Inf, or -Inf). */
static void compute_dotprods2_with_right_double_leaf(SEXP SVT1, SEXP leaf2,
		double *densebuf, int dense_len,
		double *out, int out_nrow)
{
	if (leaf2 == R_NilValue) {
		memset(densebuf, 0, sizeof(double) * dense_len);
		compute_dotprods2_with_finite_Rcol(SVT1, densebuf, dense_len,
						   out, out_nrow);
		return;
	}
	const SparseVec sv2 = leaf2SV(leaf2, dense_len);
	if (doubleSV_has_no_NaN_or_Inf(&sv2)) {
		/* Turn 'sv2' into dense vector. */
		expand_doubleSV(&sv2, densebuf);
		compute_dotprods2_with_finite_Rcol(SVT1, densebuf, dense_len,
						   out, out_nrow);
		return;
	}
	compute_dotprods2_with_Rsv(SVT1, &sv2, out, out_nrow);
	return;
}

/* Preprocesses 'leaf1' by turning it into a dense vector, but only if it
   contains no NAs. */
static void compute_dotprods2_with_left_int_leaf(SEXP leaf1, SEXP SVT2,
		int *densebuf, int dense_len,
		double *out, int out_nrow, int out_ncol)
{
	if (leaf1 == R_NilValue) {
		memset(densebuf, 0, sizeof(int) * dense_len);
		compute_dotprods2_with_noNA_int_Lcol(densebuf, dense_len, SVT2,
						     out, out_nrow, out_ncol);
		return;
	}
	const SparseVec sv1 = leaf2SV(leaf1, dense_len);
	if (intSV_has_no_NA(&sv1)) {
		/* Turn 'sv1' into dense vector. */
		expand_intSV(&sv1, densebuf);
		compute_dotprods2_with_noNA_int_Lcol(densebuf, dense_len, SVT2,
						     out, out_nrow, out_ncol);
		return;
	}
	fill_row(out, out_nrow, out_ncol, NA_REAL);
	return;
}

/* Preprocesses 'leaf2' by turning it into a dense vector, but only if it
   contains no NAs. */
static void compute_dotprods2_with_right_int_leaf(SEXP SVT1, SEXP leaf2,
		int *densebuf, int dense_len,
		double *out, int out_nrow)
{
	if (leaf2 == R_NilValue) {
		memset(densebuf, 0, sizeof(int) * dense_len);
		compute_dotprods2_with_noNA_int_Rcol(SVT1, densebuf, dense_len,
						     out, out_nrow);
		return;
	}
	const SparseVec sv2 = leaf2SV(leaf2, dense_len);
	if (intSV_has_no_NA(&sv2)) {
		/* Turn 'sv2' into dense vector. */
		expand_intSV(&sv2, densebuf);
		compute_dotprods2_with_noNA_int_Rcol(SVT1, densebuf, dense_len,
						     out, out_nrow);
		return;
	}
	fill_col(out, out_nrow, NA_REAL);
	return;
}

/* Fills 'out' row by row. Each leaf in 'SVT1' will be preprocessed
   once (see compute_dotprods2_with_left_double_leaf() for details). */
static void crossprod2_Lpp_double(SEXP SVT1, SEXP SVT2, int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	if (SVT2 == R_NilValue) {
		crossprod2_SVT_mat0_double(SVT1, in_nrow,
					   out, out_nrow, out_ncol);
		return;
	}
	double *densebuf = (double *) R_alloc(in_nrow, sizeof(double));
	for (int i = 0; i < out_nrow; i++) {
		SEXP leaf = SVT1 != R_NilValue ? VECTOR_ELT(SVT1, i) :
						 R_NilValue;
		compute_dotprods2_with_left_double_leaf(leaf, SVT2,
						densebuf, in_nrow,
						out, out_nrow, out_ncol);
		out++;
	}
	return;
}

/* Fills 'out' column by column. Each leaf in 'SVT2' will be preprocessed
   once (see compute_dotprods2_with_right_double_leaf() for details). */
static void crossprod2_Rpp_double(SEXP SVT1, SEXP SVT2, int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	if (SVT1 == R_NilValue) {
		crossprod2_mat0_SVT_double(SVT2, in_nrow,
					   out, out_nrow, out_ncol);
		return;
	}
	double *densebuf = (double *) R_alloc(in_nrow, sizeof(double));
	for (int j = 0; j < out_ncol; j++) {
		SEXP leaf = SVT2 != R_NilValue ? VECTOR_ELT(SVT2, j) :
						 R_NilValue;
		compute_dotprods2_with_right_double_leaf(SVT1, leaf,
						densebuf, in_nrow,
						out, out_nrow);
		out += out_nrow;
	}
	return;
}

/* crossprod2_Lpp_int() and crossprod2_Rpp_int() below are the workhorses
   behind binary crossprod() when the input objects are SVT_SparseMatrix
   objects of type "integer". They work **natively** on the object. However,
   it turns out that going thru crossprod2_Lpp_double() or
   crossprod2_Rpp_double() by doing:

       crossprod(`type<-`(x, "double"), `type<-`(y, "double"))

   is about as fast, if not slightly faster. So only benefit of
   crossprod2_[L|R]pp_int() over crossprod2_[L|R]pp_double() is that they
   keep memory footprint as small as possible. */

static void crossprod2_Lpp_int(SEXP SVT1, SEXP SVT2, int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	if (SVT2 == R_NilValue) {
		crossprod2_SVT_mat0_int(SVT1, in_nrow,
					out, out_nrow, out_ncol);
		return;
	}
	int *densebuf = (int *) R_alloc(in_nrow, sizeof(int));
	for (int i = 0; i < out_nrow; i++) {
		SEXP leaf = SVT1 != R_NilValue ? VECTOR_ELT(SVT1, i) :
						 R_NilValue;
		compute_dotprods2_with_left_int_leaf(leaf, SVT2,
						densebuf, in_nrow,
                				out, out_nrow, out_ncol);
		out++;
	}
	return;
}

static void crossprod2_Rpp_int(SEXP SVT1, SEXP SVT2, int in_nrow,
		double *out, int out_nrow, int out_ncol)
{
	if (SVT1 == R_NilValue) {
		crossprod2_mat0_SVT_int(SVT2, in_nrow,
					out, out_nrow, out_ncol);
		return;
	}
	int *densebuf = (int *) R_alloc(in_nrow, sizeof(int));
	for (int j = 0; j < out_ncol; j++) {
		SEXP leaf = SVT2 != R_NilValue ? VECTOR_ELT(SVT2, j) :
						 R_NilValue;
		compute_dotprods2_with_right_int_leaf(SVT1, leaf,
						densebuf, in_nrow,
						out, out_nrow);
		out += out_nrow;
	}
	return;
}


/****************************************************************************
 * Workhorses behind C_crossprod1_SVT()
 */

static void compute_sym_dotprods_double(SEXP SVT, int j,
		double *densebuf, int dense_len, double *out, int out_ncol)
{
	SEXP leaf = VECTOR_ELT(SVT, j);
	if (leaf == R_NilValue) {
		memset(densebuf, 0, sizeof(double) * dense_len);
		compute_sym_dotprods_with_finite_col(SVT, j,
						densebuf, dense_len,
						out, out_ncol);
		return;
	}
	const SparseVec sv = leaf2SV(leaf, dense_len);
	if (doubleSV_has_no_NaN_or_Inf(&sv)) {
		/* Turn 'sv' into dense vector. */
		expand_doubleSV(&sv, densebuf);
		*out = _dotprod_doubleSV_finite_doubles(&sv, densebuf);
		compute_sym_dotprods_with_finite_col(SVT, j,
						densebuf, dense_len,
						out, out_ncol);
	} else {
		*out = _dotprod_doubleSV_doubleSV(&sv, &sv);
		compute_sym_dotprods_with_doubleSV(SVT, j, &sv, out, out_ncol);
	}
	return;
}

static void compute_sym_dotprods_int(SEXP SVT, int j,
		int *densebuf, int dense_len, double *out, int out_ncol)
{
	SEXP leaf = VECTOR_ELT(SVT, j);
	if (leaf == R_NilValue) {
		memset(densebuf, 0, sizeof(int) * dense_len);
		compute_sym_dotprods_with_noNA_int_col(SVT, j,
						densebuf, dense_len,
						out, out_ncol);
		return;
	}
	const SparseVec sv = leaf2SV(leaf, dense_len);
	if (intSV_has_no_NA(&sv)) {
		/* Turn 'sv' into dense vector. */
		expand_intSV(&sv, densebuf);
		*out = _dotprod_intSV_noNA_ints(&sv, densebuf);
		compute_sym_dotprods_with_noNA_int_col(SVT, j,
						densebuf, dense_len,
						out, out_ncol);
	} else {
		sym_fill_with_NAs(out, out_ncol, j);
	}
	return;
}

static void crossprod1_double(SEXP SVT, int in_nrow, double *out, int out_ncol)
{
	if (SVT == R_NilValue)
		return;
	double *densebuf = (double *) R_alloc(in_nrow, sizeof(double));
	for (int j = 0; j < out_ncol; j++, out += out_ncol + 1)
		compute_sym_dotprods_double(SVT, j,
					    densebuf, in_nrow, out, out_ncol);
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
	if (SVT == R_NilValue)
		return;
	int *densebuf = (int *) R_alloc(in_nrow, sizeof(int));
	for (int j = 0; j < out_ncol; j++, out += out_ncol + 1)
		compute_sym_dotprods_int(SVT, j,
					 densebuf, in_nrow, out, out_ncol);
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
	R_xlen_t Lpp_nops, Rpp_nops;

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
	Lpp_nops = _REC_nzcount_SVT(y_SVT, LENGTH(y_dim)) * x_ncol;
	Rpp_nops = _REC_nzcount_SVT(x_SVT, LENGTH(x_dim)) * y_ncol;
	if (Lpp_nops < Rpp_nops) {
		/* Will preprocess cols in 'x_SVT' (left preprocessing). */
		if (x_Rtype == REALSXP) {
			crossprod2_Lpp_double(x_SVT, y_SVT, in_nrow,
					REAL(ans), x_ncol, y_ncol);
		} else {
			crossprod2_Lpp_int(x_SVT, y_SVT, in_nrow,
					REAL(ans), x_ncol, y_ncol);
		}
	} else {
		/* Will preprocess cols in 'y_SVT' (right preprocessing). */
		if (x_Rtype == REALSXP) {
			crossprod2_Rpp_double(x_SVT, y_SVT, in_nrow,
					REAL(ans), x_ncol, y_ncol);
		} else {
			crossprod2_Rpp_int(x_SVT, y_SVT, in_nrow,
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

