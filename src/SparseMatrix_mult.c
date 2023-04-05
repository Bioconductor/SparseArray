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

static int lv_has_no_NA(SEXP lv)
{
	int lv_len, k;
	SEXP lv_offs, lv_vals;
	const int *vals_p;

	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	vals_p = INTEGER(lv_vals);
	for (k = 0; k < lv_len; k++) {
		if (*vals_p == NA_INTEGER)
			return 0;
		vals_p++;
	}
	return 1;
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
 * Workhorses behind unary crossprod()
 */

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
						out, out_ncol, j);
		} else if (lv_is_finite(subSVT)) {
			expand_double_lv(subSVT, col, in_nrow);
			*out =
			    _dotprod_leaf_vector_and_finite_col(subSVT, col);
			compute_sym_dotprods_with_finite_col(SVT, col,
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
   input object is a SVT_SparseMatrix object of type "integer". It
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
	int *col;

	if (SVT == R_NilValue)
		return;
	col = (int *) R_alloc(in_nrow, sizeof(int));
	for (j = 0; j < out_ncol; j++, out += out_ncol + 1) {
		subSVT = VECTOR_ELT(SVT, j);
		if (subSVT == R_NilValue) {
			memset(col, 0, sizeof(int) * in_nrow);
			compute_sym_dotprods_with_noNA_int_col(SVT, col,
						out, out_ncol, j);
		} else if (lv_has_no_NA(subSVT)) {
			expand_int_lv(subSVT, col, in_nrow);
			*out =
			    _dotprod_leaf_vector_and_noNA_int_col(subSVT, col);
			compute_sym_dotprods_with_noNA_int_col(SVT, col,
						out, out_ncol, j);
		} else {
			sym_fill_with_NAs(out, out_ncol, j);
		}
	}
	return;
}


/****************************************************************************
 * Workhorses behind binary crossprod()
 */

/* Cross-product between a fictive matrix of zeros (on the left)
   and 'SVT2' (on the right).
   Note that 'SVT2' must contain leaf vectors of type "double", because
   that's the only type supported by _dotprod0_leaf_vector() for now. */
static void crossprod2_Lzeros_double(SEXP SVT2,
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
static void crossprod2_Rzeros_double(SEXP SVT1,
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
   'SVT2' must contain leaf vectors of type "int". */
static void crossprod2_Lzeros_int(SEXP SVT2,
		double *out, int out_nrow, int out_ncol)
{
	int j;
	SEXP subSVT2;

	if (SVT2 == R_NilValue)
		return;
	for (j = 0; j < out_ncol; j++, out += out_nrow) {
		subSVT2 = VECTOR_ELT(SVT2, j);
		if (subSVT2 == R_NilValue || lv_has_no_NA(subSVT2))
			continue;
		fill_col_with_NAs(out, out_nrow);
	}
	return;
}

/* Cross-product between 'SVT1' (on the left) and a fictive matrix
   of zeros (on the right).
   'SVT1' must contain leaf vectors of type "int". */
static void crossprod2_Rzeros_int(SEXP SVT1,
		double *out, int out_nrow, int out_ncol)
{
	int i;
	SEXP subSVT1;

	if (SVT1 == R_NilValue)
		return;
	for (i = 0; i < out_nrow; i++, out++) {
		subSVT1 = VECTOR_ELT(SVT1, i);
		if (subSVT1 == R_NilValue || lv_has_no_NA(subSVT1))
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
	double *col;

	if (SVT1 == R_NilValue) {
		crossprod2_Lzeros_double(SVT2, out, out_nrow, out_ncol);
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
		crossprod2_Rzeros_double(SVT1, out, out_nrow, out_ncol);
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
	int *col;

	if (SVT1 == R_NilValue) {
		crossprod2_Lzeros_int(SVT2, out, out_nrow, out_ncol);
		return;
	}
	col = (int *) R_alloc(in_nrow, sizeof(int));
	subSVT2 = R_NilValue;
	for (j = 0; j < out_ncol; j++, out += out_nrow) {
		if (SVT2 != R_NilValue)
			subSVT2 = VECTOR_ELT(SVT2, j);
		if (subSVT2 == R_NilValue) {
			memset(col, 0, sizeof(int) * in_nrow);
			compute_dotprods2_with_noNA_int_Rcol(SVT1, col,
						out, out_nrow);
		} else if (lv_has_no_NA(subSVT2)) {
			/* Preprocess 'subSVT2'. */
			expand_int_lv(subSVT2, col, in_nrow);
			compute_dotprods2_with_noNA_int_Rcol(SVT1, col,
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
	int *col;

	if (SVT2 == R_NilValue) {
		crossprod2_Rzeros_int(SVT1, out, out_nrow, out_ncol);
		return;
	}
	col = (int *) R_alloc(in_nrow, sizeof(int));
	subSVT1 = R_NilValue;
	for (i = 0; i < out_nrow; i++, out++) {
		if (SVT1 != R_NilValue)
			subSVT1 = VECTOR_ELT(SVT1, i);
		if (subSVT1 == R_NilValue) {
			memset(col, 0, sizeof(int) * in_nrow);
			compute_dotprods2_with_noNA_int_Lcol(col, SVT2,
						out, out_nrow, out_ncol);
		} else if (lv_has_no_NA(subSVT1)) {
			/* Preprocess 'subSVT1'. */
			expand_int_lv(subSVT1, col, in_nrow);
			compute_dotprods2_with_noNA_int_Lcol(col, SVT2,
						out, out_nrow, out_ncol);
		} else {
			fill_row_with_NAs(out, out_nrow, out_ncol);
		}
	}
	return;
}


/****************************************************************************
 * C_SVT_crossprod1() and C_SVT_crossprod2()
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
SEXP C_SVT_crossprod1(SEXP x_dim, SEXP x_type, SEXP x_SVT,
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
		      "C_SVT_crossprod1():\n"
		      "    invalid 'ans_type' value");
	if (ans_Rtype != REALSXP)
		error("SparseArray internal error in "
		      "C_SVT_crossprod1():\n"
		      "    output type \"%s\" is not supported yet",
		      type2char(ans_Rtype));

	/* Allocate 'ans' and fill with zeros. */
	ans = PROTECT(new_Rmatrix(ans_Rtype, x_ncol, x_ncol, ans_dimnames));

	if (x_Rtype == REALSXP) {
		crossprod1_double(x_SVT, x_nrow, REAL(ans), x_ncol);
	} else {
		crossprod1_int(x_SVT, x_nrow, REAL(ans), x_ncol);
	}

	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_SVT_crossprod2(SEXP x_dim, SEXP x_type, SEXP x_SVT,
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
		      "C_SVT_crossprod2():\n"
		      "    invalid 'ans_type' value");
	if (ans_Rtype != REALSXP)
		error("SparseArray internal error in "
		      "C_SVT_crossprod2():\n"
		      "    output type \"%s\" is not supported yet",
		      type2char(ans_Rtype));

	/* Allocate 'ans' and fill with zeros. */
	ans = PROTECT(new_Rmatrix(ans_Rtype, x_ncol, y_ncol, ans_dimnames));

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

