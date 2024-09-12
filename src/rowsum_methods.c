/****************************************************************************
 *     rowsum()/colsum() methods for SparseArray and dgCMatrix objects      *
 ****************************************************************************/
#include "rowsum_methods.h"

#include "S4Vectors_interface.h"

#include "argcheck_utils.h"
#include "leaf_utils.h"

#include <limits.h>  /* for INT_MAX */


/* Copied from S4Arrays/src/rowsum.c */
static void check_group(SEXP group, int x_nrow, int ngroup)
{
	if (!IS_INTEGER(group))
		error("the grouping vector must be "
		      "an integer vector or factor");
	if (LENGTH(group) != x_nrow)
		error("the grouping vector must have one element "
		      "per row in 'x' for rowsum()\n  and one element "
		      "per column in 'x' for colsum()");
	for (int i = 0; i < x_nrow; i++) {
		int g = INTEGER(group)[i];
		if (g == NA_INTEGER) {
			if (ngroup < 1)
				error("'ngroup' must be >= 1 when 'group' "
				      "contains missing values");
		} else {
			if (g < 1 || g > ngroup)
				error("all non-NA values in 'group' must "
				      "be >= 1 and <= 'ngroup'");
		}
	}
	return;
}


/****************************************************************************
 * Low-level helpers used by C_rowsum_SVT() and C_rowsum_dgCMatrix()
 */

static void compute_rowsum_doubles(
		const double *nzvals, const int *nzoffs, int nzcount,
		const int *groups, int narm, double *out, int out_len)
{
	for (int k = 0; k < nzcount; k++) {
		int g = groups[nzoffs[k]];
		if (g == NA_INTEGER)
			g = out_len;
		g--;  // from 1-base to 0-base
		double v = double1;
		if (nzvals != NULL) {
			v = nzvals[k];
			/* ISNAN(): True for *both* NA and NaN.
			   See <R_ext/Arith.h> */
			if (narm && ISNAN(v))
				continue;
		}
		out[g] += v;
	}
	return;
}

static void compute_rowsum_ints(
		const int *nzvals, const int *nzoffs, int nzcount,
		const int *groups, int narm, int *out, int out_len)
{
	for (int k = 0; k < nzcount; k++) {
		int g = groups[nzoffs[k]];
		if (g == NA_INTEGER)
			g = out_len;
		g--;  // from 1-base to 0-base
		int v = int1;
		if (nzvals != NULL) {
			v = nzvals[k];
			if (narm && v == NA_INTEGER)
				continue;
		}
		out[g] = safe_int_add(out[g], v);
	}
	return;
}

static void rowsum_SVT_double(SEXP x_SVT, int x_nrow, int x_ncol,
		const int *groups, int narm, double *out, int out_nrow)
{
	if (x_SVT == R_NilValue)
		return;
	for (int j = 0; j < x_ncol; j++, out += out_nrow) {
		SEXP subSVT = VECTOR_ELT(x_SVT, j);
		if (subSVT == R_NilValue)
			continue;
		SEXP nzvals, nzoffs;
		int nzcount = unzip_leaf(subSVT, &nzvals, &nzoffs);
		compute_rowsum_doubles(
			nzvals == R_NilValue ? NULL : REAL(nzvals),
			INTEGER(nzoffs), nzcount, groups, narm,
			out, out_nrow);
	}
	return;
}

static void rowsum_SVT_int(SEXP x_SVT, int x_nrow, int x_ncol,
		const int *groups, int narm, int *out, int out_nrow)
{
	if (x_SVT == R_NilValue)
		return;
	reset_ovflow_flag();
	for (int j = 0; j < x_ncol; j++, out += out_nrow) {
		SEXP subSVT = VECTOR_ELT(x_SVT, j);
		if (subSVT == R_NilValue)
			continue;
		SEXP nzvals, nzoffs;
		int nzcount = unzip_leaf(subSVT, &nzvals, &nzoffs);
		compute_rowsum_ints(
			nzvals == R_NilValue ? NULL : INTEGER(nzvals),
			INTEGER(nzoffs), nzcount, groups, narm,
			out, out_nrow);
	}
	if (get_ovflow_flag())
		warning("NAs produced by integer overflow");
	return;
}

static void rowsum_dgCMatrix(int x_nrow, int x_ncol,
		const double *x_slotx, const int *x_sloti, const int *x_slotp,
		const int *groups, int narm, double *out, int out_nrow)
{
	for (int j = 0; j < x_ncol; j++, out += out_nrow) {
		int offset = x_slotp[j];
		int nzcount = x_slotp[j + 1] - offset;
		compute_rowsum_doubles(
			x_slotx + offset, x_sloti + offset, nzcount,
			groups, narm, out, out_nrow);
	}
	return;
}


/****************************************************************************
 * Low-level helpers used by C_colsum_SVT() and C_colsum_dgCMatrix()
 */

/* Add sparse "double" vector represented by ('nzvals', 'nzoffs', 'nzcount')
   to dense "double" vector 'out'. */
static void add_sparse_vec_to_doubles(
		const double *nzvals, const int *nzoffs, int nzcount,
		double *out, int narm)
{
	for (int k = 0; k < nzcount; k++) {
		double x;
		if (nzvals == NULL) {
			/* lacunar leaf */
			x = double1;
		} else {
			/* regular leaf */
			x = nzvals[k];
			/* ISNAN(): True for *both* NA and NaN.
			   See <R_ext/Arith.h> */
			if (narm && ISNAN(x))
				continue;
		}
		out[nzoffs[k]] += x;
	}
	return;
}

/* Add sparse "int" vector represented by ('nzvals', 'nzoffs', 'nzcount')
   to dense "int" vector 'out'. */
static void add_sparse_vec_to_ints(
		const int *nzvals, const int *nzoffs, int nzcount,
		int *out, int narm, int *overflow)
{
	for (int k = 0; k < nzcount; k++) {
		int *out_p = out + nzoffs[k];
		if (*out_p == NA_INTEGER)
			continue;
		int x;
		if (nzvals == NULL) {
			/* lacunar leaf */
			x = int1;
		} else  {
			/* regular leaf */
			x = nzvals[k];
			if (x == NA_INTEGER) {
				if (!narm)
					*out_p = NA_INTEGER;
				continue;
			}
		}
		double y = (double) *out_p + x;
		if (-INT_MAX <= y && y <= INT_MAX) {
			*out_p = (int) y;
		} else {
			*overflow = 1;
			*out_p = NA_INTEGER;
		}
	}
	return;
}

static void colsum_SVT_double(SEXP x_SVT, int x_nrow, int x_ncol,
		const int *groups, int narm, double *out, int out_ncol)
{
	if (x_SVT == R_NilValue)
		return;
	for (int j = 0; j < x_ncol; j++) {
		SEXP subSVT = VECTOR_ELT(x_SVT, j);
		if (subSVT == R_NilValue)
			continue;
		SEXP nzvals, nzoffs;
		int nzcount = unzip_leaf(subSVT, &nzvals, &nzoffs);
		const double *nzvals_p = NULL;
		if (nzvals != R_NilValue)
			nzvals_p = REAL(nzvals);
		int g = groups[j];
		if (g == NA_INTEGER)
			g = out_ncol;
		g--;  // from 1-base to 0-base
		add_sparse_vec_to_doubles(
				nzvals_p, INTEGER(nzoffs), nzcount,
				out + g * x_nrow, narm);
	}
	return;
}

static void colsum_SVT_int(SEXP x_SVT, int x_nrow, int x_ncol,
		const int *groups, int narm, int *out, int out_ncol)
{
	if (x_SVT == R_NilValue)
		return;
	int overflow = 0;
	for (int j = 0; j < x_ncol; j++) {
		SEXP subSVT = VECTOR_ELT(x_SVT, j);
		if (subSVT == R_NilValue)
			continue;
		SEXP nzvals, nzoffs;
		int nzcount = unzip_leaf(subSVT, &nzvals, &nzoffs);
		const int *nzvals_p = NULL;
		if (nzvals != R_NilValue)
			nzvals_p = INTEGER(nzvals);
		int g = groups[j];
		if (g == NA_INTEGER)
			g = out_ncol;
		g--;  // from 1-base to 0-base
		add_sparse_vec_to_ints(
				nzvals_p, INTEGER(nzoffs), nzcount,
				out + g * x_nrow, narm, &overflow);
	}
	if (overflow)
		warning("NAs produced by integer overflow");
	return;
}

static void colsum_dgCMatrix(int x_nrow, int x_ncol,
		const double *x_slotx, const int *x_sloti, const int *x_slotp,
		const int *groups, int narm, double *out, int out_ncol)
{
	for (int j = 0; j < x_ncol; j++) {
		int offset = x_slotp[j];
		int nzcount = x_slotp[j + 1] - offset;
		int g = groups[j];
		if (g == NA_INTEGER)
			g = out_ncol;
		g--;  // from 1-base to 0-base
		add_sparse_vec_to_doubles(
				x_slotx + offset, x_sloti + offset, nzcount,
				out + g * x_nrow, narm);
	}
	return;
}


/****************************************************************************
 * C_rowsum_SVT() and C_rowsum_dgCMatrix()
 */

/* --- .Call ENTRY POINT --- */
SEXP C_rowsum_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		  SEXP group, SEXP ngroup, SEXP na_rm)
{
	if (LENGTH(x_dim) != 2)
		error("input object must have 2 dimensions");
	int x_nrow = INTEGER(x_dim)[0];
	int x_ncol = INTEGER(x_dim)[1];
	int narm = LOGICAL(na_rm)[0];

	SEXPTYPE x_Rtype = _get_and_check_Rtype_from_Rstring(x_type,
						"C_rowsum_SVT", "x_type");

	int ans_nrow = INTEGER(ngroup)[0];
	check_group(group, x_nrow, ans_nrow);

	reset_ovflow_flag();
	/* Only to detect a potential integer overflow. The returned value
	   is actually not needed so we ignore it. */
	safe_int_mult(ans_nrow, x_ncol);
	if (get_ovflow_flag())
		error("too many groups (matrix of sums will be too big)");

	/* Note that base::rowsum() only supports numeric matrices i.e.
	   matrices of type() "double" or "integer", so we do the same. */
	SEXP ans;
	if (x_Rtype == REALSXP) {
		ans = PROTECT(_new_Rmatrix0(REALSXP, ans_nrow, x_ncol,
					    R_NilValue));
		rowsum_SVT_double(x_SVT, x_nrow, x_ncol,
			INTEGER(group), narm, REAL(ans), ans_nrow);
	} else if (x_Rtype == INTSXP) {
		ans = PROTECT(_new_Rmatrix0(INTSXP, ans_nrow, x_ncol,
					    R_NilValue));
		rowsum_SVT_int(x_SVT, x_nrow, x_ncol,
			INTEGER(group), narm, INTEGER(ans), ans_nrow);
	} else {
		error("rowsum() and colsum() do not support "
		      "SVT_SparseMatrix objects of\n"
		      "  type \"%s\" at the moment",
		      type2char(x_Rtype));
	}

	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_rowsum_dgCMatrix(SEXP x, SEXP group, SEXP ngroup, SEXP na_rm)
{
	SEXP x_Dim = GET_SLOT(x, install("Dim"));
	int x_nrow = INTEGER(x_Dim)[0];
	int x_ncol = INTEGER(x_Dim)[1];
	SEXP x_slotx = GET_SLOT(x, install("x"));
	SEXP x_sloti = GET_SLOT(x, install("i"));
	SEXP x_slotp = GET_SLOT(x, install("p"));
	int narm = LOGICAL(na_rm)[0];

	int ans_nrow = INTEGER(ngroup)[0];
	check_group(group, x_nrow, ans_nrow);

	reset_ovflow_flag();
	/* Only to detect a potential integer overflow. The returned value
	   is actually not needed so we ignore it. */
	safe_int_mult(ans_nrow, x_ncol);
	if (get_ovflow_flag())
		error("too many groups (matrix of sums will be too big)");
	SEXP ans =
		PROTECT(_new_Rmatrix0(REALSXP, ans_nrow, x_ncol, R_NilValue));

	rowsum_dgCMatrix(x_nrow, x_ncol,
			 REAL(x_slotx), INTEGER(x_sloti), INTEGER(x_slotp),
			 INTEGER(group), narm, REAL(ans), ans_nrow);

	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_colsum_SVT() and C_colsum_dgCMatrix()
 */

/* --- .Call ENTRY POINT --- */
SEXP C_colsum_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		  SEXP group, SEXP ngroup, SEXP na_rm)
{
	if (LENGTH(x_dim) != 2)
		error("input object must have 2 dimensions");
	int x_nrow = INTEGER(x_dim)[0];
	int x_ncol = INTEGER(x_dim)[1];
	int narm = LOGICAL(na_rm)[0];

	SEXPTYPE x_Rtype = _get_and_check_Rtype_from_Rstring(x_type,
						"C_colsum_SVT", "x_type");

	int ans_ncol = INTEGER(ngroup)[0];
	check_group(group, x_ncol, ans_ncol);

	reset_ovflow_flag();
	/* Only to detect a potential integer overflow. The returned value
	   is actually not needed so we ignore it. */
	safe_int_mult(x_nrow, ans_ncol);
	if (get_ovflow_flag())
		error("too many groups (matrix of sums will be too big)");

	/* Note that base::rowsum() only supports numeric matrices i.e.
	   matrices of type() "double" or "integer", so we do the same. */
	SEXP ans;
	if (x_Rtype == REALSXP) {
		ans = PROTECT(_new_Rmatrix0(REALSXP, x_nrow, ans_ncol,
					    R_NilValue));
		colsum_SVT_double(x_SVT, x_nrow, x_ncol,
			INTEGER(group), narm, REAL(ans), ans_ncol);
	} else if (x_Rtype == INTSXP) {
		ans = PROTECT(_new_Rmatrix0(INTSXP, x_nrow, ans_ncol,
					    R_NilValue));
		colsum_SVT_int(x_SVT, x_nrow, x_ncol,
			INTEGER(group), narm, INTEGER(ans), ans_ncol);
	} else {
		error("rowsum() and colsum() do not support "
		      "SVT_SparseMatrix objects of\n"
		      "  type \"%s\" at the moment",
		      type2char(x_Rtype));
	}

	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_colsum_dgCMatrix(SEXP x, SEXP group, SEXP ngroup, SEXP na_rm)
{
	SEXP x_Dim = GET_SLOT(x, install("Dim"));
	int x_nrow = INTEGER(x_Dim)[0];
	int x_ncol = INTEGER(x_Dim)[1];
	SEXP x_slotx = GET_SLOT(x, install("x"));
	SEXP x_sloti = GET_SLOT(x, install("i"));
	SEXP x_slotp = GET_SLOT(x, install("p"));
	int narm = LOGICAL(na_rm)[0];

	int ans_ncol = INTEGER(ngroup)[0];
	check_group(group, x_ncol, ans_ncol);

	reset_ovflow_flag();
	/* Only to detect a potential integer overflow. The returned value
	   is actually not needed so we ignore it. */
	safe_int_mult(x_nrow, ans_ncol);
	if (get_ovflow_flag())
		error("too many groups (matrix of sums will be too big)");
	SEXP ans =
		PROTECT(_new_Rmatrix0(REALSXP, x_nrow, ans_ncol, R_NilValue));

	colsum_dgCMatrix(x_nrow, x_ncol,
			 REAL(x_slotx), INTEGER(x_sloti), INTEGER(x_slotp),
			 INTEGER(group), narm, REAL(ans), ans_ncol);

	UNPROTECT(1);
	return ans;
}

