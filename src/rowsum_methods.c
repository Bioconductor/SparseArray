/****************************************************************************
 *          rowsum() methods for SparseArray and dgCMatrix objects          *
 ****************************************************************************/
#include "rowsum_methods.h"

#include "S4Vectors_interface.h"

#include "Rvector_utils.h"
#include "leaf_vector_utils.h"

static void check_group(SEXP group, int x_nrow, int ngroup)
{
	int i, g;

	if (!IS_INTEGER(group))
		error("the grouping vector must be "
		      "an integer vector or factor");
	if (LENGTH(group) != x_nrow)
		error("the grouping vector must have "
		      "one element per row in 'x'");
	for (i = 0; i < x_nrow; i++) {
		g = INTEGER(group)[i];
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

static void compute_rowsum_doubles(const double *vals, const int *offs, int n,
		const int *groups, double *out, int out_len, int narm)
{
	int k, g;
	double v;

	for (k = 0; k < n; k++) {
		g = groups[offs[k]];
		if (g == NA_INTEGER)
			g = out_len;
		g--;  // from 1-base to 0-base
		v = vals[k];
		/* ISNAN(): True for *both* NA and NaN. See <R_ext/Arith.h> */
		if (narm && ISNAN(v))
			continue;
		out[g] += v;
	}
	return;
}

static void compute_rowsum_ints(const int *vals, const int *offs, int n,
		const int *groups, double *out, int out_len, int narm)
{
	int k, g, v;

	for (k = 0; k < n; k++) {
		g = groups[offs[k]];
		if (g == NA_INTEGER)
			g = out_len;
		g--;  // from 1-base to 0-base
		v = vals[k];
		if (narm && v == NA_INTEGER)
			continue;
		out[g] += v;
	}
	return;
}

static int rowsum_SVT(int x_nrow, int x_ncol, SEXP x_SVT,
		const int *groups, int ngroup, int narm, double *out)
{
	int j, lv_len;
	SEXP subSVT, lv_offs, lv_vals;

	for (j = 0; j < x_ncol; j++) {
		subSVT = VECTOR_ELT(x_SVT, j);
		if (subSVT == R_NilValue)
			continue;
		lv_len = _split_leaf_vector(subSVT, &lv_offs, &lv_vals);
		if (TYPEOF(lv_vals) == REALSXP) {
			compute_rowsum_doubles(
				REAL(lv_vals), INTEGER(lv_offs), lv_len,
				groups, out, ngroup, narm);
		} else if (TYPEOF(lv_vals) == INTSXP) {
			compute_rowsum_ints(
				INTEGER(lv_vals), INTEGER(lv_offs), lv_len,
				groups, out, ngroup, narm);
		} else {
			return -1;
		}
		out += ngroup;
	}
	return 0;
}

static void rowsum_dgCMatrix(int x_nrow, int x_ncol,
		const double *x_x, const int *x_i, const int *x_p,
		const int *groups, int ngroup, int narm, double *out)
{
	int j, offset, nzcount;

	for (j = 0; j < x_ncol; j++) {
		offset = x_p[j];
		nzcount = x_p[j + 1] - offset;
		compute_rowsum_doubles(
			x_x + offset, x_i + offset, nzcount,
			groups, out, ngroup, narm);
		out += ngroup;
	}
	return;
}

/* --- .Call ENTRY POINT --- */
SEXP C_rowsum_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		  SEXP group, SEXP ngroup, SEXP na_rm)
{
	int x_nrow, x_ncol, narm, ans_nrow, ret;
	SEXP ans;

	if (LENGTH(x_dim) != 2)
		error("input object must have 2 dimensions");
	x_nrow = INTEGER(x_dim)[0];
	x_ncol = INTEGER(x_dim)[1];
	narm = LOGICAL(na_rm)[0];

	ans_nrow = INTEGER(ngroup)[0];
	check_group(group, x_nrow, ans_nrow);

	reset_ovflow_flag();
	/* Only to detect a potential integer overflow. The returned value
	   is actually not needed so we ignore it. */
	safe_int_mult(ans_nrow, x_ncol);
	if (get_ovflow_flag())
		error("too many groups (matrix of sums will be too big)");
	ans = PROTECT(_new_Rmatrix0(REALSXP, ans_nrow, x_ncol, R_NilValue));

	ret = rowsum_SVT(x_nrow, x_ncol, x_SVT,
			 INTEGER(group), ans_nrow, narm, REAL(ans));
	if (ret < 0) {
		UNPROTECT(1);
		error("SparseArray internal error in "
		      "C_rowsum_SVT():\n"
		      "    rowsum_SVT() returned an error");
	}

	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_rowsum_dgCMatrix(SEXP x, SEXP group, SEXP ngroup, SEXP na_rm)
{
	SEXP x_Dim, x_x, x_i, x_p, ans;
	int x_nrow, x_ncol, narm, ans_nrow;

	x_Dim = GET_SLOT(x, install("Dim"));
	x_nrow = INTEGER(x_Dim)[0];
	x_ncol = INTEGER(x_Dim)[1];
	x_x = GET_SLOT(x, install("x"));
	x_i = GET_SLOT(x, install("i"));
	x_p = GET_SLOT(x, install("p"));
	narm = LOGICAL(na_rm)[0];

	ans_nrow = INTEGER(ngroup)[0];
	check_group(group, x_nrow, ans_nrow);

	reset_ovflow_flag();
	/* Only to detect a potential integer overflow. The returned value
	   is actually not needed so we ignore it. */
	safe_int_mult(ans_nrow, x_ncol);
	if (get_ovflow_flag())
		error("too many groups (matrix of sums will be too big)");
	ans = PROTECT(_new_Rmatrix0(REALSXP, ans_nrow, x_ncol, R_NilValue));

	rowsum_dgCMatrix(x_nrow, x_ncol,
			 REAL(x_x), INTEGER(x_i), INTEGER(x_p),
			 INTEGER(group), ans_nrow, narm, REAL(ans));

	UNPROTECT(1);
	return ans;
}

