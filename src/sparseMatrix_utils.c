/****************************************************************************
 *             Some summarization methods for dgCMatrix objects             *
 ****************************************************************************/
#include "sparseMatrix_utils.h"


/****************************************************************************
 * colMins(), colMaxs(), colRanges()
 */

typedef double (*ExtremumFUN)(const double *x, int x_len, int narm,
			      int start_on_zero);

static double min_double(const double *x, int x_len, int narm,
			 int start_on_zero)
{
	double min = start_on_zero ? 0.0 : R_PosInf;
	int min_is_NaN = 0;
	for (int i = 0; i < x_len; i++) {
		double xi = x[i];
		if (R_IsNA(xi)) {
			if (narm)
				continue;
			return NA_REAL;
		}
		if (min_is_NaN)
			continue;
		if (R_IsNaN(xi)) {
			if (narm)
				continue;
			min = xi;
			min_is_NaN = 1;
			continue;
		}
		if (xi < min)
			min = xi;
	}
	return min;
}

static double max_double(const double *x, int x_len, int narm,
			 int start_on_zero)
{
	double max = start_on_zero ? 0.0 : R_NegInf;
	int max_slotis_NaN = 0;
	for (int i = 0; i < x_len; i++) {
		double xi = x[i];
		if (R_IsNA(xi)) {
			if (narm)
				continue;
			return NA_REAL;
		}
		if (max_slotis_NaN)
			continue;
		if (R_IsNaN(xi)) {
			if (narm)
				continue;
			max = xi;
			max_slotis_NaN = 1;
			continue;
		}
		if (xi > max)
			max = xi;
	}
	return max;
}

static void minmax_double(const double *x, int x_len, int narm,
			  int start_on_zero, double *min, double *max)
{
	double tmp_min, tmp_max;
	if (start_on_zero) {
		tmp_min = tmp_max = 0.0;
	} else {
		tmp_min = R_PosInf;
		tmp_max = R_NegInf;
	}
	int minmax_slotis_NaN = 0;
	for (int i = 0; i < x_len; i++) {
		double xi = x[i];
		if (R_IsNA(xi)) {
			if (narm)
				continue;
			*min = *max = NA_REAL;
			return;
		}
		if (minmax_slotis_NaN)
			continue;
		if (R_IsNaN(xi)) {
			if (narm)
				continue;
			tmp_min = tmp_max = xi;
			minmax_slotis_NaN = 1;
			continue;
		}
		if (xi < tmp_min)
			tmp_min = xi;
		if (xi > tmp_max)
			tmp_max = xi;
	}
	*min = tmp_min;
	*max = tmp_max;
	return;
}

static SEXP C_colExtrema_dgCMatrix(ExtremumFUN extremum_FUN,
		SEXP x, SEXP na_rm)
{
	SEXP x_Dim = GET_SLOT(x, install("Dim"));
	int x_nrow = INTEGER(x_Dim)[0];
	int x_ncol = INTEGER(x_Dim)[1];
	SEXP x_slotx = GET_SLOT(x, install("x"));
	SEXP x_slotp = GET_SLOT(x, install("p"));
	int narm = LOGICAL(na_rm)[0];

	SEXP ans = PROTECT(NEW_NUMERIC(x_ncol));
	for (int j = 0; j < x_ncol; j++) {
		int offset = INTEGER(x_slotp)[j];
		int nzcount = INTEGER(x_slotp)[j + 1] - offset;
		REAL(ans)[j] = extremum_FUN(REAL(x_slotx) + offset, nzcount,
					    narm, nzcount < x_nrow);
	}
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_colMins_dgCMatrix(SEXP x, SEXP na_rm)
{
	return C_colExtrema_dgCMatrix(min_double, x, na_rm);
}

/* --- .Call ENTRY POINT --- */
SEXP C_colMaxs_dgCMatrix(SEXP x, SEXP na_rm)
{
	return C_colExtrema_dgCMatrix(max_double, x, na_rm);
}

/* --- .Call ENTRY POINT ---
   About 2x faster than the colRanges() method for dgCMatrix objects defined
   in the sparseMatrixStats package. */
SEXP C_colRanges_dgCMatrix(SEXP x, SEXP na_rm)
{
	SEXP x_Dim = GET_SLOT(x, install("Dim"));
	int x_nrow = INTEGER(x_Dim)[0];
	int x_ncol = INTEGER(x_Dim)[1];
	SEXP x_slotx = GET_SLOT(x, install("x"));
	SEXP x_slotp = GET_SLOT(x, install("p"));
	int narm = LOGICAL(na_rm)[0];

	SEXP ans = PROTECT(allocMatrix(REALSXP, x_ncol, 2));
	for (int j = 0; j < x_ncol; j++) {
		int offset = INTEGER(x_slotp)[j];
		int nzcount = INTEGER(x_slotp)[j + 1] - offset;
		minmax_double(REAL(x_slotx) + offset, nzcount, narm,
			      nzcount < x_nrow,
			      REAL(ans) + j, REAL(ans) + x_ncol + j);
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * colVars()
 */

static double col_sum(const double *x, int x_len, int nrow, int narm,
		      int *sample_size)
{
	*sample_size = nrow;
	double sum = 0.0;
	for (int i = 0; i < x_len; i++) {
		double xi = x[i];
		/* ISNAN(): True for *both* NA and NaN. See <R_ext/Arith.h> */
		if (narm && ISNAN(xi)) {
			(*sample_size)--;
			continue;
		}
		sum += xi;
	}
	return sum;
}

static double col_var(const double *x, int x_len, int nrow, int narm)
{
	int sample_size;
	double sum = col_sum(x, x_len, nrow, narm, &sample_size);
	double mean = sum / (double) sample_size;
	double sigma = mean * mean * (nrow - x_len);
	for (int i = 0; i < x_len; i++) {
		double xi = x[i];
		if (narm && ISNAN(xi))
			continue;
		double delta = xi - mean;
		sigma += delta * delta;
	}
	return sigma / (sample_size - 1.0);
}

/* --- .Call ENTRY POINT ---
   About 2.5x faster than the colVars() method for dgCMatrix objects defined
   in the sparseMatrixStats package. */
SEXP C_colVars_dgCMatrix(SEXP x, SEXP na_rm)
{
	SEXP x_Dim = GET_SLOT(x, install("Dim"));
	int x_nrow = INTEGER(x_Dim)[0];
	int x_ncol = INTEGER(x_Dim)[1];
	SEXP x_slotx = GET_SLOT(x, install("x"));
	SEXP x_slotp = GET_SLOT(x, install("p"));
	int narm = LOGICAL(na_rm)[0];

	SEXP ans = PROTECT(NEW_NUMERIC(x_ncol));
	for (int j = 0; j < x_ncol; j++) {
		int offset = INTEGER(x_slotp)[j];
		int nzcount = INTEGER(x_slotp)[j + 1] - offset;
		REAL(ans)[j] = col_var(REAL(x_slotx) + offset, nzcount,
				       x_nrow, narm);
	}
	UNPROTECT(1);
	return ans;
}

