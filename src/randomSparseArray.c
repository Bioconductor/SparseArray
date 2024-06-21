/****************************************************************************
 *                   Generate a random SparseArray object                   *
 ****************************************************************************/
#include "randomSparseArray.h"

#include "leaf_utils.h"

#include <R_ext/Random.h>
#include <math.h>  /* for exp() */
#include <string.h>  /* for memcpy() */


/****************************************************************************
 * C_simple_rpois()
 */

#define	CUMSUM_DPOIS_MAX_LENGTH 32  /* enough to support 0 <= lambda <= 4 */

static int compute_cumsum_dpois(double *cumsum_dpois, double lambda)
{
	long double csdp, dp;
	int n;

	csdp = dp = exp(-lambda);
	if (csdp >= 1.0)  /* means 'lambda' is zero or very small */
		return 0;
	cumsum_dpois[0] = (double) csdp;
	//printf("cumsum_dpois[%d]=%2.12f\n", 0, cumsum_dpois[0]);
	for (n = 1; n < CUMSUM_DPOIS_MAX_LENGTH; n++) {
		dp *= lambda / n;
		csdp += dp;
		if ((double) csdp == cumsum_dpois[n - 1])
			return n;
		cumsum_dpois[n] = (double) csdp;
		//printf("cumsum_dpois[%d]=%2.12f\n", n, cumsum_dpois[n]);
	}
	return -1;
}

/* 'ticks' must contain positive values in strictly ascending order.
   Returns an integer >= 0 and <= 'nticks'. */
static inline int binary_find_interval(double u,
		const double *ticks, int nticks)
{
	int k1, k2, k;
	double csdp;

	if (nticks == 0)
		return 0;

	/* Compare with 'cumsum_dpois[0]'. */
	if (u < ticks[0])
		return 0;

	/* Compare with 'cumsum_dpois[cumsum_dpois_len-1]'. */
	k2 = nticks - 1;
	csdp = ticks[k2];
	if (u >= csdp)
		return nticks;

	/* Binary search.
	   Seems that using >> 1 instead of / 2 is faster, even when compiling
	   with 'gcc -O2' (one would hope that the optimizer is able to do that
	   kind of optimization). */
	k1 = 0;
	while ((k = (k1 + k2) >> 1) != k1) {
		csdp = ticks[k];
		if (u < csdp)
			k2 = k;
		else
			k1 = k;
	}
	return k2;
}

/* A linear search is more efficient than a binary search for our use case
   (Poisson distribution with small lambda value) so do not try to replace
   this with a binary search.
   'ticks' must contain positive values in strictly ascending order.
   Returns an integer >= 0 and <= 'nticks'. */
static inline int find_interval(double u, const double *ticks, int nticks)
{
	int k;

	for (k = 0; k < nticks; k++)
		if (u < ticks[k])
			break;
	return k;
}

static int simple_rpois(double lambda)
{
	static double last_lambda = -1;
	static double cumsum_dpois[CUMSUM_DPOIS_MAX_LENGTH];
	static int cumsum_dpois_len;

	double u;

	if (lambda != last_lambda) {
		cumsum_dpois_len = compute_cumsum_dpois(cumsum_dpois, lambda);
		//printf("cumsum_dpois_len = %d\n", cumsum_dpois_len);
		if (cumsum_dpois_len < 0)
			error("'lambda' too big?");
		last_lambda = lambda;
	}
	u = unif_rand();
	//return binary_find_interval(u, cumsum_dpois, cumsum_dpois_len);
	return find_interval(u, cumsum_dpois, cumsum_dpois_len);
}

/* --- .Call ENTRY POINT --- */
SEXP C_simple_rpois(SEXP n, SEXP lambda)
{
	int n0, i;
	double lambda0;
	SEXP ans;

	if (!IS_INTEGER(n) || LENGTH(n) != 1)
		error("'n' must be a single integer");
	n0 = INTEGER(n)[0];
	if (n0 < 0)
		error("'n' cannot be negative");

	if (!IS_NUMERIC(lambda) || LENGTH(lambda) != 1)
		error("'lambda' must be a single numeric value");
	lambda0 = REAL(lambda)[0];
	if (lambda0 < 0.0)
		error("'lambda' cannot be negative");

	ans = PROTECT(NEW_INTEGER(n0));
	GetRNGstate();
	for (i = 0; i < n0; i++)
		INTEGER(ans)[i] = simple_rpois(lambda0);
	PutRNGstate();
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_poissonSparseArray()
 */

static SEXP build_poisson_leaf(int dim0, double lambda,
			       int *nzvals_buf, int *nzoffs_buf)
{
	int buf_len = 0;
	for (int i = 0; i < dim0; i++) {
		int val = simple_rpois(lambda);
		if (val != 0) {
			nzvals_buf[buf_len] = val;
			nzoffs_buf[buf_len] = i;
			buf_len++;
		}
	}
	return _make_leaf_from_two_arrays(INTSXP,
					  nzvals_buf, nzoffs_buf, buf_len);
}

/* Recursive. */
static SEXP REC_build_poisson_SVT(const int *dim, int ndim, double lambda,
				  int *nzvals_buf, int *nzoffs_buf)
{
	int SVT_len, is_empty, i;
	SEXP ans, ans_elt;

	if (ndim == 1)
		return build_poisson_leaf(dim[0], lambda,
					  nzvals_buf, nzoffs_buf);

	SVT_len = dim[ndim - 1];  /* cannot be 0 */
	ans = PROTECT(NEW_LIST(SVT_len));
	is_empty = 1;
	for (i = 0; i < SVT_len; i++) {
		ans_elt = REC_build_poisson_SVT(dim, ndim - 1, lambda,
						nzvals_buf, nzoffs_buf);
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

/* --- .Call ENTRY POINT --- */
SEXP C_poissonSparseArray(SEXP dim, SEXP lambda)
{
	if (!IS_NUMERIC(lambda) || LENGTH(lambda) != 1)
		error("'lambda' must be a single numeric value");
	double lambda0 = REAL(lambda)[0];
	if (lambda0 < 0.0 || lambda0 > 4.0)
		error("'lambda' must be >= 0 and <= 4");
	if (lambda0 == 0.0)
		return R_NilValue;

	const int *dim_p = INTEGER(dim);
	int ndim = LENGTH(dim);
	for (int along = 0; along < ndim; along++)
		if (dim_p[along] == 0)
			return R_NilValue;

	int *nzvals_buf = (int *) R_alloc(dim_p[0], sizeof(int));
	int *nzoffs_buf = (int *) R_alloc(dim_p[0], sizeof(int));
	GetRNGstate();
	SEXP ans = PROTECT(
		REC_build_poisson_SVT(dim_p, ndim, lambda0,
				      nzvals_buf, nzoffs_buf)
	);
	PutRNGstate();
	UNPROTECT(1);
	return ans;
}

