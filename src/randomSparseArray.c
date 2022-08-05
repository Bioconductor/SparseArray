/****************************************************************************
 *                   Generate a random SparseArray object                   *
 ****************************************************************************/
#include "randomSparseArray.h"

#include "leaf_vector_utils.h"

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

/* Returns a value >= 0 and <= 'cumsum_dpois_len'. */
static inline int bsearch_cumsum_dpois(double u,
		const double *cumsum_dpois, int cumsum_dpois_len)
{
	int k1, k2, k;
	double csdp;

	if (cumsum_dpois_len == 0)
		return 0;

	/* Compare with 'cumsum_dpois[0]'. */
	if (u < cumsum_dpois[0])
		return 0;

	/* Compare with 'cumsum_dpois[cumsum_dpois_len-1]'. */
	k2 = cumsum_dpois_len - 1;
	csdp = cumsum_dpois[k2];
	if (u >= csdp)
		return cumsum_dpois_len;

	/* Binary search.
	   Seems that using >> 1 instead of / 2 is faster, even when compiling
	   with 'gcc -O2' (one would hope that the optimizer is able to do that
	   kind of optimization). */
	k1 = 0;
	while ((k = (k1 + k2) >> 1) != k1) {
		csdp = cumsum_dpois[k];
		if (u < csdp)
			k2 = k;
		else
			k1 = k;
	}
	return k2;
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
/*
	u = 0.0;
	int k = bsearch_cumsum_dpois(u, cumsum_dpois, cumsum_dpois_len);
	printf("u = %2.12f --> k = %d\n", u, k);
	for (int n = 0; n < cumsum_dpois_len; n++) {
		printf("cumsum_dpois[%d] = %2.12f\n", n, cumsum_dpois[n]);
		u = cumsum_dpois[n] - 0.001;
		k = bsearch_cumsum_dpois(u, cumsum_dpois, cumsum_dpois_len);
		printf("u = %2.12f --> k = %d\n", u, k);
		u = cumsum_dpois[n];
		k = bsearch_cumsum_dpois(u, cumsum_dpois, cumsum_dpois_len);
		printf("u = %2.12f --> k = %d\n", u, k);
		u = cumsum_dpois[n] + 0.001;
		k = bsearch_cumsum_dpois(u, cumsum_dpois, cumsum_dpois_len);
		printf("u = %2.12f --> k = %d\n", u, k);
	}
*/
	u = unif_rand();
	return bsearch_cumsum_dpois(u, cumsum_dpois, cumsum_dpois_len);
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

/* Returns R_NilValue or a "leaf vector". */
static SEXP build_poisson_leaf_vector(int d1, double lambda,
				      int *offs_buf, int *vals_buf)
{
	int ans_len, i, val;
	SEXP ans_offs, ans_vals, ans;

	ans_len = 0;
	for (i = 0; i < d1; i++) {
		val = simple_rpois(lambda);
		if (val != 0) {
			offs_buf[ans_len] = i;
			vals_buf[ans_len] = val;
			ans_len++;
		}
	}

	if (ans_len == 0)
		return R_NilValue;

	ans_offs = PROTECT(NEW_INTEGER(ans_len));
	memcpy(INTEGER(ans_offs), offs_buf, sizeof(int) * ans_len);
	ans_vals = PROTECT(NEW_INTEGER(ans_len));
	memcpy(INTEGER(ans_vals), vals_buf, sizeof(int) * ans_len);
	ans = _new_leaf_vector(ans_offs, ans_vals);
	UNPROTECT(2);
	return ans;
}

/* Recursive. */
static SEXP REC_build_poisson_SVT(const int *dim, int ndim, double lambda,
				  int *offs_buf, int *vals_buf)
{
	int SVT_len, is_empty, i;
	SEXP ans, ans_elt;

	if (ndim == 1)
		return build_poisson_leaf_vector(dim[0], lambda,
						 offs_buf, vals_buf);

	SVT_len = dim[ndim - 1];  /* cannot be 0 */
	ans = PROTECT(NEW_LIST(SVT_len));
	is_empty = 1;
	for (i = 0; i < SVT_len; i++) {
		ans_elt = REC_build_poisson_SVT(dim, ndim - 1, lambda,
						offs_buf, vals_buf);
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
	double lambda0;
	const int *dim_p;
	int ndim, along, *offs_buf, *vals_buf;
	SEXP ans;

	if (!IS_NUMERIC(lambda) || LENGTH(lambda) != 1)
		error("'lambda' must be a single numeric value");
	lambda0 = REAL(lambda)[0];
	if (lambda0 < 0.0 || lambda0 > 4.0)
		error("'lambda' must be >= 0 and <= 4");

	dim_p = INTEGER(dim);
	ndim = LENGTH(dim);
	for (along = 0; along < ndim; along++)
		if (dim_p[along] == 0)
			return R_NilValue;

	offs_buf = (int *) R_alloc(dim_p[0], sizeof(int));
	vals_buf = (int *) R_alloc(dim_p[0], sizeof(int));
	GetRNGstate();
	ans = PROTECT(
		REC_build_poisson_SVT(dim_p, ndim, lambda0, offs_buf, vals_buf)
	);
	PutRNGstate();
	UNPROTECT(1);
	return ans;
}

