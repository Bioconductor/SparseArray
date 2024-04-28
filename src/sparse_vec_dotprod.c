/****************************************************************************
 *                      Dot product of sparse vectors                       *
 ****************************************************************************/
#include "sparse_vec_dotprod.h"

#include "leaf_vector_utils.h"


double _dotprod_sparse_vecs(const struct sparse_vec *sv1,
			    const struct sparse_vec *sv2)
{
	int k1, k2, off;
	double val1, val2;

	const double *nzvals1 = _get_double_nzvals(sv1);
	const double *nzvals2 = _get_double_nzvals(sv2);
	double ans = 0.0;
	k1 = k2 = 0;
	while (next_nzval_double_double(
		sv1->nzoffs, nzvals1, sv1->nzcount,
		sv2->nzoffs, nzvals2, sv2->nzcount,
		&k1, &k2, &off, &val1, &val2))
	{
		if (R_IsNA(val1) || R_IsNA(val2))
			return NA_REAL;
		ans += val1 * val2;
	}
	return ans;
}

/* Safe to use only if 'x2' is finite i.e. contains no NA, NaN, Inf, or -Inf,
   or if 'sv1' and 'x2' represent the same numeric vector (in sparse and
   dense form, respectively).
   If not sure, use _dotprod_sparse_vec_and_double_col() below.
   IMPORTANT: 'sv1->nzoffs' is assumed to contain valid offsets in 'x2'.
   This is NOT checked! */
double _dotprod_sparse_vec_and_finite_col(const struct sparse_vec *sv1,
					  const double *x2)
{
	const double *nzvals1 = _get_double_nzvals(sv1);
	double ans = 0.0;
	for (int k1 = 0; k1 < sv1->nzcount; k1++)
		ans += nzvals1[k1] * x2[sv1->nzoffs[k1]];
	return ans;
}

/* Like _dotprod_sparse_vec_and_finite_col() above but makes no
   assumptions about the content of 'x2'.
   Significantly slower than _dotprod_sparse_vec_and_finite_col(). */
double _dotprod_sparse_vec_and_double_col(const struct sparse_vec *sv1,
		const double *x2, int x2_len)
{
	const double *nzvals1 = _get_double_nzvals(sv1);
	double ans = 0.0;
	int k1 = 0;
	for (int i2 = 0; i2 < x2_len; i2++) {
		double v1, v2 = x2[i2];
		if (R_IsNA(v2))
			return NA_REAL;
		if (k1 < sv1->nzcount && sv1->nzoffs[k1] == i2) {
			v1 = nzvals1[k1];
			if (R_IsNA(v1))
				return NA_REAL;
			k1++;
		} else {
			v1 = 0.0;
		}
		ans += v1 * v2;
	}
	return ans;
}

/* 'sv1' must contain ints.
   Safe to use only if 'x2' contains no NA, or if 'sv1' and 'x2' represent
   the same integer vector (in sparse and dense form, respectively).
   If not sure, use _dotprod_sparse_vec_and_int_col() below.
   IMPORTANT: 'sv1->nzoffs' is assumed to contain valid offsets in 'x2'.
   This is NOT checked! */
double _dotprod_sparse_vec_and_noNA_int_col(const struct sparse_vec *sv1,
					    const int *x2)
{
	const int *nzvals1 = _get_int_nzvals(sv1);
	double ans = 0.0;
	for (int k1 = 0; k1 < sv1->nzcount; k1++) {
		int v1 = nzvals1[k1];
		if (v1 == NA_INTEGER)
			return NA_REAL;
		ans += (double) v1 * x2[sv1->nzoffs[k1]];
	}
	return ans;
}

/* 'sv1' must contain ints.
   Like _dotprod_sparse_vec_and_noNA_int_col() above but makes no
   assumptions about the content of 'x2'.
   Significantly slower than _dotprod_sparse_vec_and_noNA_int_col(). */
double _dotprod_sparse_vec_and_int_col(const struct sparse_vec *sv1,
		const int *x2, int x2_len)
{
	const int *nzvals1 = _get_int_nzvals(sv1);
	double ans = 0.0;
	int k1 = 0;
	for (int i2 = 0; i2 < x2_len; i2++) {
		int v1, v2 = x2[i2];
		if (v2 == NA_INTEGER)
			return NA_REAL;
		if (k1 < sv1->nzcount && sv1->nzoffs[k1] == i2) {
			v1 = nzvals1[k1];
			if (v1 == NA_INTEGER)
				return NA_REAL;
			k1++;
		} else {
			v1 = 0;
		}
		ans += (double) v1 * v2;
	}
	return ans;
}

double _dotprod0_int_col(const int *x, int x_len)
{
	double ans;
	int i, v;

	ans = 0.0;
	for (i = 0; i < x_len; i++) {
		v = x[i];
		if (v == NA_INTEGER)
			return NA_REAL;
		ans += (double) v * 0.0;
	}
	return ans;
}

double _dotprod0_double_col(const double *x, int x_len)
{
	double ans, v;
	int i;

	ans = 0.0;
	for (i = 0; i < x_len; i++) {
		v = x[i];
		if (R_IsNA(v))
			return NA_REAL;
		ans += v * 0.0;  /* 'v' could be Inf or NaN */
	}
	return ans;
}

double _dotprod0_sparse_vec(const struct sparse_vec *sv)
{
	return _dotprod0_double_col(_get_double_nzvals(sv), sv->nzcount);
}

