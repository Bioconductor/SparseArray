/****************************************************************************
 *                      Dot product of sparse vectors                       *
 ****************************************************************************/
#include "SparseVec_dotprod.h"

#include "SparseVec.h"


double _dotprod_doubleSV_doubleSV(const SparseVec *sv1, const SparseVec *sv2)
{
	int k1, k2, off;
	double val1, val2;

	double ans = 0.0;
	k1 = k2 = 0;
	while (next_nzvals_double_double(sv1, sv2,
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
   If not sure, use _dotprod_doubleSV_doubles() below.
   IMPORTANT: 'sv1->nzoffs' is assumed to contain valid offsets in 'x2'.
   This is NOT checked! */
double _dotprod_doubleSV_finite_doubles(const SparseVec *sv1, const double *x2)
{
	const double *nzvals1 = get_doubleSV_nzvals(sv1);
	int nzcount1 = get_SV_nzcount(sv1);
	double ans = 0.0;
	for (int k1 = 0; k1 < nzcount1; k1++)
		ans += nzvals1[k1] * x2[sv1->nzoffs[k1]];
	return ans;
}

/* Like _dotprod_doubleSV_finite_doubles() above but makes no assumptions
   about the content of 'x2'.
   Significantly slower than _dotprod_doubleSV_finite_doubles(). */
double _dotprod_doubleSV_doubles(const SparseVec *sv1, const double *x2)
{
	const double *nzvals1 = get_doubleSV_nzvals(sv1);
	int nzcount1 = get_SV_nzcount(sv1);
	double ans = 0.0;
	int k1 = 0;
	for (int i2 = 0; i2 < sv1->len; i2++) {
		double v1, v2 = x2[i2];
		if (R_IsNA(v2))
			return NA_REAL;
		if (k1 < nzcount1 && sv1->nzoffs[k1] == i2) {
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
   If not sure, use _dotprod_intSV_ints() below.
   IMPORTANT: 'sv1->nzoffs' is assumed to contain valid offsets in 'x2'.
   This is NOT checked! */
double _dotprod_intSV_noNA_ints(const SparseVec *sv1, const int *x2)
{
	const int *nzvals1 = get_intSV_nzvals(sv1);
	int nzcount1 = get_SV_nzcount(sv1);
	double ans = 0.0;
	for (int k1 = 0; k1 < nzcount1; k1++) {
		int v1 = nzvals1[k1];
		if (v1 == NA_INTEGER)
			return NA_REAL;
		ans += (double) v1 * x2[sv1->nzoffs[k1]];
	}
	return ans;
}

/* 'sv1' must contain ints.
   Like _dotprod_intSV_noNA_ints() above but makes no assumptions about the
   content of 'x2'. Significantly slower than _dotprod_intSV_noNA_ints(). */
double _dotprod_intSV_ints(const SparseVec *sv1, const int *x2)
{
	const int *nzvals1 = get_intSV_nzvals(sv1);
	int nzcount1 = get_SV_nzcount(sv1);
	double ans = 0.0;
	int k1 = 0;
	for (int i2 = 0; i2 < sv1->len; i2++) {
		int v1, v2 = x2[i2];
		if (v2 == NA_INTEGER)
			return NA_REAL;
		if (k1 < nzcount1 && sv1->nzoffs[k1] == i2) {
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

double _dotprod_doubles_zero(const double *x, int x_len)
{
	double ans = 0.0;
	for (int i = 0; i < x_len; i++) {
		double v = x[i];
		if (R_IsNA(v))
			return NA_REAL;
		ans += v * 0.0;  /* 'v' could be Inf or NaN */
	}
	return ans;
}

double _dotprod_ints_zero(const int *x, int x_len)
{
	double ans = 0.0;
	for (int i = 0; i < x_len; i++) {
		int v = x[i];
		if (v == NA_INTEGER)
			return NA_REAL;
		ans += (double) v * 0.0;
	}
	return ans;
}

double _dotprod_doubleSV_zero(const SparseVec *sv)
{
	return _dotprod_doubles_zero(get_doubleSV_nzvals(sv),
				     get_SV_nzcount(sv));
}

double _dotprod_intSV_zero(const SparseVec *sv)
{
	return _dotprod_ints_zero(get_intSV_nzvals(sv), get_SV_nzcount(sv));
}

