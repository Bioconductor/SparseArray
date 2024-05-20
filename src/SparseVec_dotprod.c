/****************************************************************************
 *                      Dot product of sparse vectors                       *
 ****************************************************************************/
#include "SparseVec_dotprod.h"

#include "SparseVec.h"


double _dotprod_doubleSV_doubleSV(const SparseVec *sv1, const SparseVec *sv2)
{
	double ans = 0.0, val1, val2;
	int k1 = 0, k2 = 0, off;
	while (next_2SV_vals_double_double(sv1, sv2,
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
	double ans = 0.0;
	int nzcount1 = get_SV_nzcount(sv1);
	if (sv1->nzvals == R_NilValue) {  /* lacunar SparseVec */
		for (int k1 = 0; k1 < nzcount1; k1++)
			ans += x2[sv1->nzoffs[k1]];
	} else {  /* regular SparseVec */
		const double *nzvals1_p = get_doubleSV_nzvals_p(sv1);
		for (int k1 = 0; k1 < nzcount1; k1++)
			ans += nzvals1_p[k1] * x2[sv1->nzoffs[k1]];
	}
	return ans;
}

/* Like _dotprod_doubleSV_finite_doubles() above but makes no assumptions
   about the content of 'x2'.
   Significantly slower than _dotprod_doubleSV_finite_doubles(). */
double _dotprod_doubleSV_doubles(const SparseVec *sv1, const double *x2)
{
	double ans = 0.0;
	int k1 = 0;
	for (int i2 = 0; i2 < sv1->len; i2++) {
		double val1 = double0, val2 = x2[i2];
		if (R_IsNA(val2))
			return NA_REAL;
		if (k1 < get_SV_nzcount(sv1) && sv1->nzoffs[k1] == i2) {
			val1 = get_doubleSV_nzval(sv1, k1);
			if (R_IsNA(val1))
				return NA_REAL;
			k1++;
		}
		ans += val1 * val2;
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
	double ans = 0.0;
	int nzcount1 = get_SV_nzcount(sv1);
	if (sv1->nzvals == R_NilValue) {  /* lacunar SparseVec */
		for (int k1 = 0; k1 < nzcount1; k1++)
			ans += (double) x2[sv1->nzoffs[k1]];
	} else {  /* regular SparseVec */
		const int *nzvals1_p = get_intSV_nzvals_p(sv1);
		for (int k1 = 0; k1 < nzcount1; k1++) {
			int v1 = nzvals1_p[k1];
			if (v1 == NA_INTEGER)
				return NA_REAL;
			ans += (double) v1 * x2[sv1->nzoffs[k1]];
		}
	}
	return ans;
}

/* 'sv1' must contain ints.
   Like _dotprod_intSV_noNA_ints() above but makes no assumptions about the
   content of 'x2'. Significantly slower than _dotprod_intSV_noNA_ints(). */
double _dotprod_intSV_ints(const SparseVec *sv1, const int *x2)
{
	double ans = 0.0;
	int k1 = 0;
	for (int i2 = 0; i2 < sv1->len; i2++) {
		int val1 = int0, val2 = x2[i2];
		if (val2 == NA_INTEGER)
			return NA_REAL;
		if (k1 < get_SV_nzcount(sv1) && sv1->nzoffs[k1] == i2) {
			val1 = get_intSV_nzval(sv1, k1);
			if (val1 == NA_INTEGER)
				return NA_REAL;
			k1++;
		}
		ans += (double) val1 * val2;
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
		ans += v * double0;  /* 'v' could be Inf or NaN */
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
		ans += (double) v * double0;
	}
	return ans;
}

double _dotprod_doubleSV_zero(const SparseVec *sv)
{
	if (sv->nzvals == R_NilValue)  /* lacunar SparseVec */
		return 0.0;
	/* regular SparseVec */
	return _dotprod_doubles_zero(get_doubleSV_nzvals_p(sv),
				     get_SV_nzcount(sv));
}

double _dotprod_intSV_zero(const SparseVec *sv)
{
	if (sv->nzvals == R_NilValue)  /* lacunar SparseVec */
		return 0.0;
	/* regular SparseVec */
	return _dotprod_ints_zero(get_intSV_nzvals_p(sv), get_SV_nzcount(sv));
}

