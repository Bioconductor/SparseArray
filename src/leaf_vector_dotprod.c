/****************************************************************************
 *                      Dot product of "leaf vectors"                       *
 ****************************************************************************/
#include "leaf_vector_dotprod.h"

#include "leaf_vector_utils.h"


double _dotprod_leaf_vectors(SEXP lv1, SEXP lv2)
{
	int lv1_len, lv2_len, k1, k2, off;
	SEXP lv1_offs, lv1_vals, lv2_offs, lv2_vals;
	double ans, v1, v2;

	lv1_len = _split_leaf_vector(lv1, &lv1_offs, &lv1_vals);
	lv2_len = _split_leaf_vector(lv2, &lv2_offs, &lv2_vals);
	k1 = k2 = 0;
	ans = 0.0;
	while (next_nzval_double_double(
			INTEGER(lv1_offs), REAL(lv1_vals), lv1_len,
			INTEGER(lv2_offs), REAL(lv2_vals), lv2_len,
			&k1, &k2, &off, &v1, &v2))
	{
		if (R_IsNA(v1) || R_IsNA(v2))
			return NA_REAL;
		ans += v1 * v2;
	}
	return ans;
}

/* 'lv1' must be a leaf vector containing doubles.
   Safe to use if 'x2' is finite i.e. contains no NA, NaN, Inf, or -Inf,
   or if 'lv1' and 'x2' represent the same numeric vector.
   If not sure, use _dotprod_leaf_vector_and_double_col() below.
   The offsets in 'lv1' are assumed to be valid offsets in 'x2'. This
   is NOT checked! */
double _dotprod_leaf_vector_and_finite_col(SEXP lv1, const double *x2)
{
	int lv1_len, k1;
	SEXP lv1_offs, lv1_vals;
	const int *offs1_p;
	const double *vals1_p;
	double ans;

	lv1_len = _split_leaf_vector(lv1, &lv1_offs, &lv1_vals);
	offs1_p = INTEGER(lv1_offs);
	vals1_p = REAL(lv1_vals);
	ans = 0.0;
	for (k1 = 0; k1 < lv1_len; k1++) {
		ans += *vals1_p * x2[*offs1_p];
		offs1_p++;
		vals1_p++;
	}
	return ans;
}

/* Makes no assumption about the content of 'x2'.
   Significantly slower than _dotprod_leaf_vector_and_finite_col() above. */
double _dotprod_leaf_vector_and_double_col(SEXP lv1,
		const double *x2, int x2_len)
{
	int lv1_len, k1;
	SEXP lv1_offs, lv1_vals;
	const int *offs1_p;
	const double *vals1_p;
	double ans, v1, v2;
	int i2;

	lv1_len = _split_leaf_vector(lv1, &lv1_offs, &lv1_vals);
	offs1_p = INTEGER(lv1_offs);
	vals1_p = REAL(lv1_vals);
	ans = 0.0;
	for (i2 = k1 = 0; i2 < x2_len; i2++) {
		v2 = x2[i2];
		if (R_IsNA(v2))
			return NA_REAL;
		if (k1 < lv1_len && i2 == offs1_p[k1]) {
			v1 = vals1_p[k1];
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

/* 'lv1' must be a leaf vector containing ints.
   Safe to use if 'x2' contains no NA, or if 'lv1' and 'x2' represent
   the same integer vector.
   If not sure, use _dotprod_leaf_vector_and_int_col() below.
   The offsets in 'lv1' are assumed to be valid offsets in 'x2'. This
   is NOT checked! */
double _dotprod_leaf_vector_and_noNA_int_col(SEXP lv1, const int *x2)
{
	int lv1_len, k1;
	SEXP lv1_offs, lv1_vals;
	const int *offs1_p;
	const int *vals1_p;
	double ans;
	int v1;

	lv1_len = _split_leaf_vector(lv1, &lv1_offs, &lv1_vals);
	offs1_p = INTEGER(lv1_offs);
	vals1_p = INTEGER(lv1_vals);
	ans = 0.0;
	for (k1 = 0; k1 < lv1_len; k1++) {
		v1 = *vals1_p;
		if (v1 == NA_INTEGER)
			return NA_REAL;
		ans += (double) v1 * x2[*offs1_p];
		offs1_p++;
		vals1_p++;
	}
	return ans;
}

/* Makes no assumption about the content of 'x2'.
   Significantly slower than _dotprod_leaf_vector_and_noNA_int_col() above. */
double _dotprod_leaf_vector_and_int_col(SEXP lv1,
		const int *x2, int x2_len)
{
	SEXP lv1_offs, lv1_vals;
	const int *offs1_p, *vals1_p;
	double ans;
	int i2, v1, v2;

	_split_leaf_vector(lv1, &lv1_offs, &lv1_vals);
	offs1_p = INTEGER(lv1_offs);
	vals1_p = INTEGER(lv1_vals);
	ans = 0.0;
	for (i2 = 0; i2 < x2_len; i2++) {
		v2 = x2[i2];
		if (v2 == NA_INTEGER)
			return NA_REAL;
		if (*offs1_p > i2) {
			v1 = 0;
		} else {
			v1 = *vals1_p;
			if (v1 == NA_INTEGER)
				return NA_REAL;
			offs1_p++;
			vals1_p++;
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

double _dotprod0_leaf_vector(SEXP lv)
{
	int lv_len;
	SEXP lv_offs, lv_vals;

	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	return _dotprod0_double_col(REAL(lv_vals), lv_len);
}

