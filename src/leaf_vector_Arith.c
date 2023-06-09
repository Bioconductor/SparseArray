/****************************************************************************
 *                   'Arith' operations on "leaf vectors"                   *
 ****************************************************************************/
#include "leaf_vector_Arith.h"

#include "leaf_vector_utils.h"

#include <limits.h>  /* for INT_MAX */
#include <math.h>    /* for trunc(), pow(), fmod(), floor() */


int _get_Arith_opcode(SEXP op)
{
	const char *s;

	if (!IS_CHARACTER(op) || LENGTH(op) != 1)
		error("SparseArray internal error in _get_Arith_opcode():\n"
		      "    'op' must be a single string");
	op = STRING_ELT(op, 0);
	if (op == NA_STRING)
		error("SparseArray internal error in _get_Arith_opcode():\n"
		      "    'op' cannot be NA");
	s = CHAR(op);
	if (strcmp(s, "+") == 0)
		return ADD_OPCODE;
	if (strcmp(s, "-") == 0)
		return SUB_OPCODE;
	if (strcmp(s, "*") == 0)
		return MULT_OPCODE;
	if (strcmp(s, "/") == 0)
		return DIV_OPCODE;
	if (strcmp(s, "^") == 0)
		return POW_OPCODE;
	if (strcmp(s, "%%") == 0)
		return MOD_OPCODE;
	if (strcmp(s, "%/%") == 0)
		return IDIV_OPCODE;
	error("SparseArray internal error in _get_Arith_opcode():\n"
	      "    invalid op: \"%s\"", s);
	return 0;  /* will never reach this */
}

/* Assumes that 'ans_Rtype' is equal or bigger than the type of the
   vals in 'lv'. Performs **in-place** replacement if 'ans_Rtype' is 0!
   Note that this could also have been achieved by calling:

     _Arith_lv1_v2(lv, -1, MULT_OPCODE, ...)

   but _unary_minus_leaf_vector() takes a lot of shortcuts so is way
   more efficient. */
SEXP _unary_minus_leaf_vector(SEXP lv, SEXPTYPE ans_Rtype)
{
	int lv_len, supported, k;
	SEXP lv_offs, lv_vals, ans_vals, ans;

	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	if (ans_Rtype == 0) {
		/* In-place replacement! */
		ans_vals = lv_vals;
	} else {
		ans_vals = PROTECT(allocVector(ans_Rtype, lv_len));
	}
	supported = 0;
	if (TYPEOF(lv_vals) == INTSXP) {
		if (ans_Rtype == INTSXP || ans_Rtype == 0) {
			const int *lv_vals_p = INTEGER(lv_vals);
			int *ans_vals_p = INTEGER(ans_vals);
			int v;
			for (k = 0; k < lv_len; k++) {
				v = lv_vals_p[k];
				if (v != NA_INTEGER)
					v = -v;
				ans_vals_p[k] = v;
			}
			supported = 1;
		} else if (ans_Rtype == REALSXP) {
			const int *lv_vals_p = INTEGER(lv_vals);
			double *ans_vals_p = REAL(ans_vals);
			int v;
			for (k = 0; k < lv_len; k++) {
				v = lv_vals_p[k];
				if (v == NA_INTEGER) {
					ans_vals_p[k] = NA_REAL;
				} else {
					v = -v;
					ans_vals_p[k] = (double) v;
				}
			}
			supported = 1;
		}
	} else if (TYPEOF(lv_vals) == REALSXP) {
		if (ans_Rtype == REALSXP || ans_Rtype == 0) {
			const double *lv_vals_p = REAL(lv_vals);
			double *ans_vals_p = REAL(ans_vals);
			for (k = 0; k < lv_len; k++) {
				ans_vals_p[k] = - lv_vals_p[k];
			}
			supported = 1;
		}
	}
	if (!supported) {
		if (ans_Rtype != 0)
			UNPROTECT(1);
		error("_unary_minus_leaf_vector() only supports input "
		      "of type \"integer\" or \"double\" at the moment");
	}
	if (ans_Rtype == 0)
		return lv;
	ans = _new_leaf_vector(lv_offs, ans_vals);
	UNPROTECT(1);
	return ans;
}

/* Does not support 'opcode' values DIV_OPCODE ("/") or POW_OPCODE ("^"). */
static inline int Arith_int(int x, int y, int opcode, int *ovflow)
{
	double vv;
	int zz;

	if (x == NA_INTEGER || y == NA_INTEGER)
		return NA_INTEGER;
	switch (opcode) {
	    case ADD_OPCODE:  vv = (double) x + y; break;
	    case SUB_OPCODE:  vv = (double) x - y; break;
	    case MULT_OPCODE: vv = (double) x * y; break;
	    case MOD_OPCODE:
		if (y == 0)
			return NA_INTEGER;
		zz = x % y;
		/* R's "%%" wants the result to be the same sign as 'y' (this
		   deviates from C's modulo operator %), so we adjust to
		   provide R's "%%" behavior. */
		if ((y > 0 && zz < 0) || (y < 0 && zz > 0))
			zz += y;
		return zz;
	    case IDIV_OPCODE:
		if (y == 0)
			return NA_INTEGER;
		zz = x / y;
		if (((y > 0 && x < 0) || (y < 0 && x > 0)) && zz * y != x)
			zz--;
		return zz;
	    default:
		error("SparseArray internal error in Arith_int():\n"
		      "    unsupported 'opcode'");
	}
	if (vv <= INT_MIN || vv > INT_MAX) {
		*ovflow = 1;
		return NA_INTEGER;
	}
	return (int) vv;
}

/* When x is negative, the power operator in R ("^") and the pow() function
   from the C standard library behave differently. Two differences:
   a. If y is R_PosInf, then x ^ y always returns NaN (even when -1 < x < 0),
      whereas pow(x, y) returns R_PosInf if x < -1, o1r 1 if x == -1, or 0
      if x > -1.
   b. x ^ y is expected to return NaN for any noninteger y. However, in
      the specific case where x is R_NegInf, pow(x, y) returns R_PosInf
      for any noninteger y on my Intel Ubuntu 22.04 laptop. */
static inline double Rpow_double(double x, double y)
{
	if (x < 0 && y == R_PosInf || x == R_NegInf && y != trunc(y))
		return R_NaN;
	return pow(x, y);
}

static inline double Rmod_double(double x, double y)
{
	if (y == 0.0 || x == R_PosInf || x == R_NegInf)
		return R_NaN;
	if (x == 0.0)
		return 0.0;
	if (y == R_PosInf)
		return ISNAN(x) || x > 0 ? x : R_PosInf;
	if (y == R_NegInf)
		return ISNAN(x) || x < 0 ? x : R_NegInf;
	return x - y * floor(x / y);
}

static inline double Ridiv_double(double x, double y)
{
	if (y == R_PosInf) {
		if (x == R_NegInf)
			return R_NaN;
		if (x < 0)
			return -1.0;
	} else if (y == R_NegInf) {
		if (x == R_PosInf)
			return R_NaN;
		if (x > 0)
			return -1.0;
	}
	return floor(x / y);
}

static inline double Arith_double(double x, double y, int opcode)
{
	switch (opcode) {
	    case ADD_OPCODE:  return x + y;
	    case SUB_OPCODE:  return x - y;
	    case MULT_OPCODE: return x * y;
	    case DIV_OPCODE:  return x / y;
	    case POW_OPCODE:  return Rpow_double(x, y);
	    case MOD_OPCODE:  return Rmod_double(x, y);
	    case IDIV_OPCODE: return Ridiv_double(x, y);
	}
	error("SparseArray internal error in Arith_double():\n"
	      "    unsupported 'opcode'");
	return 0.0;  /* will never reach this */
}

static int sparse_Arith_ints_1int(
		const int *offs1, const int *vals1, int n1, int v2,
		int opcode, int *offs_buf, int *vals_buf, int *ovflow)
{
	int ans_len, k, v;

	for (ans_len = k = 0; k < n1; k++) {
		v = Arith_int(vals1[k], v2, opcode, ovflow);
		if (v != 0) {
			offs_buf[ans_len] = offs1[k];
			vals_buf[ans_len] = v;
			ans_len++;
		}
	}
	return ans_len;
}

static int sparse_Arith_ints_ints(
		const int *offs1, const int *vals1, int n1,
		const int *offs2, const int *vals2, int n2,
		int opcode, int *offs_buf, int *vals_buf, int *ovflow)
{
	int ans_len, k1, k2, off, v1, v2, v;

	ans_len = k1 = k2 = 0;
	while (next_nzval_int_int(offs1, vals1, n1,
				  offs2, vals2, n2,
				  &k1, &k2, &off, &v1, &v2))
	{
		v = Arith_int(v1, v2, opcode, ovflow);
		if (v != 0) {
			offs_buf[ans_len] = off;
			vals_buf[ans_len] = v;
			ans_len++;
		}
	}
	return ans_len;
}

static int sparse_Arith_ints_1double(
		const int *offs1, const int *vals1, int n1, double v2,
		int opcode, int *offs_buf, double *vals_buf)
{
	int ans_len, k, v1;
	double v;

	for (ans_len = k = 0; k < n1; k++) {
		v1 = vals1[k];
		if (v1 == NA_INTEGER) {
			v = NA_REAL;
		} else {
			v = Arith_double((double) v1, v2, opcode);
		}
		if (v != 0) {
			offs_buf[ans_len] = offs1[k];
			vals_buf[ans_len] = v;
			ans_len++;
		}
	}
	return ans_len;
}

static int sparse_Arith_ints_doubles(
		const int *offs1, const int *vals1, int n1,
		const int *offs2, const double *vals2, int n2,
		int opcode, int *offs_buf, double *vals_buf)
{
	int ans_len, k1, k2, off, v1;
	double v2, v;

	ans_len = k1 = k2 = 0;
	while (next_nzval_int_double(offs1, vals1, n1,
				     offs2, vals2, n2,
				     &k1, &k2, &off, &v1, &v2))
	{
		if (v1 == NA_INTEGER) {
			v = NA_REAL;
		} else {
			v = Arith_double((double) v1, v2, opcode);
		}
		if (v != 0.0) {
			offs_buf[ans_len] = off;
			vals_buf[ans_len] = v;
			ans_len++;
		}
	}
	return ans_len;
}

static int sparse_Arith_doubles_ints(
		const int *offs1, const double *vals1, int n1,
		const int *offs2, const int *vals2, int n2,
		int opcode, int *offs_buf, double *vals_buf)
{
	int ans_len, k1, k2, off, v2;
	double v1, v;

	ans_len = k1 = k2 = 0;
	while (next_nzval_double_int(offs1, vals1, n1,
				     offs2, vals2, n2,
				     &k1, &k2, &off, &v1, &v2))
	{
		if (v2 == NA_INTEGER) {
			v = NA_REAL;
		} else {
			v = Arith_double(v1, (double) v2, opcode);
		}
		if (v != 0.0) {
			offs_buf[ans_len] = off;
			vals_buf[ans_len] = v;
			ans_len++;
		}
	}
	return ans_len;
}

static int sparse_Arith_doubles_1double(
		const int *offs1, const double *vals1, int n1, double v2,
		int opcode, int *offs_buf, double *vals_buf)
{
	int ans_len, k;
	double v;

	for (ans_len = k = 0; k < n1; k++) {
		v = Arith_double(vals1[k], v2, opcode);
		if (v != 0.0) {
			offs_buf[ans_len] = offs1[k];
			vals_buf[ans_len] = v;
			ans_len++;
		}
	}
	return ans_len;
}

static int sparse_Arith_doubles_doubles(
		const int *offs1, const double *vals1, int n1,
		const int *offs2, const double *vals2, int n2,
		int opcode, int *offs_buf, double *vals_buf)
{
	int ans_len, k1, k2, off;
	double v1, v2, v;

	ans_len = k1 = k2 = 0;
	while (next_nzval_double_double(offs1, vals1, n1,
					offs2, vals2, n2,
					&k1, &k2, &off, &v1, &v2))
	{
		v = Arith_double(v1, v2, opcode);
		if (v != 0.0) {
			offs_buf[ans_len] = off;
			vals_buf[ans_len] = v;
			ans_len++;
		}
	}
	return ans_len;
}

/* 'v2' is assumed to be an atomic vector of length 1. This is NOT checked! */
SEXP _Arith_lv1_v2(SEXP lv1, SEXP v2, int opcode, SEXPTYPE ans_Rtype,
		   int *offs_buf, void *vals_buf, int *ovflow)
{
	int lv1_len, ans_len;
	SEXP lv1_offs, lv1_vals;
	const int *offs1_p;
	SEXPTYPE buf_Rtype;

	lv1_len = _split_leaf_vector(lv1, &lv1_offs, &lv1_vals);
	offs1_p = INTEGER(lv1_offs);
	buf_Rtype = REALSXP;
	ans_len = -1;
	if (TYPEOF(lv1_vals) == INTSXP) {
		if (TYPEOF(v2) == INTSXP) {
			buf_Rtype = INTSXP;
			ans_len = sparse_Arith_ints_1int(
				offs1_p, INTEGER(lv1_vals), lv1_len,
				INTEGER(v2)[0],
				opcode, offs_buf, (int *) vals_buf, ovflow);
		} else if (TYPEOF(v2) == REALSXP) {
			ans_len = sparse_Arith_ints_1double(
				offs1_p, INTEGER(lv1_vals), lv1_len,
				REAL(v2)[0],
				opcode, offs_buf, (double *) vals_buf);
		}
	} else if (TYPEOF(lv1_vals) == REALSXP) {
		if (TYPEOF(v2) == REALSXP) {
			ans_len = sparse_Arith_doubles_1double(
				offs1_p, REAL(lv1_vals), lv1_len,
				REAL(v2)[0],
				opcode, offs_buf, (double *) vals_buf);
		}
	}
	if (ans_len == -1)
		error("_Arith_lv1_v2() only supports input of "
		      "type \"integer\" or \"double\" at the moment");
	if (ans_Rtype != buf_Rtype)
		error("SparseArray internal error in "
		      "_Arith_lv1_v2():\n"
		      "    ans_Rtype != buf_Rtype");
	return _make_leaf_vector_from_bufs(ans_Rtype,
					   offs_buf, vals_buf, ans_len);
}

/* Multiply the vals in 'lv' with zero. Will return R_NilValue if all
   the vals in 'lv' are finite (i.e. no NA, NaN, Inf, or -Inf).
   Note that this could simply be achieved by calling:

     _Arith_lv1_v2(lv, 0, MULT_OPCODE, ...)

   but mult0_leaf_vector() takes a lot of shortcuts so is way more
   efficient.
   Assumes that 'ans_Rtype' is equal or bigger than the type of the
   vals in 'lv'. */
static SEXP mult0_leaf_vector(SEXP lv, SEXPTYPE ans_Rtype,
			      int *offs_buf, void *vals_buf)
{
	int lv_len, ans_len, k;
	SEXP lv_offs, lv_vals;
	const int *lv_offs_p;

	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	lv_offs_p = INTEGER(lv_offs);
	ans_len = -1;
	if (TYPEOF(lv_vals) == INTSXP) {
		if (ans_Rtype == INTSXP) {
			/* We only keep NAs. */
			const int *lv_vals_p = INTEGER(lv_vals);
			int *ans_vals_p = (int *) vals_buf;
			int v;
			for (k = ans_len = 0; k < lv_len; k++) {
				v = lv_vals_p[k];
				if (v == NA_INTEGER) {
					offs_buf[ans_len] = lv_offs_p[k];
					ans_vals_p[ans_len] = NA_INTEGER;
					ans_len++;
				}
			}
		} else if (ans_Rtype == REALSXP) {
			/* We only keep NAs. */
			const int *lv_vals_p = INTEGER(lv_vals);
			double *ans_vals_p = (double *) vals_buf;
			int v;
			for (k = ans_len = 0; k < lv_len; k++) {
				v = lv_vals_p[k];
				if (v == NA_INTEGER) {
					offs_buf[ans_len] = lv_offs_p[k];
					ans_vals_p[ans_len] = NA_REAL;
					ans_len++;
				}
			}
		}
	} else if (TYPEOF(lv_vals) == REALSXP) {
		if (ans_Rtype == REALSXP) {
			ans_len = sparse_Arith_doubles_1double(
				lv_offs_p, REAL(lv_vals), lv_len, 0.0,
				MULT_OPCODE,
				offs_buf, (double *) vals_buf);
		}
	}
	if (ans_len == -1)
		error("mult0_leaf_vector() only supports input "
		      "of type \"integer\" or \"double\" at the moment");
	return _make_leaf_vector_from_bufs(ans_Rtype,
					   offs_buf, vals_buf, ans_len);
}

/* 'lv1' and 'lv2' must be "leaf vectors", with the following exceptions:
   - 'lv1' can be NULL if 'opcode' is SUB_OPCODE;
   - either 'lv1' or 'lv2' (but not both) can be NULL if 'opcode'
     is MULT_OPCODE. */
SEXP _Arith_lv1_lv2(SEXP lv1, SEXP lv2, int opcode, SEXPTYPE ans_Rtype,
		    int *offs_buf, void *vals_buf, int *ovflow)
{
	int lv1_len, lv2_len, ans_len;
	SEXP lv1_offs, lv1_vals, lv2_offs, lv2_vals;
	const int *offs1_p, *offs2_p;
	SEXPTYPE buf_Rtype;

	if (lv1 == R_NilValue) {
		if (lv2 == R_NilValue)
			error("SparseArray internal error in "
			      "_Arith_lv1_lv2():\n"
			      "    'lv1' and 'lv2' cannot both be NULL");
		if (opcode == SUB_OPCODE)
			return _unary_minus_leaf_vector(lv2, ans_Rtype);
		if (opcode == MULT_OPCODE)
			return mult0_leaf_vector(lv2, ans_Rtype,
						 offs_buf, vals_buf);
		error("SparseArray internal error in "
		      "_Arith_lv1_lv2():\n"
		      "    'op' must be \"-\" or \"*\" when 'lv1' is NULL");
	}
	if (lv2 == R_NilValue) {
		if (opcode == MULT_OPCODE)
			return mult0_leaf_vector(lv1, ans_Rtype,
						 offs_buf, vals_buf);
		error("SparseArray internal error in "
		      "_Arith_lv1_lv2():\n"
		      "    'op' must be \"*\" when 'lv2' is NULL");
	}
	lv1_len = _split_leaf_vector(lv1, &lv1_offs, &lv1_vals);
	lv2_len = _split_leaf_vector(lv2, &lv2_offs, &lv2_vals);
	offs1_p = INTEGER(lv1_offs);
	offs2_p = INTEGER(lv2_offs);
	buf_Rtype = REALSXP;
	ans_len = -1;
	if (TYPEOF(lv1_vals) == INTSXP) {
		if (TYPEOF(lv2_vals) == INTSXP) {
			buf_Rtype = INTSXP;
			ans_len = sparse_Arith_ints_ints(
				offs1_p, INTEGER(lv1_vals), lv1_len,
				offs2_p, INTEGER(lv2_vals), lv2_len,
				opcode, offs_buf, (int *) vals_buf, ovflow);
		} else if (TYPEOF(lv2_vals) == REALSXP) {
			ans_len = sparse_Arith_ints_doubles(
				offs1_p, INTEGER(lv1_vals), lv1_len,
				offs2_p, REAL(lv2_vals), lv2_len,
				opcode, offs_buf, (double *) vals_buf);
		}
	} else if (TYPEOF(lv1_vals) == REALSXP) {
		if (TYPEOF(lv2_vals) == INTSXP) {
			ans_len = sparse_Arith_doubles_ints(
				offs1_p, REAL(lv1_vals), lv1_len,
				offs2_p, INTEGER(lv2_vals), lv2_len,
				opcode, offs_buf, (double *) vals_buf);
		} else if (TYPEOF(lv2_vals) == REALSXP) {
			ans_len = sparse_Arith_doubles_doubles(
				offs1_p, REAL(lv1_vals), lv1_len,
				offs2_p, REAL(lv2_vals), lv2_len,
				opcode, offs_buf, (double *) vals_buf);
		}
	}
	if (ans_len == -1)
		error("_Arith_lv1_lv2() only supports input of "
		      "type \"integer\" or \"double\" at the moment");
	if (ans_Rtype != buf_Rtype)
		error("SparseArray internal error in "
		      "_Arith_lv1_lv2():\n"
		      "    ans_Rtype != buf_Rtype");
	return _make_leaf_vector_from_bufs(ans_Rtype,
					   offs_buf, vals_buf, ans_len);
}

