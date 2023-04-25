/****************************************************************************
 *                  'Compare' operations on "leaf vectors"                  *
 ****************************************************************************/
#include "leaf_vector_Compare.h"

#include "leaf_vector_utils.h"

int _get_Compare_opcode(SEXP op)
{
	const char *s;

	if (!IS_CHARACTER(op) || LENGTH(op) != 1)
		error("'op' must be a single string");
	op = STRING_ELT(op, 0);
	if (op == NA_STRING)
		error("'op' cannot be NA");
	s = CHAR(op);
	if (strcmp(s, "==") == 0)
		return EQ_OPCODE;
	if (strcmp(s, "!=") == 0)
		return NE_OPCODE;
	if (strcmp(s, "<=") == 0)
		return LE_OPCODE;
	if (strcmp(s, ">=") == 0)
		return GE_OPCODE;
	if (strcmp(s, "<") == 0)
		return LT_OPCODE;
	if (strcmp(s, ">") == 0)
		return GT_OPCODE;
	error("SparseArray internal error in "
	      "_get_Compare_opcode():\n"
	      "    invalid op: \"%s\"", s);
	return 0;  /* will never reach this */
}

static int flip_opcode(int opcode)
{
	switch (opcode) {
	    case NE_OPCODE: return opcode;
	    case LT_OPCODE: return GT_OPCODE;
	    case GT_OPCODE: return LT_OPCODE;
	}
	error("SparseArray internal error in flip_opcode():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static inline int Compare_int(int x, int y, int opcode)
{
	if (x == NA_INTEGER || y == NA_INTEGER)
		return NA_INTEGER;
	switch (opcode) {
	    case EQ_OPCODE: return x == y;
	    case NE_OPCODE: return x != y;
	    case LE_OPCODE: return x <= y;
	    case GE_OPCODE: return x >= y;
	    case LT_OPCODE: return x < y;
	    case GT_OPCODE: return x > y;
	}
	error("SparseArray internal error in Compare_int():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static inline int Compare_double(double x, double y, int opcode)
{
	if (ISNAN(x) || ISNAN(y))
		return NA_INTEGER;
	switch (opcode) {
	    case EQ_OPCODE: return x == y;
	    case NE_OPCODE: return x != y;
	    case LE_OPCODE: return x <= y;
	    case GE_OPCODE: return x >= y;
	    case LT_OPCODE: return x < y;
	    case GT_OPCODE: return x > y;
	}
	error("SparseArray internal error in Compare_double():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static int sparse_Compare_ints_1int(
		const int *offs1, const int *vals1, int n1, int v2,
		int opcode, int *offs_buf, int *vals_buf)
{
	int ans_len, k, v;

	for (ans_len = k = 0; k < n1; k++) {
		v = Compare_int(vals1[k], v2, opcode);
		if (v != 0) {
			offs_buf[ans_len] = offs1[k];
			vals_buf[ans_len] = v;
			ans_len++;
		}
	}
	return ans_len;
}

static int sparse_Compare_ints_1double(
		const int *offs1, const int *vals1, int n1, double v2,
		int opcode, int *offs_buf, int *vals_buf)
{
	int ans_len, k, v1, v;

	for (ans_len = k = 0; k < n1; k++) {
		v1 = vals1[k];
		if (v1 == NA_INTEGER) {
			v = NA_INTEGER;
		} else {
			v = Compare_double((double) v1, v2, opcode);
		}
		if (v != 0) {
			offs_buf[ans_len] = offs1[k];
			vals_buf[ans_len] = v;
			ans_len++;
		}
	}
	return ans_len;
}

static int sparse_Compare_ints_ints(
		const int *offs1, const int *vals1, int n1,
		const int *offs2, const int *vals2, int n2,
		int opcode, int *offs_buf, int *vals_buf)
{
	int ans_len, k1, k2, off, v1, v2, v;

	ans_len = k1 = k2 = 0;
	while (next_nzvals_int_int(offs1, vals1, n1,
				   offs2, vals2, n2,
				   &k1, &k2, &off, &v1, &v2))
	{
		v = Compare_int(v1, v2, opcode);
		if (v != 0) {
			offs_buf[ans_len] = off;
			vals_buf[ans_len] = v;
			ans_len++;
		}
	}
	return ans_len;
}

static int sparse_Compare_doubles_ints(
		const int *offs1, const double *vals1, int n1,
		const int *offs2, const int *vals2, int n2,
		int opcode, int *offs_buf, int *vals_buf)
{
	int ans_len, k1, k2, off, v2, v;
	double v1;

	ans_len = k1 = k2 = 0;
	while (next_nzvals_double_int(offs1, vals1, n1,
				      offs2, vals2, n2,
				      &k1, &k2, &off, &v1, &v2))
	{
		if (v2 == NA_INTEGER) {
			v = NA_INTEGER;
		} else {
			v = Compare_double(v1, (double) v2, opcode);
		}
		if (v != 0) {
			offs_buf[ans_len] = off;
			vals_buf[ans_len] = v;
			ans_len++;
		}
	}
	return ans_len;
}

static int sparse_Compare_doubles_1double(
		const int *offs1, const double *vals1, int n1, double v2,
		int opcode, int *offs_buf, int *vals_buf)
{
	int ans_len, k, v;

	for (ans_len = k = 0; k < n1; k++) {
		v = Compare_double(vals1[k], v2, opcode);
		if (v != 0) {
			offs_buf[ans_len] = offs1[k];
			vals_buf[ans_len] = v;
			ans_len++;
		}
	}
	return ans_len;
}

static int sparse_Compare_doubles_doubles(
		const int *offs1, const double *vals1, int n1,
		const int *offs2, const double *vals2, int n2,
		int opcode, int *offs_buf, int *vals_buf)
{
        int ans_len, k1, k2, off, v;
        double v1, v2;

	ans_len = k1 = k2 = 0;
	while (next_nzvals_double_double(offs1, vals1, n1,
					 offs2, vals2, n2,
					 &k1, &k2, &off, &v1, &v2))
	{
		v = Compare_double(v1, v2, opcode);
		if (v != 0) {
			offs_buf[ans_len] = off;
			vals_buf[ans_len] = v;
			ans_len++;
		}
	}
	return ans_len;
}

/* 'v2' is assumed to be != NA_INTEGER. This is NOT checked! */
static SEXP Compare_lv1_1int(SEXP lv1, int v2, int opcode,
			     int *offs_buf, int *vals_buf)
{
	int lv1_len, ans_len;
	SEXP lv1_offs, lv1_vals;
	const int *offs1_p;

	lv1_len = _split_leaf_vector(lv1, &lv1_offs, &lv1_vals);
	offs1_p = INTEGER(lv1_offs);
	ans_len = -1;
	switch (TYPEOF(lv1_vals)) {
	    case LGLSXP: case INTSXP:
		ans_len = sparse_Compare_ints_1int(
				offs1_p, INTEGER(lv1_vals), lv1_len, v2,
				opcode, offs_buf, vals_buf);
	    break;
	    case REALSXP:
		ans_len = sparse_Compare_doubles_1double(
				offs1_p, REAL(lv1_vals), lv1_len, (double) v2,
				opcode, offs_buf, vals_buf);
	    break;
	}
	if (ans_len == -1)
		error("Compare_lv1_1int() only supports input of "
		      "type \"logical\", \"integer\", or \"double\" "
		      "at the moment");
	return _new_leaf_vector_from_bufs(LGLSXP,
				offs_buf, vals_buf, ans_len);
}

static SEXP Compare_lv1_1double(SEXP lv1, double v2, int opcode,
				int *offs_buf, int *vals_buf)
{
	int lv1_len, ans_len;
	SEXP lv1_offs, lv1_vals;
	const int *offs1_p;

	lv1_len = _split_leaf_vector(lv1, &lv1_offs, &lv1_vals);
	offs1_p = INTEGER(lv1_offs);
	ans_len = -1;
	switch (TYPEOF(lv1_vals)) {
	    case LGLSXP: case INTSXP:
		ans_len = sparse_Compare_ints_1double(
				offs1_p, INTEGER(lv1_vals), lv1_len, v2,
				opcode, offs_buf, vals_buf);
	    break;
	    case REALSXP:
		ans_len = sparse_Compare_doubles_1double(
				offs1_p, REAL(lv1_vals), lv1_len, v2,
				opcode, offs_buf, vals_buf);
	    break;
	}
	if (ans_len == -1)
		error("Compare_lv1_1double() only supports input of "
		      "type \"logical\", \"integer\", or \"double\" "
		      "at the moment");
	return _new_leaf_vector_from_bufs(LGLSXP,
				offs_buf, vals_buf, ans_len);
}

/* 'v2' is assumed to be an atomic vector of length 1. This is NOT checked! */
SEXP _Compare_lv1_v2(SEXP lv1, SEXP v2, int opcode,
		     int *offs_buf, int *vals_buf)
{
	if (TYPEOF(v2) == LGLSXP || TYPEOF(v2) == INTSXP)
		return Compare_lv1_1int(lv1, INTEGER(v2)[0], opcode,
					offs_buf, vals_buf);
	if (TYPEOF(v2) == REALSXP)
		return Compare_lv1_1double(lv1, REAL(v2)[0], opcode,
					offs_buf, vals_buf);
	error("_Compare_lv1_v2() only supports input of "
	      "type \"logical\", \"integer\", or \"double\" "
	      "at the moment");
	return R_NilValue;  /* will never reach this */
}

/* Each of 'SVT1' and 'SVT2' must be a "leaf vector" or NULL. */
SEXP _Compare_lv1_lv2(SEXP lv1, SEXP lv2, int opcode,
		      int *offs_buf, int *vals_buf)
{
	int lv1_len, lv2_len, ans_len;
	SEXP lv1_offs, lv1_vals, lv2_offs, lv2_vals;
	const int *offs1_p, *offs2_p;

	if (lv1 == R_NilValue) {
		if (lv2 == R_NilValue)
			return R_NilValue;
		return Compare_lv1_1int(lv1, 0, flip_opcode(opcode),
					offs_buf, vals_buf);
	}
	if (lv2 == R_NilValue)
		return Compare_lv1_1int(lv1, 0, opcode,
					offs_buf, vals_buf);
	lv1_len = _split_leaf_vector(lv1, &lv1_offs, &lv1_vals);
	lv2_len = _split_leaf_vector(lv2, &lv2_offs, &lv2_vals);
	offs1_p = INTEGER(lv1_offs);
	offs2_p = INTEGER(lv2_offs);
	ans_len = -1;
	if (TYPEOF(lv1_vals) == INTSXP) {
		if (TYPEOF(lv2_vals) == INTSXP) {
			ans_len = sparse_Compare_ints_ints(
				offs1_p, INTEGER(lv1_vals), lv1_len,
				offs2_p, INTEGER(lv2_vals), lv2_len,
				opcode, offs_buf, vals_buf);
		} else if (TYPEOF(lv2_vals) == REALSXP) {
			ans_len = sparse_Compare_doubles_ints(
				offs2_p, REAL(lv2_vals), lv2_len,
				offs1_p, INTEGER(lv1_vals), lv1_len,
				flip_opcode(opcode), offs_buf, vals_buf);
		}
	} else if (TYPEOF(lv1_vals) == REALSXP) {
		if (TYPEOF(lv2_vals) == INTSXP) {
			ans_len = sparse_Compare_doubles_ints(
				offs1_p, REAL(lv1_vals), lv1_len,
				offs2_p, INTEGER(lv2_vals), lv2_len,
				opcode, offs_buf, vals_buf);
		} else if (TYPEOF(lv2_vals) == REALSXP) {
			ans_len = sparse_Compare_doubles_doubles(
				offs1_p, REAL(lv1_vals), lv1_len,
				offs2_p, REAL(lv2_vals), lv2_len,
				opcode, offs_buf, vals_buf);
		}
	}
	if (ans_len == -1)
		error("_Compare_lv1_lv2() only supports input of "
		      "type \"integer\" or \"double\" at the moment");
	return _new_leaf_vector_from_bufs(LGLSXP,
				offs_buf, vals_buf, ans_len);
}

