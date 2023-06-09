/****************************************************************************
 *                   'Logic' operations on "leaf vectors"                   *
 ****************************************************************************/
#include "leaf_vector_Logic.h"

#include "leaf_vector_utils.h"

int _get_Logic_opcode(SEXP op)
{
	const char *s;

	if (!IS_CHARACTER(op) || LENGTH(op) != 1)
		error("SparseArray internal error in _get_Logic_opcode():\n"
		      "    'op' must be a single string");
	op = STRING_ELT(op, 0);
	if (op == NA_STRING)
		error("SparseArray internal error in _get_Logic_opcode():\n"
		      "    'op' cannot be NA");
	s = CHAR(op);
	if (strcmp(s, "&") == 0)
		return AND_OPCODE;
	if (strcmp(s, "|") == 0)
		return OR_OPCODE;
	error("SparseArray internal error in _get_Logic_opcode():\n"
	      "    invalid op: \"%s\"", s);
	return 0;  /* will never reach this */
}

static inline int Logic_int_int(int x, int y, int opcode)
{
	switch (opcode) {
	    case AND_OPCODE:
		if (x == 0 || y == 0)
			return 0;
		if (x == NA_INTEGER || y == NA_INTEGER)
			return NA_INTEGER;
		return 1;
	    case OR_OPCODE:
		if (x == 1 || y == 1)
			return 1;
		if (x == NA_INTEGER || y == NA_INTEGER)
			return NA_INTEGER;
		return 0;
	}
	error("SparseArray internal error in Logic_int_int():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static int sparse_Logic_ints_ints(
		const int *offs1, const int *vals1, int n1,
		const int *offs2, const int *vals2, int n2,
		int opcode, int *offs_buf, int *vals_buf)
{
	int ans_len, k1, k2, off, x, y, v;

	ans_len = k1 = k2 = 0;
	while (next_nzval_int_int(offs1, vals1, n1,
				  offs2, vals2, n2,
				  &k1, &k2, &off, &x, &y))
	{
		v = Logic_int_int(x, y, opcode);
		if (v != 0) {
			offs_buf[ans_len] = off;
			vals_buf[ans_len] = v;
			ans_len++;
		}
	}
	return ans_len;
}

/* Each of 'lv1' and 'lv2' must be a "leaf vector" or NULL. */
SEXP _Logic_lv1_lv2(SEXP lv1, SEXP lv2, int opcode,
		    int *offs_buf, int *vals_buf)
{
	int lv1_len, lv2_len, ans_len;
	SEXP lv1_offs, lv1_vals, lv2_offs, lv2_vals;

	if (lv1 == R_NilValue || lv2 == R_NilValue) {
		if (opcode == AND_OPCODE)
			return R_NilValue;
		return lv1 == R_NilValue ? lv2 : lv1;
	}
	lv1_len = _split_leaf_vector(lv1, &lv1_offs, &lv1_vals);
	lv2_len = _split_leaf_vector(lv2, &lv2_offs, &lv2_vals);
	ans_len = sparse_Logic_ints_ints(
			INTEGER(lv1_offs), INTEGER(lv1_vals), lv1_len,
			INTEGER(lv2_offs), INTEGER(lv2_vals), lv2_len,
			opcode, offs_buf, vals_buf);
	return _make_leaf_vector_from_bufs(LGLSXP,
					   offs_buf, vals_buf, ans_len);
}

