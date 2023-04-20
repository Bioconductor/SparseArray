/****************************************************************************
 *                     Ops operations on "leaf vectors"                     *
 ****************************************************************************/
#include "leaf_vector_Ops.h"

#include "leaf_vector_utils.h"

#include <string.h>  /* for memcpy() */
#include <limits.h>  /* for INT_MAX */


int _get_Arith_opcode(SEXP op, SEXPTYPE x_Rtype, SEXPTYPE y_Rtype)
{
	const char *s;

	if (!IS_CHARACTER(op) || LENGTH(op) != 1)
		error("'op' must be a single string");
	op = STRING_ELT(op, 0);
	if (op == NA_STRING)
		error("'op' cannot be NA");
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
	error("SparseArray internal error in "
	      "_get_Arith_opcode():\n"
	      "    invalid op: \"%s\"", s);
	return 0;
}

int _get_Compare_opcode(SEXP op, SEXPTYPE x_Rtype, SEXPTYPE y_Rtype)
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
	return 0;
}

int _get_Logic_opcode(SEXP op, SEXPTYPE x_Rtype, SEXPTYPE y_Rtype)
{
	const char *s;

	if (!IS_CHARACTER(op) || LENGTH(op) != 1)
		error("'op' must be a single string");
	op = STRING_ELT(op, 0);
	if (op == NA_STRING)
		error("'op' cannot be NA");
	s = CHAR(op);
	if (strcmp(s, "&") == 0)
		return AND_OPCODE;
	if (strcmp(s, "|") == 0)
		return OR_OPCODE;
	error("SparseArray internal error in "
	      "_get_Logic_opcode():\n"
	      "    invalid op: \"%s\"", s);
	return 0;
}

#define	ARGS_AND_BODY_OF_NEXT_NZVALS_FUNCTION(vals1_type, vals2_type)(	\
		const int *offs1, const vals1_type *vals1, int n1,	\
		const int *offs2, const vals2_type *vals2, int n2,	\
		int *k1, int *k2,					\
		int *off, vals1_type *v1, vals2_type *v2)		\
{									\
	int off1, off2;							\
									\
	if (*k1 < n1 && *k2 < n2) {					\
		off1 = offs1[*k1];					\
		off2 = offs2[*k2];					\
		if (off1 < off2) {					\
			*off = off1;					\
			*v1 = vals1[(*k1)++];				\
			*v2 = 0;					\
			return 1;					\
		}							\
		if (off1 > off2) {					\
			*off = off2;					\
			*v1 = 0;					\
			*v2 = vals2[(*k2)++];				\
			return 2;					\
		}							\
		*off = off1;						\
		*v1 = vals1[(*k1)++];					\
		*v2 = vals2[(*k2)++];					\
		return 3;						\
	}								\
	if (*k1 < n1) {							\
		*off = offs1[*k1];					\
		*v1 = vals1[(*k1)++];					\
		*v2 = 0;						\
		return 1;						\
	}								\
	if (*k2 < n2) {							\
		*off = offs2[*k2];					\
		*v1 = 0;						\
		*v2 = vals2[(*k2)++];					\
		return 2;						\
	}								\
	return 0;							\
}

static inline int next_nzvals_int
	ARGS_AND_BODY_OF_NEXT_NZVALS_FUNCTION(int, int)

static inline int next_nzvals_double
	ARGS_AND_BODY_OF_NEXT_NZVALS_FUNCTION(double, double)

static int Arith_sparse_ints(
		const int *offs1, const int *vals1, int n1,
		const int *offs2, const int *vals2, int n2,
		int opcode, int *offs_buf, int *vals_buf, int *ovflow)
{
	int ans_len, k1, k2, off, v1, v2;
	double v;

	ans_len = k1 = k2 = 0;
	while (next_nzvals_int(offs1, vals1, n1,
			       offs2, vals2, n2,
			       &k1, &k2, &off, &v1, &v2))
	{
		if (v1 == NA_INTEGER || v2 == NA_INTEGER) {
			offs_buf[ans_len] = off;
			vals_buf[ans_len] = NA_INTEGER;
			ans_len++;
			continue;
		}
		switch (opcode) {
		    case ADD_OPCODE:  v = (double) v1 + v2; break;
		    case SUB_OPCODE:  v = (double) v1 - v2; break;
		    case MULT_OPCODE: v = (double) v1 * v2; break;
		    default:
			error("SparseArray internal error in "
			      "Arith_sparse_ints():\n"
			      "    unsupported 'opcode'");
		}
		if (v != 0.0) {
			offs_buf[ans_len] = off;
			if (v <= INT_MAX && v > INT_MIN) {
				vals_buf[ans_len] = (int) v;
			} else {
				vals_buf[ans_len] = NA_INTEGER;
				*ovflow = 1;
			}
			ans_len++;
		}
	}
	return ans_len;
}

static int Arith_sparse_doubles(
		const int *offs1, const double *vals1, int n1,
		const int *offs2, const double *vals2, int n2,
		int opcode, int *offs_buf, double *vals_buf)
{
	int ans_len, k1, k2, off;
	double v1, v2, v;

	ans_len = k1 = k2 = 0;
	while (next_nzvals_double(offs1, vals1, n1,
				  offs2, vals2, n2,
				  &k1, &k2, &off, &v1, &v2))
	{
		switch (opcode) {
		    case ADD_OPCODE:  v = v1 + v2; break;
		    case SUB_OPCODE:  v = v1 - v2; break;
		    case MULT_OPCODE: v = v1 * v2; break;
		    default:
			error("SparseArray internal error in "
			      "Arith_sparse_doubles():\n"
			      "    unsupported 'opcode'");
		}
		if (v != 0.0) {
			offs_buf[ans_len] = off;
			vals_buf[ans_len] = v;
			ans_len++;
		}
	}
	return ans_len;
}

static SEXP mult0_leaf_vector(SEXP lv)
{
	error("mult0_leaf_vector() not ready yet!");
	return lv;
}

/* Empty leaf vectors are ILLEGAL so 'n' is NOT allowed to be 0. This is NOT
   checked! TODO: Should this go to leaf_vector_utils.c? */
static SEXP make_leaf_vector_from_offs_and_vals(SEXPTYPE Rtype,
		const int *offs, const void * vals, int n)
{
	SEXP ans, ans_offs, ans_vals;

	ans = PROTECT(_alloc_and_split_leaf_vector(n, Rtype,
						   &ans_offs, &ans_vals));
	memcpy(INTEGER(ans_offs), offs, sizeof(int) * n);
	switch (Rtype) {
	    case INTSXP:
		memcpy(INTEGER(ans_vals), (int *) vals, sizeof(int) * n);
		break;
	    case REALSXP:
		memcpy(REAL(ans_vals), (double *) vals, sizeof(double) * n);
		break;
	    default:
		UNPROTECT(1);
		error("SparseArray internal error in "
		      "make_leaf_vector_from_offs_and_vals():\n"
		      "    type \"%s\" is not supported",
		      type2char(Rtype));
	}
	UNPROTECT(1);
	return ans;
}

/* 'lv1' and 'lv2' must be "leaf vectors", with the exception that one of
   them can be NULL if 'opcode' is MULT_OPCODE. */
SEXP _Arith_leaf_vectors(SEXP lv1, SEXP lv2, int opcode,
			 int *offs_buf, void *vals_buf, int *ovflow)
{
	int lv1_len, lv2_len, ans_len;
	SEXP lv1_offs, lv1_vals, lv2_offs, lv2_vals;
	const int *offs1_p, *offs2_p;

	if (opcode == MULT_OPCODE) {
		if (lv1 == R_NilValue)
			return mult0_leaf_vector(lv2);
		if (lv2 == R_NilValue)
			return mult0_leaf_vector(lv1);
	}
	lv1_len = _split_leaf_vector(lv1, &lv1_offs, &lv1_vals);
	lv2_len = _split_leaf_vector(lv2, &lv2_offs, &lv2_vals);
	if (TYPEOF(lv1_vals) != TYPEOF(lv2_vals))
		error("_Arith_leaf_vectors() only supports "
		      "inputs of same type for now");
	offs1_p = INTEGER(lv1_offs);
	offs2_p = INTEGER(lv2_offs);
	if (TYPEOF(lv1_vals) == INTSXP) {
		ans_len = Arith_sparse_ints(
				offs1_p, INTEGER(lv1_vals), lv1_len,
				offs2_p, INTEGER(lv2_vals), lv2_len,
				opcode, offs_buf, (int *) vals_buf, ovflow);
	} else if (TYPEOF(lv1_vals) == REALSXP) {
		ans_len = Arith_sparse_doubles(
				offs1_p, REAL(lv1_vals), lv1_len,
				offs2_p, REAL(lv2_vals), lv2_len,
				opcode, offs_buf, (double *) vals_buf);
	} else {
		error("_Arith_leaf_vectors() only supports "
		      "input of type \"int\" or \"double\" for now");
	}
	if (ans_len == 0)
		return R_NilValue;
	return make_leaf_vector_from_offs_and_vals(TYPEOF(lv1_vals),
				offs_buf, vals_buf, ans_len);
}

