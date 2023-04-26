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

#define ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN1(Ltype, Rtype)(		\
		const int *offs1, const Ltype *vals1, int n1,		\
		Rtype v2,						\
		int opcode, int *offs_buf, int *vals_buf)		\
{									\
	int ans_len, k, v;						\
									\
	for (ans_len = k = 0; k < n1; k++) {				\
		v = Compare_ ## Ltype ## _ ## Rtype			\
					(vals1[k], v2, opcode);		\
		if (v != 0) {						\
			offs_buf[ans_len] = offs1[k];			\
			vals_buf[ans_len] = v;				\
			ans_len++;					\
		}							\
	}								\
	return ans_len;							\
}

#define ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN2(Ltype, Rtype)(		\
		const int *offs1, const Ltype *vals1, int n1,		\
		const int *offs2, const Rtype *vals2, int n2,		\
		int opcode, int *offs_buf, int *vals_buf)		\
{									\
	int ans_len, k1, k2, off, v;					\
	Ltype v1;							\
	Rtype v2;							\
									\
	ans_len = k1 = k2 = 0;						\
	while (next_nzvals_ ## Ltype ## _ ## Rtype			\
					(offs1, vals1, n1,		\
					 offs2, vals2, n2,		\
					 &k1, &k2, &off, &v1, &v2))	\
	{								\
		v = Compare_ ## Ltype ## _ ## Rtype			\
					(v1, v2, opcode);		\
		if (v != 0) {						\
			offs_buf[ans_len] = off;			\
			vals_buf[ans_len] = v;				\
			ans_len++;					\
		}							\
	}								\
	return ans_len;							\
}


/****************************************************************************
 * sparse_Compare_Rbytes_*() functions
 */

static inline int Compare_Rbyte_Rbyte(Rbyte x, Rbyte y, int opcode)
{
	switch (opcode) {
	    case EQ_OPCODE: return x == y;
	    case NE_OPCODE: return x != y;
	    case LE_OPCODE: return x <= y;
	    case GE_OPCODE: return x >= y;
	    case LT_OPCODE: return x < y;
	    case GT_OPCODE: return x > y;
	}
	error("SparseArray internal error in Compare_Rbyte_Rbyte():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static int sparse_Compare_Rbytes_Rbyte
	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN1(Rbyte, Rbyte)
static int sparse_Compare_Rbytes_Rbytes
	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN2(Rbyte, Rbyte)

static inline int Compare_Rbyte_int(Rbyte x, int y, int opcode)
{
	int x1;

	x1 = (int) x;
	switch (opcode) {
	    case EQ_OPCODE: return x1 == y;
	    case NE_OPCODE: return x1 != y;
	    case LE_OPCODE: return x1 <= y;
	    case GE_OPCODE: return x1 >= y;
	    case LT_OPCODE: return x1 < y;
	    case GT_OPCODE: return x1 > y;
	}
	error("SparseArray internal error in Compare_Rbyte_int():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static int sparse_Compare_Rbytes_int
	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN1(Rbyte, int)
static int sparse_Compare_Rbytes_ints
	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN2(Rbyte, int)

static inline int Compare_Rbyte_double(Rbyte x, double y, int opcode)
{
	double x1;

	x1 = (double) x;
	switch (opcode) {
	    case EQ_OPCODE: return x1 == y;
	    case NE_OPCODE: return x1 != y;
	    case LE_OPCODE: return x1 <= y;
	    case GE_OPCODE: return x1 >= y;
	    case LT_OPCODE: return x1 < y;
	    case GT_OPCODE: return x1 > y;
	}
	error("SparseArray internal error in Compare_Rbyte_double():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static int sparse_Compare_Rbytes_double
	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN1(Rbyte, double)
static int sparse_Compare_Rbytes_doubles
	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN2(Rbyte, double)

static inline int Compare_Rbyte_Rcomplex(Rbyte x, Rcomplex y, int opcode)
{
	double x1;

	x1 = (double) x;
	switch (opcode) {
	    case EQ_OPCODE: return x1 == y.r && 0.0 == y.i;
	    case NE_OPCODE: return x1 != y.r || 0.0 != y.i;
	}
	error("SparseArray internal error in Compare_Rbyte_Rcomplex():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static int sparse_Compare_Rbytes_Rcomplex
	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN1(Rbyte, Rcomplex)
//static int sparse_Compare_Rbytes_Rcomplexes
//	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN2(Rbyte, Rcomplex)


/****************************************************************************
 * sparse_Compare_ints_*() functions
 */

static inline int Compare_int_int(int x, int y, int opcode)
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
	error("SparseArray internal error in Compare_int_int():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static int sparse_Compare_ints_int
	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN1(int, int)
static int sparse_Compare_ints_ints
	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN2(int, int)

static inline int Compare_int_double(int x, double y, int opcode)
{
	double x1;

	if (x == NA_INTEGER || ISNAN(y))
		return NA_INTEGER;
	x1 = (double) x;
	switch (opcode) {
	    case EQ_OPCODE: return x1 == y;
	    case NE_OPCODE: return x1 != y;
	    case LE_OPCODE: return x1 <= y;
	    case GE_OPCODE: return x1 >= y;
	    case LT_OPCODE: return x1 < y;
	    case GT_OPCODE: return x1 > y;
	}
	error("SparseArray internal error in Compare_int_double():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static int sparse_Compare_ints_double
	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN1(int, double)
static int sparse_Compare_ints_doubles
	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN2(int, double)

static inline int Compare_int_Rcomplex(int x, Rcomplex y, int opcode)
{
	double x1;

	if (x == NA_INTEGER || ISNAN(y.r) || ISNAN(y.i))
		return NA_INTEGER;
	x1 = (double) x;
	switch (opcode) {
	    case EQ_OPCODE: return x1 == y.r && 0.0 == y.i;
	    case NE_OPCODE: return x1 != y.r || 0.0 != y.i;
	}
	error("SparseArray internal error in Compare_int_Rcomplex():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static int sparse_Compare_ints_Rcomplex
	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN1(int, Rcomplex)
//static int sparse_Compare_ints_Rcomplexes
//	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN2(int, Rcomplex)


/****************************************************************************
 * sparse_Compare_doubles_*() functions
 */

static inline int Compare_double_double(double x, double y, int opcode)
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
	error("SparseArray internal error in Compare_double_double():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static int sparse_Compare_doubles_double
	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN1(double, double)
static int sparse_Compare_doubles_doubles
	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN2(double, double)

static inline int Compare_double_Rcomplex(double x, Rcomplex y, int opcode)
{
	if (ISNAN(x) || ISNAN(y.r) || ISNAN(y.i))
		return NA_INTEGER;
	switch (opcode) {
	    case EQ_OPCODE: return x == y.r && 0.0 == y.i;
	    case NE_OPCODE: return x != y.r || 0.0 != y.i;
	}
	error("SparseArray internal error in Compare_double_Rcomplex():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static int sparse_Compare_doubles_Rcomplex
	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN1(double, Rcomplex)
//static int sparse_Compare_doubles_Rcomplexes
//	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN2(double, Rcomplex)


/****************************************************************************
 * sparse_Compare_Rcomplexes_*() functions
 */

static inline int Compare_Rcomplex_Rcomplex(Rcomplex x, Rcomplex y, int opcode)
{
	if (ISNAN(x.r) || ISNAN(x.i) || ISNAN(y.r) || ISNAN(y.i))
		return NA_INTEGER;
	switch (opcode) {
	    case EQ_OPCODE: return x.r == y.r && x.i == y.i;
	    case NE_OPCODE: return x.r != y.r || x.i != y.i;
	}
	error("SparseArray internal error in Compare_Rcomplex_Rcomplex():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static int sparse_Compare_Rcomplexes_Rcomplex
	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN1(Rcomplex, Rcomplex)
//static int sparse_Compare_Rcomplexes_Rcomplexes
//	ARGS_AND_BODY_OF_SPARSE_COMPARE_FUN2(Rcomplex, Rcomplex)


/****************************************************************************
 * Compare_lv1_1raw()
 * Compare_lv1_1int()
 * Compare_lv1_1double()
 * _Compare_lv1_v2()
 * _Compare_lv1_lv2()
 */

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
		ans_len = sparse_Compare_ints_int(
				offs1_p, INTEGER(lv1_vals), lv1_len, v2,
				opcode, offs_buf, vals_buf);
	    break;
	    case REALSXP:
		ans_len = sparse_Compare_doubles_double(
				offs1_p, REAL(lv1_vals), lv1_len, (double) v2,
				opcode, offs_buf, vals_buf);
	    break;
	    case RAWSXP:
		ans_len = sparse_Compare_Rbytes_double(
				offs1_p, RAW(lv1_vals), lv1_len, (double) v2,
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
		ans_len = sparse_Compare_ints_double(
				offs1_p, INTEGER(lv1_vals), lv1_len, v2,
				opcode, offs_buf, vals_buf);
	    break;
	    case REALSXP:
		ans_len = sparse_Compare_doubles_double(
				offs1_p, REAL(lv1_vals), lv1_len, v2,
				opcode, offs_buf, vals_buf);
	    break;
	    case RAWSXP:
		ans_len = sparse_Compare_Rbytes_double(
				offs1_p, RAW(lv1_vals), lv1_len, v2,
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
	//if (TYPEOF(v2) == RAWSXP)
	//	return Compare_lv1_1raw(lv1, RAW(v2)[0], opcode,
	//				offs_buf, vals_buf);
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
			ans_len = sparse_Compare_ints_doubles(
				offs1_p, INTEGER(lv1_vals), lv1_len,
				offs2_p, REAL(lv2_vals), lv2_len,
				opcode, offs_buf, vals_buf);
		}
	} else if (TYPEOF(lv1_vals) == REALSXP) {
		if (TYPEOF(lv2_vals) == INTSXP) {
			ans_len = sparse_Compare_ints_doubles(
				offs2_p, INTEGER(lv2_vals), lv2_len,
				offs1_p, REAL(lv1_vals), lv1_len,
				flip_opcode(opcode), offs_buf, vals_buf);
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

