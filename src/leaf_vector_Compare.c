/****************************************************************************
 *                  'Compare' operations on "leaf vectors"                  *
 ****************************************************************************/
#include "leaf_vector_Compare.h"

#include "leaf_vector_utils.h"

int _get_Compare_opcode(SEXP op)
{
	const char *s;

	if (!IS_CHARACTER(op) || LENGTH(op) != 1)
		error("SparseArray internal error in _get_Compare_opcode():\n"
		      "    'op' must be a single string");
	op = STRING_ELT(op, 0);
	if (op == NA_STRING)
		error("SparseArray internal error in _get_Compare_opcode():\n"
		      "    'op' cannot be NA");
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
	error("SparseArray internal error in _get_Compare_opcode():\n"
	      "    invalid op: \"%s\"", s);
	return 0;  /* will never reach this */
}

static inline int flip_opcode(int opcode)
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

#define FUNDEF_sparse_Compare1(Ltype, Rtype)				\
	(const int *offs1, const Ltype *vals1, int n1, Rtype y,		\
	 int opcode, int *offs_buf, int *vals_buf)			\
{									\
	int ans_len, k, v;						\
									\
	for (ans_len = k = 0; k < n1; k++) {				\
		v = Compare_ ## Ltype ## _ ## Rtype			\
					(vals1[k], y, opcode);		\
		if (v != 0) {						\
			offs_buf[ans_len] = offs1[k];			\
			vals_buf[ans_len] = v;				\
			ans_len++;					\
		}							\
	}								\
	return ans_len;							\
}

#define FUNDEF_sparse_Compare2(Ltype, Rtype)				\
	(const int *offs1, const Ltype *vals1, int n1,			\
	 const int *offs2, const Rtype *vals2, int n2,			\
	 int opcode, int *offs_buf, int *vals_buf)			\
{									\
	int ans_len, k1, k2, off, v;					\
	Ltype x;							\
	Rtype y;							\
									\
	ans_len = k1 = k2 = 0;						\
	while (next_nzvals_ ## Ltype ## _ ## Rtype			\
					(offs1, vals1, n1,		\
					 offs2, vals2, n2,		\
					 &k1, &k2, &off, &x, &y))	\
	{								\
		v = Compare_ ## Ltype ## _ ## Rtype			\
					(x, y, opcode);			\
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
	FUNDEF_sparse_Compare1(Rbyte, Rbyte)
static int sparse_Compare_Rbytes_Rbytes
	FUNDEF_sparse_Compare2(Rbyte, Rbyte)

/* WARNING: Only valid to use on the int values of an integer vector (INTSXP).
   Should NOT be used on the int values of a logical vector (LGLSEXP), on
   which Compare_Rbyte_int() will generally produce wrong results!
   For example if Rbyte value 'x' is >= 2 and int value 'y' is 1 (TRUE),
   then == and <= will both return 0 when they are expected to return 1. */
static inline int Compare_Rbyte_int(Rbyte x, int y, int opcode)
{
	int x1;

	if (y == NA_INTEGER)
		return NA_INTEGER;
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
	FUNDEF_sparse_Compare1(Rbyte, int)
static int sparse_Compare_Rbytes_ints
	FUNDEF_sparse_Compare2(Rbyte, int)

static inline int Compare_Rbyte_double(Rbyte x, double y, int opcode)
{
	double x1;

	if (ISNAN(y))
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
	error("SparseArray internal error in Compare_Rbyte_double():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static int sparse_Compare_Rbytes_double
	FUNDEF_sparse_Compare1(Rbyte, double)
static int sparse_Compare_Rbytes_doubles
	FUNDEF_sparse_Compare2(Rbyte, double)

static inline int Compare_Rbyte_Rcomplex(Rbyte x, Rcomplex y, int opcode)
{
	double x1;

	if (ISNAN(y.r) || ISNAN(y.i))
		return NA_INTEGER;
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
	FUNDEF_sparse_Compare1(Rbyte, Rcomplex)
static int sparse_Compare_Rbytes_Rcomplexes
	FUNDEF_sparse_Compare2(Rbyte, Rcomplex)

/* 'v2' is assumed to be an atomic vector of length 1. This is NOT checked!
   Also its type is expected to be "raw" or bigger. */
static int sparse_Compare_Rbytes_v2(
		const int *offs1, const Rbyte *vals1, int n1, SEXP v2,
		int opcode, int *offs_buf, int *vals_buf)
{
	SEXPTYPE Rtype2;

	Rtype2 = TYPEOF(v2);
	switch (Rtype2) {
	    case RAWSXP:
		return sparse_Compare_Rbytes_Rbyte(
				offs1, vals1, n1, RAW(v2)[0],
				opcode, offs_buf, vals_buf);
	    case LGLSXP: case INTSXP:
		return sparse_Compare_Rbytes_int(
				offs1, vals1, n1, INTEGER(v2)[0],
				opcode, offs_buf, vals_buf);
	    case REALSXP:
		return sparse_Compare_Rbytes_double(
				offs1, vals1, n1, REAL(v2)[0],
				opcode, offs_buf, vals_buf);
	    case CPLXSXP:
		return sparse_Compare_Rbytes_Rcomplex(
				offs1, vals1, n1, COMPLEX(v2)[0],
				opcode, offs_buf, vals_buf);
	}
	error("SparseArray internal error in "
	      "sparse_Compare_Rbytes_v2():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
	return 0;  /* will never reach this */
}


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
	FUNDEF_sparse_Compare1(int, int)
static int sparse_Compare_ints_ints
	FUNDEF_sparse_Compare2(int, int)

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
	FUNDEF_sparse_Compare1(int, double)
static int sparse_Compare_ints_doubles
	FUNDEF_sparse_Compare2(int, double)

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
	FUNDEF_sparse_Compare1(int, Rcomplex)
static int sparse_Compare_ints_Rcomplexes
	FUNDEF_sparse_Compare2(int, Rcomplex)

/* 'v2' is assumed to be an atomic vector of length 1. This is NOT checked!
   Also its type is expected to be "integer" or bigger. */
static int sparse_Compare_ints_v2(
		const int *offs1, const int *vals1, int n1, SEXP v2,
		int opcode, int *offs_buf, int *vals_buf)
{
	SEXPTYPE Rtype2;

	Rtype2 = TYPEOF(v2);
	switch (Rtype2) {
	    case LGLSXP: case INTSXP:
		return sparse_Compare_ints_int(
				offs1, vals1, n1, INTEGER(v2)[0],
				opcode, offs_buf, vals_buf);
	    case REALSXP:
		return sparse_Compare_ints_double(
				offs1, vals1, n1, REAL(v2)[0],
				opcode, offs_buf, vals_buf);
	    case CPLXSXP:
		return sparse_Compare_ints_Rcomplex(
				offs1, vals1, n1, COMPLEX(v2)[0],
				opcode, offs_buf, vals_buf);
	}
	error("SparseArray internal error in "
	      "sparse_Compare_ints_v2():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
	return 0;  /* will never reach this */
}


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
	FUNDEF_sparse_Compare1(double, double)
static int sparse_Compare_doubles_doubles
	FUNDEF_sparse_Compare2(double, double)

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
	FUNDEF_sparse_Compare1(double, Rcomplex)
static int sparse_Compare_doubles_Rcomplexes
	FUNDEF_sparse_Compare2(double, Rcomplex)

/* 'v2' is assumed to be an atomic vector of length 1. This is NOT checked!
   Also its type is expected to be "double" or bigger. */
static int sparse_Compare_doubles_v2(
		const int *offs1, const double *vals1, int n1, SEXP v2,
		int opcode, int *offs_buf, int *vals_buf)
{
	SEXPTYPE Rtype2;

	Rtype2 = TYPEOF(v2);
	switch (Rtype2) {
	    case REALSXP:
		return sparse_Compare_doubles_double(
				offs1, vals1, n1, REAL(v2)[0],
				opcode, offs_buf, vals_buf);
	    case CPLXSXP:
		return sparse_Compare_doubles_Rcomplex(
				offs1, vals1, n1, COMPLEX(v2)[0],
				opcode, offs_buf, vals_buf);
	}
	error("SparseArray internal error in "
	      "sparse_Compare_doubles_v2():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
	return 0;  /* will never reach this */
}


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
	FUNDEF_sparse_Compare1(Rcomplex, Rcomplex)
static int sparse_Compare_Rcomplexes_Rcomplexes
	FUNDEF_sparse_Compare2(Rcomplex, Rcomplex)

/* 'v2' is assumed to be an atomic vector of length 1. This is NOT checked!
   Also its type is expected to be "complex". */
static int sparse_Compare_Rcomplexes_v2(
		const int *offs1, const Rcomplex *vals1, int n1, SEXP v2,
		int opcode, int *offs_buf, int *vals_buf)
{
	SEXPTYPE Rtype2;

	Rtype2 = TYPEOF(v2);
	switch (Rtype2) {
	    case CPLXSXP:
		return sparse_Compare_Rcomplexes_Rcomplex(
				offs1, vals1, n1, COMPLEX(v2)[0],
				opcode, offs_buf, vals_buf);
	}
	error("SparseArray internal error in "
	      "sparse_Compare_Rcomplexes_v2():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
	return 0;  /* will never reach this */
}


/****************************************************************************
 * sparse_Compare_Rbytes_lv2()
 * sparse_Compare_ints_lv2()
 * sparse_Compare_doubles_lv2()
 * sparse_Compare_Rcomplexes_lv2()
 */

static int sparse_Compare_Rbytes_lv2(
		const int *offs1, const Rbyte *vals1, int n1, SEXP lv2,
		int opcode, int *offs_buf, int *vals_buf)
{
	int lv2_len;
	SEXP lv2_offs, lv2_vals;
	const int *offs2;
	SEXPTYPE Rtype2;

	lv2_len = _split_leaf_vector(lv2, &lv2_offs, &lv2_vals);
	offs2 = INTEGER(lv2_offs);
	Rtype2 = TYPEOF(lv2_vals);
	switch (Rtype2) {
	    case RAWSXP:
		return sparse_Compare_Rbytes_Rbytes(
				offs1, vals1, n1,
				offs2, RAW(lv2_vals), lv2_len,
				opcode, offs_buf, vals_buf);
	    case LGLSXP: case INTSXP:
		return sparse_Compare_Rbytes_ints(
				offs1, vals1, n1,
				offs2, INTEGER(lv2_vals), lv2_len,
				opcode, offs_buf, vals_buf);
	    case REALSXP:
		return sparse_Compare_Rbytes_doubles(
				offs1, vals1, n1,
				offs2, REAL(lv2_vals), lv2_len,
				opcode, offs_buf, vals_buf);
	    case CPLXSXP:
		return sparse_Compare_Rbytes_Rcomplexes(
				offs1, vals1, n1,
				offs2, COMPLEX(lv2_vals), lv2_len,
				opcode, offs_buf, vals_buf);
	}
	error("SparseArray internal error in "
	      "sparse_Compare_Rbytes_lv2():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
	return 0;  /* will never reach this */
}

static int sparse_Compare_ints_lv2(
		const int *offs1, const int *vals1, int n1, SEXP lv2,
		int opcode, int *offs_buf, int *vals_buf)
{
	int lv2_len;
	SEXP lv2_offs, lv2_vals;
	const int *offs2;
	SEXPTYPE Rtype2;

	lv2_len = _split_leaf_vector(lv2, &lv2_offs, &lv2_vals);
	offs2 = INTEGER(lv2_offs);
	Rtype2 = TYPEOF(lv2_vals);
	switch (Rtype2) {
	    case RAWSXP:
		return sparse_Compare_Rbytes_ints(
				offs2, RAW(lv2_vals), lv2_len,
				offs1, vals1, n1,
				flip_opcode(opcode), offs_buf, vals_buf);
	    case LGLSXP: case INTSXP:
		return sparse_Compare_ints_ints(
				offs1, vals1, n1,
				offs2, INTEGER(lv2_vals), lv2_len,
				opcode, offs_buf, vals_buf);
	    case REALSXP:
		return sparse_Compare_ints_doubles(
				offs1, vals1, n1,
				offs2, REAL(lv2_vals), lv2_len,
				opcode, offs_buf, vals_buf);
	    case CPLXSXP:
		return sparse_Compare_ints_Rcomplexes(
				offs1, vals1, n1,
				offs2, COMPLEX(lv2_vals), lv2_len,
				opcode, offs_buf, vals_buf);
	}
	error("SparseArray internal error in "
	      "sparse_Compare_ints_lv2():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
	return 0;  /* will never reach this */
}

static int sparse_Compare_doubles_lv2(
		const int *offs1, const double *vals1, int n1, SEXP lv2,
		int opcode, int *offs_buf, int *vals_buf)
{
	int lv2_len;
	SEXP lv2_offs, lv2_vals;
	const int *offs2;
	SEXPTYPE Rtype2;

	lv2_len = _split_leaf_vector(lv2, &lv2_offs, &lv2_vals);
	offs2 = INTEGER(lv2_offs);
	Rtype2 = TYPEOF(lv2_vals);
	switch (Rtype2) {
	    case RAWSXP:
		return sparse_Compare_Rbytes_doubles(
				offs2, RAW(lv2_vals), lv2_len,
				offs1, vals1, n1,
				flip_opcode(opcode), offs_buf, vals_buf);
	    case LGLSXP: case INTSXP:
		return sparse_Compare_ints_doubles(
				offs2, INTEGER(lv2_vals), lv2_len,
				offs1, vals1, n1,
				flip_opcode(opcode), offs_buf, vals_buf);
	    case REALSXP:
		return sparse_Compare_doubles_doubles(
				offs1, vals1, n1,
				offs2, REAL(lv2_vals), lv2_len,
				opcode, offs_buf, vals_buf);
	    case CPLXSXP:
		return sparse_Compare_doubles_Rcomplexes(
				offs1, vals1, n1,
				offs2, COMPLEX(lv2_vals), lv2_len,
				opcode, offs_buf, vals_buf);
	}
	error("SparseArray internal error in "
	      "sparse_Compare_doubles_lv2():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
	return 0;  /* will never reach this */
}

static int sparse_Compare_Rcomplexes_lv2(
		const int *offs1, const Rcomplex *vals1, int n1, SEXP lv2,
		int opcode, int *offs_buf, int *vals_buf)
{
	int lv2_len;
	SEXP lv2_offs, lv2_vals;
	const int *offs2;
	SEXPTYPE Rtype2;

	lv2_len = _split_leaf_vector(lv2, &lv2_offs, &lv2_vals);
	offs2 = INTEGER(lv2_offs);
	Rtype2 = TYPEOF(lv2_vals);
	switch (Rtype2) {
	    case RAWSXP:
		return sparse_Compare_Rbytes_Rcomplexes(
				offs2, RAW(lv2_vals), lv2_len,
				offs1, vals1, n1,
				flip_opcode(opcode), offs_buf, vals_buf);
	    case LGLSXP: case INTSXP:
		return sparse_Compare_ints_Rcomplexes(
				offs2, INTEGER(lv2_vals), lv2_len,
				offs1, vals1, n1,
				flip_opcode(opcode), offs_buf, vals_buf);
	    case REALSXP:
		return sparse_Compare_doubles_Rcomplexes(
				offs2, REAL(lv2_vals), lv2_len,
				offs1, vals1, n1,
				flip_opcode(opcode), offs_buf, vals_buf);
	    case CPLXSXP:
		return sparse_Compare_Rcomplexes_Rcomplexes(
				offs1, vals1, n1,
				offs2, COMPLEX(lv2_vals), lv2_len,
				opcode, offs_buf, vals_buf);
	}
	error("SparseArray internal error in "
	      "sparse_Compare_Rcomplexes_lv2():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
	return 0;  /* will never reach this */
}


/****************************************************************************
 * Compare_lv1_zero()
 * _Compare_lv1_v2()
 * _Compare_lv1_lv2()
 */

static SEXP Compare_lv1_zero(SEXP lv1, int opcode,
			     int *offs_buf, int *vals_buf)
{
	int lv1_len, ans_len;
	SEXP lv1_offs, lv1_vals;
	const int *offs1;

	lv1_len = _split_leaf_vector(lv1, &lv1_offs, &lv1_vals);
	offs1 = INTEGER(lv1_offs);
	switch (TYPEOF(lv1_vals)) {
	    case RAWSXP:
		ans_len = sparse_Compare_Rbytes_Rbyte(
				offs1, RAW(lv1_vals), lv1_len, Rbyte0,
				opcode, offs_buf, vals_buf);
	    break;
	    case LGLSXP: case INTSXP:
		ans_len = sparse_Compare_ints_int(
				offs1, INTEGER(lv1_vals), lv1_len, int0,
				opcode, offs_buf, vals_buf);
	    break;
	    case REALSXP:
		ans_len = sparse_Compare_doubles_double(
				offs1, REAL(lv1_vals), lv1_len, double0,
				opcode, offs_buf, vals_buf);
	    break;
	    case CPLXSXP:
		ans_len = sparse_Compare_Rcomplexes_Rcomplex(
				offs1, COMPLEX(lv1_vals), lv1_len, Rcomplex0,
				opcode, offs_buf, vals_buf);
	    break;
	    default:
		error("SparseArray internal error in "
		      "Compare_lv1_zero():\n"
		      "    unsupported 'TYPEOF(lv1_vals)': \"%s\"",
		      type2char(TYPEOF(lv1_vals)));
	}
	return _make_leaf_vector_from_bufs(LGLSXP,
					   offs_buf, vals_buf, ans_len);
}

/* 'v2' is assumed to be an atomic vector of length 1. This is NOT checked!
   Also its type is expected to be the same size as lv1's type or bigger. */
SEXP _Compare_lv1_v2(SEXP lv1, SEXP v2, int opcode,
		     int *offs_buf, int *vals_buf)
{
	int lv1_len, ans_len;
	SEXP lv1_offs, lv1_vals;
	const int *offs1;

	lv1_len = _split_leaf_vector(lv1, &lv1_offs, &lv1_vals);
	offs1 = INTEGER(lv1_offs);
	switch (TYPEOF(lv1_vals)) {
	    case RAWSXP:
		ans_len = sparse_Compare_Rbytes_v2(
				offs1, RAW(lv1_vals), lv1_len, v2,
				opcode, offs_buf, vals_buf);
	    break;
	    case LGLSXP: case INTSXP:
		ans_len = sparse_Compare_ints_v2(
				offs1, INTEGER(lv1_vals), lv1_len, v2,
				opcode, offs_buf, vals_buf);
	    break;
	    case REALSXP:
		ans_len = sparse_Compare_doubles_v2(
				offs1, REAL(lv1_vals), lv1_len, v2,
				opcode, offs_buf, vals_buf);
	    break;
	    case CPLXSXP:
		ans_len = sparse_Compare_Rcomplexes_v2(
				offs1, COMPLEX(lv1_vals), lv1_len, v2,
				opcode, offs_buf, vals_buf);
	    break;
	    default:
		error("SparseArray internal error in "
		      "_Compare_lv1_v2():\n"
		      "    unsupported 'TYPEOF(lv1_vals)': \"%s\"",
		      type2char(TYPEOF(lv1_vals)));
	}
	return _make_leaf_vector_from_bufs(LGLSXP,
					   offs_buf, vals_buf, ans_len);
}

/* Each of 'lv1' and 'lv2' must be a "leaf vector" or NULL. */
SEXP _Compare_lv1_lv2(SEXP lv1, SEXP lv2, int opcode,
		      int *offs_buf, int *vals_buf)
{
	int lv1_len, ans_len;
	SEXP lv1_offs, lv1_vals;
	const int *offs1;
	SEXPTYPE Rtype1;

	if (lv1 == R_NilValue) {
		if (lv2 == R_NilValue)
			return R_NilValue;
		return Compare_lv1_zero(lv2, flip_opcode(opcode),
					offs_buf, vals_buf);
	}
	if (lv2 == R_NilValue)
		return Compare_lv1_zero(lv1, opcode,
					offs_buf, vals_buf);
	lv1_len = _split_leaf_vector(lv1, &lv1_offs, &lv1_vals);
	offs1 = INTEGER(lv1_offs);
	Rtype1 = TYPEOF(lv1_vals);
	switch (Rtype1) {
	    case RAWSXP:
		ans_len = sparse_Compare_Rbytes_lv2(
				offs1, RAW(lv1_vals), lv1_len, lv2,
				opcode, offs_buf, vals_buf);
		break;
	    case LGLSXP: case INTSXP:
		ans_len = sparse_Compare_ints_lv2(
				offs1, INTEGER(lv1_vals), lv1_len, lv2,
				opcode, offs_buf, vals_buf);
		break;
	    case REALSXP:
		ans_len = sparse_Compare_doubles_lv2(
				offs1, REAL(lv1_vals), lv1_len, lv2,
				opcode, offs_buf, vals_buf);
		break;
	    case CPLXSXP:
		ans_len = sparse_Compare_Rcomplexes_lv2(
				offs1, COMPLEX(lv1_vals), lv1_len, lv2,
				opcode, offs_buf, vals_buf);
		break;
	    default:
		error("SparseArray internal error in "
		      "_Compare_lv1_lv2():\n"
		      "    unsupported 'Rtype1': \"%s\"", type2char(Rtype1));
	}
	return _make_leaf_vector_from_bufs(LGLSXP,
					   offs_buf, vals_buf, ans_len);
}

