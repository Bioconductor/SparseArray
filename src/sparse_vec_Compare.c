/****************************************************************************
 *                  'Compare' operations on sparse vectors                  *
 ****************************************************************************/
#include "sparse_vec_Compare.h"

#include "sparse_vec.h"


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


/****************************************************************************
 * Core Compare_<Ltype>_<Rtype>() functions (10 in total)
 */

static inline int Compare_Rbyte_Rbyte(int opcode, Rbyte x, Rbyte y)
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

/* WARNING: Only valid to use on the int values of an integer vector (INTSXP).
   Should NOT be used on the int values of a logical vector (LGLSEXP), on
   which Compare_Rbyte_int() will generally produce wrong results!
   For example if Rbyte value 'x' is >= 2 and int value 'y' is 1 (TRUE),
   then == and <= will both return 0 when they are expected to return 1. */
static inline int Compare_Rbyte_int(int opcode, Rbyte x, int y)
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

static inline int Compare_Rbyte_double(int opcode, Rbyte x, double y)
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

static inline int Compare_Rbyte_Rcomplex(int opcode, Rbyte x, Rcomplex y)
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

static inline int Compare_int_int(int opcode, int x, int y)
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

static inline int Compare_int_double(int opcode, int x, double y)
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

static inline int Compare_int_Rcomplex(int opcode, int x, Rcomplex y)
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

static inline int Compare_double_double(int opcode, double x, double y)
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

static inline int Compare_double_Rcomplex(int opcode, double x, Rcomplex y)
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

static inline int Compare_Rcomplex_Rcomplex(int opcode, Rcomplex x, Rcomplex y)
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


/****************************************************************************
 * Two macros to generate the sv_Compare_*() functions (20 in total)
 */

/* Expands to def of sv_Compare_<Ltype>s_<Rtype>() function. */
#define FUNDEF_sv_Compare1(Ltype, Rtype)(int opcode,			\
		const struct sparse_vec *sv1,				\
		Rtype y,						\
		int *out_nzoffs, int *out_nzvals)			\
{									\
	const Ltype *nzvals1 = get_ ## Ltype ## _nzvals(sv1);		\
	int nzcount = 0;						\
	for (int k = 0; k < sv1->nzcount; k++) {			\
		Ltype x = nzvals1[k];					\
		int v = Compare_ ## Ltype ## _ ## Rtype			\
					(opcode, x, y);			\
		if (v != 0) {						\
			out_nzoffs[nzcount] = sv1->nzoffs[k];		\
			out_nzvals[nzcount] = v;			\
			nzcount++;					\
		}							\
	}								\
	return nzcount;							\
}

/* Expands to def of sv_Compare_<Ltype>s_<Rtype>s() function. */
#define FUNDEF_sv_Compare2(Ltype, Rtype)(int opcode,			\
		const struct sparse_vec *sv1,				\
		const struct sparse_vec *sv2,				\
		int *out_nzoffs, int *out_nzvals)			\
{									\
	int k1, k2, off;						\
	Ltype x;							\
	Rtype y;							\
									\
	int nzcount = 0;						\
	k1 = k2 = 0;							\
	while (next_nzvals_ ## Ltype ## _ ## Rtype			\
		(sv1, sv2, &k1, &k2, &off, &x, &y))			\
	{								\
		int v = Compare_ ## Ltype ## _ ## Rtype			\
					(opcode, x, y);			\
		if (v != 0) {						\
			out_nzoffs[nzcount] = off;			\
			out_nzvals[nzcount] = v;			\
			nzcount++;					\
		}							\
	}								\
	return nzcount;							\
}

static int sv_Compare_Rbytes_Rbyte
	FUNDEF_sv_Compare1(Rbyte, Rbyte)
static int sv_Compare_Rbytes_Rbytes
	FUNDEF_sv_Compare2(Rbyte, Rbyte)

static int sv_Compare_Rbytes_int
	FUNDEF_sv_Compare1(Rbyte, int)
static int sv_Compare_Rbytes_ints
	FUNDEF_sv_Compare2(Rbyte, int)

static int sv_Compare_Rbytes_double
	FUNDEF_sv_Compare1(Rbyte, double)
static int sv_Compare_Rbytes_doubles
	FUNDEF_sv_Compare2(Rbyte, double)

static int sv_Compare_Rbytes_Rcomplex
	FUNDEF_sv_Compare1(Rbyte, Rcomplex)
static int sv_Compare_Rbytes_Rcomplexes
	FUNDEF_sv_Compare2(Rbyte, Rcomplex)

static int sv_Compare_ints_int
	FUNDEF_sv_Compare1(int, int)
static int sv_Compare_ints_ints
	FUNDEF_sv_Compare2(int, int)

static int sv_Compare_ints_double
	FUNDEF_sv_Compare1(int, double)
static int sv_Compare_ints_doubles
	FUNDEF_sv_Compare2(int, double)

static int sv_Compare_ints_Rcomplex
	FUNDEF_sv_Compare1(int, Rcomplex)
static int sv_Compare_ints_Rcomplexes
	FUNDEF_sv_Compare2(int, Rcomplex)

static int sv_Compare_doubles_double
	FUNDEF_sv_Compare1(double, double)
static int sv_Compare_doubles_doubles
	FUNDEF_sv_Compare2(double, double)

static int sv_Compare_doubles_Rcomplex
	FUNDEF_sv_Compare1(double, Rcomplex)
static int sv_Compare_doubles_Rcomplexes
	FUNDEF_sv_Compare2(double, Rcomplex)

static int sv_Compare_Rcomplexes_Rcomplex
	FUNDEF_sv_Compare1(Rcomplex, Rcomplex)
static int sv_Compare_Rcomplexes_Rcomplexes
	FUNDEF_sv_Compare2(Rcomplex, Rcomplex)


/****************************************************************************
 * sv_Compare_Rbytes_scalar()
 * sv_Compare_ints_scalar()
 * sv_Compare_doubles_scalar()
 * sv_Compare_Rcomplexes_scalar()
 *
 * All these functions assume that 'scalar' is an atomic vector of length 1.
 * They do NOT check it!
 */

/* Type of 'scalar' is expected to be "raw" or bigger. */
static int sv_Compare_Rbytes_scalar(int opcode,
		const struct sparse_vec *sv1, SEXP scalar,
		int *out_nzoffs, int *out_nzvals)
{
	SEXPTYPE Rtype2 = TYPEOF(scalar);
	switch (Rtype2) {
	    case RAWSXP:
		return sv_Compare_Rbytes_Rbyte(opcode,
					sv1, RAW(scalar)[0],
					out_nzoffs, out_nzvals);
	    case LGLSXP: case INTSXP:
		return sv_Compare_Rbytes_int(opcode,
					sv1, INTEGER(scalar)[0],
					out_nzoffs, out_nzvals);
	    case REALSXP:
		return sv_Compare_Rbytes_double(opcode,
					sv1, REAL(scalar)[0],
					out_nzoffs, out_nzvals);
	    case CPLXSXP:
		return sv_Compare_Rbytes_Rcomplex(opcode,
					sv1, COMPLEX(scalar)[0],
					out_nzoffs, out_nzvals);
	}
	error("SparseArray internal error in "
	      "sv_Compare_Rbytes_scalar():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
	return 0;  /* will never reach this */
}

/* Type of 'scalar' is expected to be "integer" or bigger. */
static int sv_Compare_ints_scalar(int opcode,
		const struct sparse_vec *sv1, SEXP scalar,
		int *out_nzoffs, int *out_nzvals)
{
	SEXPTYPE Rtype2 = TYPEOF(scalar);
	switch (Rtype2) {
	    case LGLSXP: case INTSXP:
		return sv_Compare_ints_int(opcode,
					sv1, INTEGER(scalar)[0],
					out_nzoffs, out_nzvals);
	    case REALSXP:
		return sv_Compare_ints_double(opcode,
					sv1, REAL(scalar)[0],
					out_nzoffs, out_nzvals);
	    case CPLXSXP:
		return sv_Compare_ints_Rcomplex(opcode,
					sv1, COMPLEX(scalar)[0],
					out_nzoffs, out_nzvals);
	}
	error("SparseArray internal error in "
	      "sv_Compare_ints_scalar():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
	return 0;  /* will never reach this */
}

/* Type of 'scalar' is expected to be "double" or bigger. */
static int sv_Compare_doubles_scalar(int opcode,
		const struct sparse_vec *sv1, SEXP scalar,
		int *out_nzoffs, int *out_nzvals)
{
	SEXPTYPE Rtype2 = TYPEOF(scalar);
	switch (Rtype2) {
	    case REALSXP:
		return sv_Compare_doubles_double(opcode,
					sv1, REAL(scalar)[0],
					out_nzoffs, out_nzvals);
	    case CPLXSXP:
		return sv_Compare_doubles_Rcomplex(opcode,
					sv1, COMPLEX(scalar)[0],
					out_nzoffs, out_nzvals);
	}
	error("SparseArray internal error in "
	      "sv_Compare_doubles_scalar():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
	return 0;  /* will never reach this */
}

/* Type of 'scalar' is expected to be "complex". */
static int sv_Compare_Rcomplexes_scalar(int opcode,
		const struct sparse_vec *sv1, SEXP scalar,
		int *out_nzoffs, int *out_nzvals)
{
	SEXPTYPE Rtype2 = TYPEOF(scalar);
	switch (Rtype2) {
	    case CPLXSXP:
		return sv_Compare_Rcomplexes_Rcomplex(opcode,
					sv1, COMPLEX(scalar)[0],
					out_nzoffs, out_nzvals);
	}
	error("SparseArray internal error in "
	      "sv_Compare_Rcomplexes_scalar():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
	return 0;  /* will never reach this */
}


/****************************************************************************
 * sv_Compare_Rbytes_ANY()
 * sv_Compare_ints_ANY()
 * sv_Compare_doubles_ANY()
 * sv_Compare_Rcomplexes_ANY()
 */

static int sv_Compare_Rbytes_ANY(int opcode,
		const struct sparse_vec *sv1,
		const struct sparse_vec *sv2,
		int *out_nzoffs, int *out_nzvals)
{
	switch (sv2->Rtype) {
	    case RAWSXP:
		return sv_Compare_Rbytes_Rbytes(opcode,
				sv1, sv2, out_nzoffs, out_nzvals);
	    case LGLSXP: case INTSXP:
		return sv_Compare_Rbytes_ints(opcode,
				sv1, sv2, out_nzoffs, out_nzvals);
	    case REALSXP:
		return sv_Compare_Rbytes_doubles(opcode,
				sv1, sv2, out_nzoffs, out_nzvals);
	    case CPLXSXP:
		return sv_Compare_Rbytes_Rcomplexes(opcode,
				sv1, sv2, out_nzoffs, out_nzvals);
	}
	error("SparseArray internal error in "
	      "sv_Compare_Rbytes_ANY():\n"
	      "    unsupported 'sv2->Rtype': \"%s\"", type2char(sv2->Rtype));
	return 0;  /* will never reach this */
}

static int sv_Compare_ints_ANY(int opcode,
		const struct sparse_vec *sv1,
		const struct sparse_vec *sv2,
		int *out_nzoffs, int *out_nzvals)
{
	switch (sv2->Rtype) {
	    case RAWSXP:
		return sv_Compare_Rbytes_ints(flip_opcode(opcode),
				sv2, sv1, out_nzoffs, out_nzvals);
	    case LGLSXP: case INTSXP:
		return sv_Compare_ints_ints(opcode,
				sv1, sv2, out_nzoffs, out_nzvals);
	    case REALSXP:
		return sv_Compare_ints_doubles(opcode,
				sv1, sv2, out_nzoffs, out_nzvals);
	    case CPLXSXP:
		return sv_Compare_ints_Rcomplexes(opcode,
				sv1, sv2, out_nzoffs, out_nzvals);
	}
	error("SparseArray internal error in "
	      "sv_Compare_ints_ANY():\n"
	      "    unsupported 'sv2->Rtype': \"%s\"", type2char(sv2->Rtype));
	return 0;  /* will never reach this */
}

static int sv_Compare_doubles_ANY(int opcode,
		const struct sparse_vec *sv1,
		const struct sparse_vec *sv2,
		int *out_nzoffs, int *out_nzvals)
{
	switch (sv2->Rtype) {
	    case RAWSXP:
		return sv_Compare_Rbytes_doubles(flip_opcode(opcode),
				sv2, sv1, out_nzoffs, out_nzvals);
	    case LGLSXP: case INTSXP:
		return sv_Compare_ints_doubles(flip_opcode(opcode),
				sv2, sv1, out_nzoffs, out_nzvals);
	    case REALSXP:
		return sv_Compare_doubles_doubles(opcode,
				sv1, sv2, out_nzoffs, out_nzvals);
	    case CPLXSXP:
		return sv_Compare_doubles_Rcomplexes(opcode,
				sv1, sv2, out_nzoffs, out_nzvals);
	}
	error("SparseArray internal error in "
	      "sv_Compare_doubles_ANY():\n"
	      "    unsupported 'sv2->Rtype': \"%s\"", type2char(sv2->Rtype));
	return 0;  /* will never reach this */
}

static int sv_Compare_Rcomplexes_ANY(int opcode,
		const struct sparse_vec *sv1,
		const struct sparse_vec *sv2,
		int *out_nzoffs, int *out_nzvals)
{
	switch (sv2->Rtype) {
	    case RAWSXP:
		return sv_Compare_Rbytes_Rcomplexes(flip_opcode(opcode),
				sv2, sv1, out_nzoffs, out_nzvals);
	    case LGLSXP: case INTSXP:
		return sv_Compare_ints_Rcomplexes(flip_opcode(opcode),
				sv2, sv1, out_nzoffs, out_nzvals);
	    case REALSXP:
		return sv_Compare_doubles_Rcomplexes(flip_opcode(opcode),
				sv2, sv1, out_nzoffs, out_nzvals);
	    case CPLXSXP:
		return sv_Compare_Rcomplexes_Rcomplexes(opcode,
				sv1, sv2, out_nzoffs, out_nzvals);
	}
	error("SparseArray internal error in "
	      "sv_Compare_Rcomplexes_ANY():\n"
	      "    unsupported 'sv2->Rtype': \"%s\"", type2char(sv2->Rtype));
	return 0;  /* will never reach this */
}


/****************************************************************************
 * _sparse_vec_Compare_sv1_zero()
 * _sparse_vec_Compare_sv1_scalar()
 * _sparse_vec_Compare_sv1_sv2()
 */

int _sparse_vec_Compare_sv1_zero(int opcode,
		const struct sparse_vec *sv1,
		int *out_nzoffs, int *out_nzvals)
{
	switch (sv1->Rtype) {
	    case RAWSXP:
		return sv_Compare_Rbytes_Rbyte(opcode,
				sv1, Rbyte0, out_nzoffs, out_nzvals);
	    case LGLSXP: case INTSXP:
		return sv_Compare_ints_int(opcode,
				sv1, int0, out_nzoffs, out_nzvals);
	    case REALSXP:
		return sv_Compare_doubles_double(opcode,
				sv1, double0, out_nzoffs, out_nzvals);
	    case CPLXSXP:
		return sv_Compare_Rcomplexes_Rcomplex(opcode,
				sv1, Rcomplex0, out_nzoffs, out_nzvals);
	}
	error("SparseArray internal error in "
	      "_sparse_vec_Compare_sv1_zero():\n"
	      "    unsupported 'sv1->Rtype': \"%s\"", type2char(sv1->Rtype));
	return 0;  /* will never reach this */
}

/* 'scalar' is assumed to be an atomic vector of length 1. This is NOT checked!
   Also its type is expected to be the same size as sv1's type or bigger. */
int _sparse_vec_Compare_sv1_scalar(int opcode,
		const struct sparse_vec *sv1,
		SEXP scalar,
		int *out_nzoffs, int *out_nzvals)
{
	switch (sv1->Rtype) {
	    case RAWSXP:
		return sv_Compare_Rbytes_scalar(opcode,
				sv1, scalar, out_nzoffs, out_nzvals);
	    case LGLSXP: case INTSXP:
		return sv_Compare_ints_scalar(opcode,
				sv1, scalar, out_nzoffs, out_nzvals);
	    case REALSXP:
		return sv_Compare_doubles_scalar(opcode,
				sv1, scalar, out_nzoffs, out_nzvals);
	    case CPLXSXP:
		return sv_Compare_Rcomplexes_scalar(opcode,
				sv1, scalar, out_nzoffs, out_nzvals);
	}
	error("SparseArray internal error in "
	      "_sparse_vec_Compare_sv1_scalar():\n"
	      "    unsupported 'sv1->Rtype': \"%s\"", type2char(sv1->Rtype));
	return 0;  /* will never reach this */
}

int _sparse_vec_Compare_sv1_sv2(int opcode,
		const struct sparse_vec *sv1,
		const struct sparse_vec *sv2,
		int *out_nzoffs, int *out_nzvals)
{
	switch (sv1->Rtype) {
	    case RAWSXP:
		return sv_Compare_Rbytes_ANY(opcode,
				sv1, sv2, out_nzoffs, out_nzvals);
	    case LGLSXP: case INTSXP:
		return sv_Compare_ints_ANY(opcode,
				sv1, sv2, out_nzoffs, out_nzvals);
	    case REALSXP:
		return sv_Compare_doubles_ANY(opcode,
				sv1, sv2, out_nzoffs, out_nzvals);
	    case CPLXSXP:
		return sv_Compare_Rcomplexes_ANY(opcode,
				sv1, sv2, out_nzoffs, out_nzvals);
	}
	error("SparseArray internal error in "
	      "_sparse_vec_Compare_sv1_sv2():\n"
	      "    unsupported 'sv1->Rtype': \"%s\"", type2char(sv1->Rtype));
	return 0;  /* will never reach this */
}

