/****************************************************************************
 *                  'Compare' operations on sparse vectors                  *
 ****************************************************************************/
#include "SparseVec_Compare.h"

#include "SparseVec.h"


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
 *
 * All these functions must return 'int0' (FALSE), 'int1' (TRUE),
 * or 'intNA' (NA).
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
   Should NOT be used on the int values of a logical vector (LGLSXP), on
   which Compare_Rbyte_int() will generally produce wrong results!
   For example if Rbyte value 'x' is >= 2 and int value 'y' is 1 (TRUE),
   then == and <= will both return 0 when they are expected to return 1. */
static inline int Compare_Rbyte_int(int opcode, Rbyte x, int y)
{
	int x1;

	if (y == NA_INTEGER)
		return intNA;
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
		return intNA;
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
		return intNA;
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
		return intNA;
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
		return intNA;
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
		return intNA;
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
		return intNA;
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
		return intNA;
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
		return intNA;
	switch (opcode) {
	    case EQ_OPCODE: return x.r == y.r && x.i == y.i;
	    case NE_OPCODE: return x.r != y.r || x.i != y.i;
	}
	error("SparseArray internal error in Compare_Rcomplex_Rcomplex():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}


/****************************************************************************
 * Two macros to generate the code (arg list + body) of the following
 * functions:
 *   - Compare_<Ltype>SV_<Rtype>() (10 functions)
 *   - Compare_<Ltype>SV_<Rtype>SV() (10 functions)
 */

/* Generate code of Compare_<Ltype>SV_<Rtype>() functions.
   These functions should never be called on 'sv1' and 'y' when 'sv1'
   has a background set to zero and 'y' is NA or NaN. */
#define FUNDEF_Compare_LtypeSV_Rtype(Ltype, Rtype)(int opcode,		\
		const SparseVec *sv1, Rtype y,				\
		SparseVec *out_sv)					\
{									\
        if (out_sv->len != sv1->len)					\
                error("SparseArray internal error in "			\
                      "Compare_<Ltype>SV_<Rtype>():\n"			\
                      "    'sv1' and 'out_sv' are incompatible");	\
	int *out_nzvals = (int *) out_sv->nzvals;			\
	out_sv->nzcount = 0;						\
	int out_background = out_sv->na_background ? intNA : int0;	\
	const Ltype *nzvals1_p = get_ ## Ltype ## SV_nzvals_p(sv1);	\
	if (nzvals1_p == NULL) {  /* lacunar SparseVec */		\
		int out_val = Compare_ ## Ltype ## _ ## Rtype		\
					(opcode, Ltype ## 1, y);	\
		if (out_val == out_background)				\
			return;						\
		/* What 'out_val' is expected to be at this point    */	\
		/* depends on 'out_background':                      */	\
		/* - If 'out_background' is 'int0' then 'sv1' also   */	\
		/*   has a background set to zero so 'y' cannot      */	\
		/*   be NA or NaN (i.e. is.na(y) must be FALSE).     */	\
		/*   This means that 'out_val' can only be TRUE      */	\
		/*   (i.e. 'int1'). In particular 'out_val' cannot   */	\
		/*   be NA (i.e. 'intNA') or FALSE (i.e. 'int0').    */	\
		/* - If background is NA then 'out_val' can be TRUE  */	\
		/*   or FALSE. It cannot be NA.                      */	\
		out_nzvals[0] = out_val;				\
		out_sv->nzcount = PROPAGATE_NZOFFS;			\
		return;							\
	}								\
	/* regular SparseVec */						\
	int nzcount1 = get_SV_nzcount(sv1);				\
	for (int k = 0; k < nzcount1; k++) {				\
		int out_val = Compare_ ## Ltype ## _ ## Rtype		\
					(opcode, nzvals1_p[k], y);	\
		if (out_val == out_background)				\
			continue;					\
		APPEND_TO_NZVALS_NZOFFS(out_val, sv1->nzoffs[k],	\
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);	\
	}								\
	return;								\
}

/* Generate code of Compare_<Ltype>SV_<Rtype>SV() functions. */
#define FUNDEF_Compare_LtypeSV_RtypeSV(Ltype, Rtype)(int opcode,	\
		const SparseVec *sv1, const SparseVec *sv2,		\
		SparseVec *out_sv)					\
{									\
        if (out_sv->len != sv1->len || out_sv->len != sv2->len)		\
                error("SparseArray internal error in "			\
                      "Compare_<Ltype>SV_<Rtype>SV()():\n"		\
                      "    'sv1', 'sv2', and 'out_sv' are incompatible"); \
	int *out_nzvals = (int *) out_sv->nzvals;			\
	out_sv->nzcount = 0;						\
	int out_background = out_sv->na_background ? intNA : int0;	\
	int k1 = 0, k2 = 0;						\
	int off;							\
	Ltype x;							\
	Rtype y;							\
	while (next_ ## Ltype ## _ ## Rtype ## _vals			\
		(sv1, sv2, &k1, &k2, &off, &x, &y))			\
	{								\
		int out_val = Compare_ ## Ltype ## _ ## Rtype		\
					(opcode, x, y);			\
		if (out_val == out_background)				\
			continue;					\
		APPEND_TO_NZVALS_NZOFFS(out_val, off,			\
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);   \
	}								\
	return;								\
}

static void Compare_RbyteSV_Rbyte
	FUNDEF_Compare_LtypeSV_Rtype(Rbyte, Rbyte)
static void Compare_RbyteSV_RbyteSV
	FUNDEF_Compare_LtypeSV_RtypeSV(Rbyte, Rbyte)

static void Compare_RbyteSV_int
	FUNDEF_Compare_LtypeSV_Rtype(Rbyte, int)
static void Compare_RbyteSV_intSV
	FUNDEF_Compare_LtypeSV_RtypeSV(Rbyte, int)

static void Compare_RbytesSV_double
	FUNDEF_Compare_LtypeSV_Rtype(Rbyte, double)
static void Compare_RbytesSV_doubleSV
	FUNDEF_Compare_LtypeSV_RtypeSV(Rbyte, double)

static void Compare_RbyteSV_Rcomplex
	FUNDEF_Compare_LtypeSV_Rtype(Rbyte, Rcomplex)
static void Compare_RbyteSV_RcomplexSV
	FUNDEF_Compare_LtypeSV_RtypeSV(Rbyte, Rcomplex)

static void Compare_intSV_int
	FUNDEF_Compare_LtypeSV_Rtype(int, int)
static void Compare_intSV_intSV
	FUNDEF_Compare_LtypeSV_RtypeSV(int, int)

static void Compare_intSV_double
	FUNDEF_Compare_LtypeSV_Rtype(int, double)
static void Compare_intSV_doubleSV
	FUNDEF_Compare_LtypeSV_RtypeSV(int, double)

static void Compare_intSV_Rcomplex
	FUNDEF_Compare_LtypeSV_Rtype(int, Rcomplex)
static void Compare_intSV_RcomplexSV
	FUNDEF_Compare_LtypeSV_RtypeSV(int, Rcomplex)

static void Compare_doubleSV_double
	FUNDEF_Compare_LtypeSV_Rtype(double, double)
static void Compare_doubleSV_doubleSV
	FUNDEF_Compare_LtypeSV_RtypeSV(double, double)

static void Compare_doubleSV_Rcomplex
	FUNDEF_Compare_LtypeSV_Rtype(double, Rcomplex)
static void Compare_doubleSV_RcomplexSV
	FUNDEF_Compare_LtypeSV_RtypeSV(double, Rcomplex)

static void Compare_RcomplexSV_Rcomplex
	FUNDEF_Compare_LtypeSV_Rtype(Rcomplex, Rcomplex)
static void Compare_RcomplexSV_RcomplexSV
	FUNDEF_Compare_LtypeSV_RtypeSV(Rcomplex, Rcomplex)


/****************************************************************************
 * Compare_RbyteSV_scalar()
 * Compare_intSV_scalar()
 * Compare_doubleSV_scalar()
 * Compare_RcomplexSV_scalar()
 *
 * All these functions assume that 'scalar' is an atomic SEXP of length 1.
 * They do NOT check it!
 */

/* Type of 'scalar' is expected to be "raw" or bigger. */
static void Compare_RbyteSV_scalar(int opcode,
		const SparseVec *sv1, SEXP scalar, SparseVec *out_sv)
{
	SEXPTYPE Rtype2 = TYPEOF(scalar);
	switch (Rtype2) {
	    case RAWSXP:
		Compare_RbyteSV_Rbyte(opcode,
				sv1, RAW(scalar)[0], out_sv);
		return;
	    case INTSXP: case LGLSXP:
		Compare_RbyteSV_int(opcode,
				sv1, INTEGER(scalar)[0], out_sv);
		return;
	    case REALSXP:
		Compare_RbytesSV_double(opcode,
				sv1, REAL(scalar)[0], out_sv);
		return;
	    case CPLXSXP:
		Compare_RbyteSV_Rcomplex(opcode,
				sv1, COMPLEX(scalar)[0], out_sv);
		return;
	}
	error("SparseArray internal error in "
	      "Compare_RbyteSV_scalar():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
}

/* Type of 'scalar' is expected to be "integer" or bigger. */
static void Compare_intSV_scalar(int opcode,
		const SparseVec *sv1, SEXP scalar, SparseVec *out_sv)
{
	SEXPTYPE Rtype2 = TYPEOF(scalar);
	switch (Rtype2) {
	    case INTSXP: case LGLSXP:
		Compare_intSV_int(opcode,
				sv1, INTEGER(scalar)[0], out_sv);
		return;
	    case REALSXP:
		Compare_intSV_double(opcode,
				sv1, REAL(scalar)[0], out_sv);
		return;
	    case CPLXSXP:
		Compare_intSV_Rcomplex(opcode,
				sv1, COMPLEX(scalar)[0], out_sv);
		return;
	}
	error("SparseArray internal error in "
	      "Compare_intSV_scalar():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
}

/* Type of 'scalar' is expected to be "double" or bigger. */
static void Compare_doubleSV_scalar(int opcode,
		const SparseVec *sv1, SEXP scalar, SparseVec *out_sv)
{
	SEXPTYPE Rtype2 = TYPEOF(scalar);
	switch (Rtype2) {
	    case REALSXP:
		Compare_doubleSV_double(opcode,
				sv1, REAL(scalar)[0], out_sv);
		return;
	    case CPLXSXP:
		Compare_doubleSV_Rcomplex(opcode,
				sv1, COMPLEX(scalar)[0], out_sv);
		return;
	}
	error("SparseArray internal error in "
	      "Compare_doubleSV_scalar():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
}

/* Type of 'scalar' is expected to be "complex". */
static void Compare_RcomplexSV_scalar(int opcode,
		const SparseVec *sv1, SEXP scalar, SparseVec *out_sv)
{
	SEXPTYPE Rtype2 = TYPEOF(scalar);
	switch (Rtype2) {
	    case CPLXSXP:
		Compare_RcomplexSV_Rcomplex(opcode, sv1, COMPLEX(scalar)[0],
					out_sv);
		return;
	}
	error("SparseArray internal error in "
	      "Compare_RcomplexSV_scalar():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
}


/****************************************************************************
 * Compare_RbyteSV_SV()
 * Compare_intSV_SV()
 * Compare_doubleSV_SV()
 * Compare_RcomplexSV_SV()
 */

static void Compare_RbyteSV_SV(int opcode,
		const SparseVec *sv1, const SparseVec *sv2,
		SparseVec *out_sv)
{
	SEXPTYPE Rtype2 = get_SV_Rtype(sv2);
	switch (Rtype2) {
	    case RAWSXP:
		Compare_RbyteSV_RbyteSV(opcode, sv1, sv2, out_sv);
		return;
	    case INTSXP: case LGLSXP:
		Compare_RbyteSV_intSV(opcode, sv1, sv2, out_sv);
		return;
	    case REALSXP:
		Compare_RbytesSV_doubleSV(opcode, sv1, sv2, out_sv);
		return;
	    case CPLXSXP:
		Compare_RbyteSV_RcomplexSV(opcode, sv1, sv2, out_sv);
		return;
	}
	error("SparseArray internal error in "
	      "Compare_RbyteSV_SV():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
}

static void Compare_intSV_SV(int opcode,
		const SparseVec *sv1, const SparseVec *sv2,
		SparseVec *out_sv)
{
	SEXPTYPE Rtype2 = get_SV_Rtype(sv2);
	switch (Rtype2) {
	    case RAWSXP:
		Compare_RbyteSV_intSV(flip_Compare_opcode(opcode),
				sv2, sv1, out_sv);
		return;
	    case INTSXP: case LGLSXP:
		Compare_intSV_intSV(opcode, sv1, sv2, out_sv);
		return;
	    case REALSXP:
		Compare_intSV_doubleSV(opcode, sv1, sv2, out_sv);
		return;
	    case CPLXSXP:
		Compare_intSV_RcomplexSV(opcode, sv1, sv2, out_sv);
		return;
	}
	error("SparseArray internal error in "
	      "Compare_intSV_SV():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
}

static void Compare_doubleSV_SV(int opcode,
		const SparseVec *sv1, const SparseVec *sv2,
		SparseVec *out_sv)
{
	SEXPTYPE Rtype2 = get_SV_Rtype(sv2);
	switch (Rtype2) {
	    case RAWSXP:
		Compare_RbytesSV_doubleSV(flip_Compare_opcode(opcode),
				sv2, sv1, out_sv);
		return;
	    case INTSXP: case LGLSXP:
		Compare_intSV_doubleSV(flip_Compare_opcode(opcode),
				sv2, sv1, out_sv);
		return;
	    case REALSXP:
		Compare_doubleSV_doubleSV(opcode, sv1, sv2, out_sv);
		return;
	    case CPLXSXP:
		Compare_doubleSV_RcomplexSV(opcode, sv1, sv2, out_sv);
		return;
	}
	error("SparseArray internal error in "
	      "Compare_doubleSV_SV():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
}

static void Compare_RcomplexSV_SV(int opcode,
		const SparseVec *sv1, const SparseVec *sv2,
		SparseVec *out_sv)
{
	SEXPTYPE Rtype2 = get_SV_Rtype(sv2);
	switch (Rtype2) {
	    case RAWSXP:
		Compare_RbyteSV_RcomplexSV(flip_Compare_opcode(opcode),
				sv2, sv1, out_sv);
		return;
	    case INTSXP: case LGLSXP:
		Compare_intSV_RcomplexSV(flip_Compare_opcode(opcode),
				sv2, sv1, out_sv);
		return;
	    case REALSXP:
		Compare_doubleSV_RcomplexSV(flip_Compare_opcode(opcode),
				sv2, sv1, out_sv);
		return;
	    case CPLXSXP:
		Compare_RcomplexSV_RcomplexSV(opcode, sv1, sv2, out_sv);
		return;
	}
	error("SparseArray internal error in "
	      "Compare_RcomplexSV_SV():\n"
	      "    unsupported 'Rtype2': \"%s\"", type2char(Rtype2));
}


/****************************************************************************
 * _Compare_sv1_zero()
 * _Compare_sv1_scalar()
 * _Compare_sv1_sv2()
 */

void _Compare_sv1_zero(int opcode, const SparseVec *sv1,
		SparseVec *out_sv)
{
	SEXPTYPE Rtype1 = get_SV_Rtype(sv1);
	switch (Rtype1) {
	    case RAWSXP:
		Compare_RbyteSV_Rbyte(opcode, sv1, Rbyte0, out_sv);
		return;
	    case INTSXP: case LGLSXP:
		Compare_intSV_int(opcode, sv1, int0, out_sv);
		return;
	    case REALSXP:
		Compare_doubleSV_double(opcode, sv1, double0, out_sv);
		return;
	    case CPLXSXP:
		Compare_RcomplexSV_Rcomplex(opcode, sv1, Rcomplex0, out_sv);
		return;
	}
	error("SparseArray internal error in "
	      "_Compare_sv1_zero():\n"
	      "    unsupported 'Rtype1': \"%s\"", type2char(Rtype1));
}

/* 'scalar' is assumed to be an atomic vector of length 1. This is NOT checked!
   Also its type is expected to be the same as sv1's type or bigger. */
void _Compare_sv1_scalar(int opcode, const SparseVec *sv1, SEXP scalar,
		SparseVec *out_sv)
{
	SEXPTYPE Rtype1 = get_SV_Rtype(sv1);
	switch (Rtype1) {
	    case RAWSXP:
		Compare_RbyteSV_scalar(opcode, sv1, scalar, out_sv);
		return;
	    case INTSXP: case LGLSXP:
		Compare_intSV_scalar(opcode, sv1, scalar, out_sv);
		return;
	    case REALSXP:
		Compare_doubleSV_scalar(opcode, sv1, scalar, out_sv);
		return;
	    case CPLXSXP:
		Compare_RcomplexSV_scalar(opcode, sv1, scalar, out_sv);
		return;
	}
	error("SparseArray internal error in "
	      "_Compare_sv1_scalar():\n"
	      "    unsupported 'Rtype1': \"%s\"", type2char(Rtype1));
}

void _Compare_sv1_sv2(int opcode, const SparseVec *sv1, const SparseVec *sv2,
		SparseVec *out_sv)
{
	SEXPTYPE Rtype1 = get_SV_Rtype(sv1);
	switch (Rtype1) {
	    case RAWSXP:
		Compare_RbyteSV_SV(opcode, sv1, sv2, out_sv);
		return;
	    case INTSXP: case LGLSXP:
		Compare_intSV_SV(opcode, sv1, sv2, out_sv);
		return;
	    case REALSXP:
		Compare_doubleSV_SV(opcode, sv1, sv2, out_sv);
		return;
	    case CPLXSXP:
		Compare_RcomplexSV_SV(opcode, sv1, sv2, out_sv);
		return;
	}
	error("SparseArray internal error in "
	      "_Compare_sv1_sv2():\n"
	      "    unsupported 'Rtype1': \"%s\"", type2char(Rtype1));
}

