/****************************************************************************
 *                   'Arith' operations on sparse vectors                   *
 ****************************************************************************/
#include "SparseVec_Arith.h"

#include "SparseVec.h"

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

/* Does not support 'opcode' values DIV_OPCODE ("/") or POW_OPCODE ("^"). */
static inline int iarith(int opcode, int x, int y, int *ovflow)
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
		   deviates from C modulo operator %), so we adjust to
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
		error("SparseArray internal error in iarith():\n"
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
   a. In R, x ^ 0 and 1 ^ y are always 1, even when x or y is NA.
   b. If x < 0 and y is infinite (i.e. R_PosInf or R_NegInf), then x ^ y
      always returns NaN (even when -1 < x < 0), whereas pow(x, y) returns
      the following:
      - pow(x, R_PosInf) returns R_PosInf if x < -1, or 1 if x == -1, or 0
        if -1 < x < 1;
      - pow(x, R_NegInf) returns 0 if x < -1, or 1 if x == -1, or R_PosInf
        if -1 < x < 1.
   c. x ^ y is expected to return NaN for any noninteger y. However, in
      the specific case where x is R_NegInf, pow(x, y) returns R_PosInf
      for any noninteger y on my Intel Ubuntu 22.04 laptop. */
static inline double Rpow_double(double x, double y)
{
	if (x == 1.0 || y == 0.0)
		return 1.0;
	if (R_IsNaN(y) ||
	    (x < 0.0 && (y == R_PosInf || y == R_NegInf)) ||
	    (x == R_NegInf && y != trunc(y)))
		return R_NaN;
	return pow(x, y);
}

static inline double Rmod_double(double x, double y)
{
	if (y == 0.0 || x == R_PosInf || x == R_NegInf)
		return R_NaN;
	if (x == 0.0)
		return ISNAN(y) ? y : 0.0;
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

static inline double darith_double_double(int opcode, double x, double y)
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
	error("SparseArray internal error in darith_double_double():\n"
	      "    unsupported 'opcode'");
	return 0.0;  /* will never reach this */
}

static inline double darith_double_int(int opcode, double x, int y)
{
	double yy = y == intNA ? doubleNA : (double) y;
	return darith_double_double(opcode, x, yy);
}

static inline double darith_int_double(int opcode, int x, double y)
{
	double xx = x == intNA ? doubleNA : (double) x;
	return darith_double_double(opcode, xx, y);
}

static inline double darith_int_int(int opcode, int x, int y)
{
	double xx = x == intNA ? doubleNA : (double) x;
	double yy = y == intNA ? doubleNA : (double) y;
	return darith_double_double(opcode, xx, yy);
}

static inline void check_outRtype(SEXPTYPE expected_outRtype,
				  SEXPTYPE effective_outRtype,
				  const char *fun)
{
	if (expected_outRtype == effective_outRtype)
		return;
	error("SparseArray internal error in %s():\n"
	      "    expected_outRtype (\"%s\") != effective_outRtype (\"%s\")",
	      fun, type2char(expected_outRtype), type2char(effective_outRtype));
}


/****************************************************************************
 * dArith_intSV_ints()
 * dArith_intSV_doubles()
 * dArith_doubleSV_ints()
 * dArith_doubleSV_doubles()
 * iArith_intSV_ints()
 *
 * dArith_ints_intSV()
 * dArith_doubles_intSV()
 * dArith_ints_doubleSV()
 * dArith_doubles_doubleSV()
 * iArith_ints_intSV()
 */

#define DEFINE_dArith_SV_y_FUN(funname, Ltype, Rtype)			    \
static void dArith_ ## Ltype ## SV_ ## Rtype ## s(int opcode,		    \
	const SparseVec *sv1, const Rtype *y, int y_len,		    \
	SparseVec *out_sv)						    \
{									    \
	if (sv1->len != out_sv->len)					    \
		error("SparseArray internal error in %s():\n"		    \
		      "    'sv1' and 'out_sv' are incompatible",	    \
		      funname);						    \
	if (sv1->len != 0 && y_len == 0)				    \
		error("SparseArray internal error in %s():\n"		    \
		      "    'y_len' cannot be 0 unless 'sv1->len' is 0",	    \
		      funname);						    \
	check_outRtype(get_SV_Rtype(out_sv), REALSXP, funname);		    \
	double (*darith_fun)(int, Ltype, Rtype);			    \
	darith_fun = &darith_ ## Ltype ## _ ## Rtype;			    \
	double *out_nzvals = (double *) out_sv->nzvals;			    \
	out_sv->nzcount = 0;						    \
	const Ltype *nzvals1_p = get_ ## Ltype ## SV_nzvals_p(sv1);	    \
	if (nzvals1_p == NULL && y_len == 1) {				    \
		/* shortcut for "lacunar SparseVec <op> scalar" case */	    \
		double out_val = darith_fun(opcode, Ltype ## 1, y[0]);	    \
		if (IS_BG_DOUBLE(out_val, out_sv->na_background))	    \
			return;						    \
		out_nzvals[0] = out_val;				    \
		out_sv->nzcount = PROPAGATE_NZOFFS;			    \
		return;							    \
	}								    \
	int nzcount1 = get_SV_nzcount(sv1);				    \
	for (int k = 0; k < nzcount1; k++) {				    \
		Ltype x = nzvals1_p == NULL ? Ltype ## 1 : nzvals1_p[k];    \
		int nzoff1 = sv1->nzoffs[k];				    \
		Rtype yy = y[nzoff1 % y_len];				    \
		double out_val = darith_fun(opcode, x, yy);		    \
		if (IS_BG_DOUBLE(out_val, out_sv->na_background))	    \
			continue;					    \
		APPEND_TO_NZVALS_NZOFFS(out_val, nzoff1,		    \
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);	    \
	}								    \
}

DEFINE_dArith_SV_y_FUN("dArith_intSV_ints", int, int)
DEFINE_dArith_SV_y_FUN("dArith_intSV_doubles", int, double)
DEFINE_dArith_SV_y_FUN("dArith_doubleSV_ints", double, int)
DEFINE_dArith_SV_y_FUN("dArith_doubleSV_doubles", double, double)

static void iArith_intSV_ints(int opcode,
		const SparseVec *sv1, const int *y, int y_len,
		SparseVec *out_sv, int *ovflow)
{
	if (sv1->len != out_sv->len)
		error("SparseArray internal error in "
		      "iArith_intSV_ints():\n"
		      "    'sv1' and 'out_sv' are incompatible");
	if (sv1->len != 0 && y_len == 0)
		error("SparseArray internal error in "
		      "iArith_intSV_ints():\n"
		      "    'y_len' cannot be 0 unless 'sv1->len' is 0");
	check_outRtype(out_sv->Rtype, INTSXP, "iArith_intSV_ints");
	int *out_nzvals = (int *) out_sv->nzvals;
	out_sv->nzcount = 0;
	int out_bg_val = out_sv->na_background ? intNA : int0;
	const int *nzvals1_p = get_intSV_nzvals_p(sv1);
	if (nzvals1_p == NULL && y_len == 1) {
		/* shortcut for "lacunar SparseVec <op> scalar" case */
		int out_val = iarith(opcode, int1, y[0], ovflow);
		if (out_val == out_bg_val)
			return;
		out_nzvals[0] = out_val;
		out_sv->nzcount = PROPAGATE_NZOFFS;
		return;
	}
	/* regular SparseVec */
	int nzcount1 = get_SV_nzcount(sv1);
	for (int k = 0; k < nzcount1; k++) {
		int x = nzvals1_p == NULL ? int1 : nzvals1_p[k];
		int nzoff1 = sv1->nzoffs[k];
		int yy = y[nzoff1 % y_len];
		int out_val = iarith(opcode, x, yy, ovflow);
		if (out_val == out_bg_val)
			continue;
		APPEND_TO_NZVALS_NZOFFS(out_val, nzoff1,
				out_nzvals, out_sv->nzoffs, out_sv->nzcount);
	}
	return;
}

#define DEFINE_dArith_x_SV_FUN(funname, Ltype, Rtype)			    \
static void dArith_ ## Ltype ## s_ ## Rtype ## SV(int opcode,		    \
	const Ltype *x, int x_len, const SparseVec *sv2,		    \
	SparseVec *out_sv)						    \
{									    \
	if (sv2->len != out_sv->len)					    \
		error("SparseArray internal error in %s():\n"		    \
		      "    'sv2' and 'out_sv' are incompatible",	    \
		      funname);						    \
	if (sv2->len != 0 && x_len == 0)				    \
		error("SparseArray internal error in %s():\n"		    \
		      "    'x_len' cannot be 0 unless 'sv2->len' is 0",	    \
		      funname);						    \
	check_outRtype(get_SV_Rtype(out_sv), REALSXP, funname);		    \
	double (*darith_fun)(int, Ltype, Rtype);			    \
	darith_fun = &darith_ ## Ltype ## _ ## Rtype;			    \
	double *out_nzvals = (double *) out_sv->nzvals;			    \
	out_sv->nzcount = 0;						    \
	const Rtype *nzvals2_p = get_ ## Rtype ## SV_nzvals_p(sv2);	    \
	if (nzvals2_p == NULL && x_len == 1) {				    \
		/* shortcut for "scalar <op> lacunar SparseVec" case */	    \
		double out_val = darith_fun(opcode, x[0], Rtype ## 1);	    \
		if (IS_BG_DOUBLE(out_val, out_sv->na_background))	    \
			return;						    \
		out_nzvals[0] = out_val;				    \
		out_sv->nzcount = PROPAGATE_NZOFFS;			    \
		return;							    \
	}								    \
	int nzcount2 = get_SV_nzcount(sv2);				    \
	for (int k = 0; k < nzcount2; k++) {				    \
		int nzoff2 = sv2->nzoffs[k];				    \
		Ltype xx = x[nzoff2 % x_len];				    \
		Rtype y = nzvals2_p == NULL ? Rtype ## 1 : nzvals2_p[k];    \
		double out_val = darith_fun(opcode, xx, y);		    \
		if (IS_BG_DOUBLE(out_val, out_sv->na_background))	    \
			continue;					    \
		APPEND_TO_NZVALS_NZOFFS(out_val, nzoff2,		    \
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);	    \
	}								    \
}

DEFINE_dArith_x_SV_FUN("dArith_ints_intSV", int, int)
DEFINE_dArith_x_SV_FUN("dArith_doubles_intSV", double, int)
DEFINE_dArith_x_SV_FUN("dArith_ints_doubleSV", int, double)
DEFINE_dArith_x_SV_FUN("dArith_doubles_doubleSV", double, double)

static void iArith_ints_intSV(int opcode,
		const int *x, int x_len, const SparseVec *sv2,
		SparseVec *out_sv, int *ovflow)
{
	if (sv2->len != out_sv->len)
		error("SparseArray internal error in "
		      "iArith_ints_intSV():\n"
		      "    'sv2' and 'out_sv' are incompatible");
	if (sv2->len != 0 && x_len == 0)
		error("SparseArray internal error in "
		      "iArith_ints_intSV():\n"
		      "    'x_len' cannot be 0 unless 'sv2->len' is 0");
	check_outRtype(out_sv->Rtype, INTSXP, "iArith_ints_intSV");
	int *out_nzvals = (int *) out_sv->nzvals;
	out_sv->nzcount = 0;
	int out_bg_val = out_sv->na_background ? intNA : int0;
	const int *nzvals2_p = get_intSV_nzvals_p(sv2);
	if (nzvals2_p == NULL && x_len == 1) {
		/* shortcut for "scalar <op> lacunar SparseVec" case */
		int out_val = iarith(opcode, x[0], int1, ovflow);
		if (out_val == out_bg_val)
			return;
		out_nzvals[0] = out_val;
		out_sv->nzcount = PROPAGATE_NZOFFS;
		return;
	}
	/* regular SparseVec */
	int nzcount2 = get_SV_nzcount(sv2);
	for (int k = 0; k < nzcount2; k++) {
		int nzoff2 = sv2->nzoffs[k];
		int xx = x[nzoff2 % x_len];
		int y = nzvals2_p == NULL ? int1 : nzvals2_p[k];
		int out_val = iarith(opcode, xx, y, ovflow);
		if (out_val == out_bg_val)
			continue;
		APPEND_TO_NZVALS_NZOFFS(out_val, nzoff2,
				out_nzvals, out_sv->nzoffs, out_sv->nzcount);
	}
	return;
}


/****************************************************************************
 * dArith_intSV_intSV()
 * dArith_intSV_doubleSV()
 * dArith_doubleSV_intSV()
 * dArith_doubleSV_doubleSV()
 * iArith_intSV_intSV()
 */

#define DEFINE_dArith_SV_SV_FUN(funname, Ltype, Rtype)			    \
static void dArith_ ## Ltype ## SV_ ## Rtype ## SV(int opcode,		    \
	const SparseVec *sv1, const SparseVec *sv2, SparseVec *out_sv)	    \
{									    \
	if (out_sv->len != sv1->len || out_sv->len != sv2->len)		    \
		error("SparseArray internal error in %s():\n"		    \
		      "    'sv1', 'sv2', and 'out_sv' are incompatible",    \
		      funname);						    \
	check_outRtype(get_SV_Rtype(out_sv), REALSXP, funname);		    \
	double (*darith_fun)(int, Ltype, Rtype);			    \
	darith_fun = &darith_ ## Ltype ## _ ## Rtype;			    \
	double *out_nzvals = (double *) out_sv->nzvals;			    \
	int out_nzcount = 0, k1 = 0, k2 = 0, off;			    \
	Ltype x;							    \
	Rtype y;							    \
	while (next_ ## Ltype ## SV_ ## Rtype ## SV_vals(sv1, sv2,	    \
			&k1, &k2, &off, &x, &y))			    \
	{								    \
		double out_val = darith_fun(opcode, x, y);		    \
		if (IS_BG_DOUBLE(out_val, out_sv->na_background))	    \
			continue;					    \
		APPEND_TO_NZVALS_NZOFFS(out_val, off,			    \
				out_nzvals, out_sv->nzoffs, out_nzcount);   \
	}								    \
	out_sv->nzcount = out_nzcount;					    \
	return;								    \
}

DEFINE_dArith_SV_SV_FUN("dArith_intSV_intSV", int, int)
DEFINE_dArith_SV_SV_FUN("dArith_intSV_doubleSV", int, double)
DEFINE_dArith_SV_SV_FUN("dArith_doubleSV_intSV", double, int)
DEFINE_dArith_SV_SV_FUN("dArith_doubleSV_doubleSV", double, double)

static void iArith_intSV_intSV(int opcode,
		const SparseVec *sv1, const SparseVec *sv2, SparseVec *out_sv,
		int *ovflow)
{
	if (out_sv->len != sv1->len || out_sv->len != sv2->len)
		error("SparseArray internal error in "
		      "iArith_intSV_intSV():\n"
		      "    'sv1', 'sv2', and 'out_sv' are incompatible");
	check_outRtype(out_sv->Rtype, INTSXP, "iArith_intSV_intSV");
	int *out_nzvals = (int *) out_sv->nzvals;
	int out_bg_val = out_sv->na_background ? intNA : int0;
	int out_nzcount = 0, k1 = 0, k2 = 0, off, x, y;
	while (next_intSV_intSV_vals(sv1, sv2, &k1, &k2, &off, &x, &y)) {
		int out_val = iarith(opcode, x, y, ovflow);
		if (out_val == out_bg_val)
			continue;
		APPEND_TO_NZVALS_NZOFFS(out_val, off,
				out_nzvals, out_sv->nzoffs, out_nzcount);
	}
	out_sv->nzcount = out_nzcount;
	return;
}


/****************************************************************************
 * Arith_SV_ints()
 * Arith_SV_doubles()
 * Arith_ints_SV()
 * Arith_doubles_SV()
 */

static void Arith_SV_ints(int opcode,
		const SparseVec *sv1, const int *y, int y_len,
		SparseVec *out_sv, int *ovflow)
{
	SEXPTYPE Rtype1 = get_SV_Rtype(sv1);
	switch (Rtype1) {
	    case INTSXP:
		if (opcode == DIV_OPCODE || opcode == POW_OPCODE) {
			dArith_intSV_ints(opcode, sv1, y, y_len, out_sv);
		} else {
			iArith_intSV_ints(opcode, sv1, y, y_len, out_sv,
					  ovflow);
		}
		return;
	    case REALSXP:
		dArith_doubleSV_ints(opcode, sv1, y, y_len, out_sv);
		return;
	}
	error("SparseArray internal error in Arith_SV_ints():\n"
	      "    'sv1' of type \"%s\" not supported yet",
	      type2char(Rtype1));
}

static void Arith_SV_doubles(int opcode,
		const SparseVec *sv1, const double *y, int y_len,
		SparseVec *out_sv)
{
	SEXPTYPE Rtype1 = get_SV_Rtype(sv1);
	switch (Rtype1) {
	    case INTSXP:
		dArith_intSV_doubles(opcode, sv1, y, y_len, out_sv);
		return;
	    case REALSXP:
		dArith_doubleSV_doubles(opcode, sv1, y, y_len, out_sv);
		return;
	}
	error("SparseArray internal error in Arith_SV_doubles():\n"
	      "    'sv1' of type \"%s\" not supported yet",
	      type2char(Rtype1));
}

static void Arith_ints_SV(int opcode,
		const int *x, int x_len, const SparseVec *sv2,
		SparseVec *out_sv, int *ovflow)
{
	SEXPTYPE Rtype2 = get_SV_Rtype(sv2);
	switch (Rtype2) {
	    case INTSXP:
		if (opcode == DIV_OPCODE || opcode == POW_OPCODE) {
			dArith_ints_intSV(opcode, x, x_len, sv2, out_sv);
		} else {
			iArith_ints_intSV(opcode, x, x_len, sv2, out_sv,
					  ovflow);
		}
		return;
	    case REALSXP:
		dArith_ints_doubleSV(opcode, x, x_len, sv2, out_sv);
		return;
	}
	error("SparseArray internal error in Arith_ints_SV():\n"
	      "    'sv2' of type \"%s\" not supported yet",
	      type2char(Rtype2));
}

static void Arith_doubles_SV(int opcode,
		const double *x, int x_len, const SparseVec *sv2,
		SparseVec *out_sv)
{
	SEXPTYPE Rtype2 = get_SV_Rtype(sv2);
	switch (Rtype2) {
	    case INTSXP:
		dArith_doubles_intSV(opcode, x, x_len, sv2, out_sv);
		return;
	    case REALSXP:
		dArith_doubles_doubleSV(opcode, x, x_len, sv2, out_sv);
		return;
	}
	error("SparseArray internal error in Arith_doubles_SV():\n"
	      "    'sv2' of type \"%s\" not supported yet",
	      type2char(Rtype2));
}


/****************************************************************************
 * _Arith_sv1_v2()
 * _Arith_v1_sv2()
 */

void _Arith_sv1_v2(int opcode, const SparseVec *sv1, SEXP v2, int i2,
		   SparseVec *out_sv, int *ovflow)
{
	if (out_sv->na_background != sv1->na_background)
		error("SparseArray internal error in "
		      "_Arith_sv1_v2():\n"
		      "    out_sv->na_background != sv1->na_background");
	SEXPTYPE Rtype2 = TYPEOF(v2);
	int y_len = LENGTH(v2);
	switch (Rtype2) {
	    case INTSXP: {
		const int *y = INTEGER(v2);
		if (i2 >= 0) {
			y += i2 % y_len;
			y_len = 1;
		}
		Arith_SV_ints(opcode, sv1, y, y_len, out_sv, ovflow);
		return;
	    }
	    case REALSXP: {
		const double *y = REAL(v2);
		if (i2 >= 0) {
			y += i2 % y_len;
			y_len = 1;
		}
		Arith_SV_doubles(opcode, sv1, y, y_len, out_sv);
		return;
	    }
	}
	error("SparseArray internal error in _Arith_sv1_v2():\n"
	      "    'v2' of type \"%s\" not supported yet",
	      type2char(Rtype2));
}

void _Arith_v1_sv2(int opcode, SEXP v1, const SparseVec *sv2,
		   SparseVec *out_sv, int *ovflow)
{
	if (out_sv->na_background != sv2->na_background)
		error("SparseArray internal error in "
		      "_Arith_v1_sv2():\n"
		      "    out_sv->na_background != sv2->na_background");
	SEXPTYPE Rtype1 = TYPEOF(v1);
	switch (Rtype1) {
	    case INTSXP:
		Arith_ints_SV(opcode, INTEGER(v1), LENGTH(v1), sv2,
			      out_sv, ovflow);
		return;
	    case REALSXP:
		Arith_doubles_SV(opcode, REAL(v1), LENGTH(v1), sv2,
			      out_sv);
		return;
	}
	error("SparseArray internal error in _Arith_v1_sv2():\n"
	      "    'v1' of type \"%s\" not supported yet",
	      type2char(Rtype1));
}


/****************************************************************************
 * _Arith_sv1_zero()
 * _Arith_sv1_na()
 */

/* Multiplies the vals in 'sv1' with zero. Will return 0 (i.e. no output) if
   the nonzero values in 'sv1' are finite (i.e. no NA, NaN, Inf, or -Inf).
   Note that this could also be achieved with something like:

     Arith_SV_ints(MULT_OPCODE, sv1, &int0, 1, ...);

   but mult_sv1_zero() takes a lot of shortcuts so is A LOT more efficient.
   Assumes that 'out_sv->Rtype' is equal or bigger than the type of the
   nonzero values in 'sv1'. */
static void mult_sv1_zero(const SparseVec *sv1, SparseVec *out_sv)
{
	const int *nzvals_p = get_intSV_nzvals_p(sv1);
	if (nzvals_p == NULL)  { /* lacunar SparseVec */
		out_sv->nzcount = 0;
		return;
	}
	/* regular SparseVec */
	SEXPTYPE Rtype = get_SV_Rtype(sv1);
	if (Rtype == INTSXP) {
		int nzcount, in_nzcount = get_SV_nzcount(sv1);
		if (out_sv->Rtype == INTSXP) {
			/* We only keep NAs. */
			int *out_nzvals = (int *) out_sv->nzvals;
			for (int k = nzcount = 0; k < in_nzcount; k++) {
				int x = nzvals_p[k];
				if (x == intNA) {
					out_nzvals[nzcount] = intNA;
					out_sv->nzoffs[nzcount] =
						sv1->nzoffs[k];
					nzcount++;
				}
			}
			out_sv->nzcount = nzcount;
			return;
		}
		if (out_sv->Rtype == REALSXP) {
			/* We only keep NAs. */
			double *out_nzvals = (double *) out_sv->nzvals;
			for (int k = nzcount = 0; k < in_nzcount; k++) {
				int x = nzvals_p[k];
				if (x == intNA) {
					out_nzvals[nzcount] = doubleNA;
					out_sv->nzoffs[nzcount] =
						sv1->nzoffs[k];
					nzcount++;
				}
			}
			out_sv->nzcount = nzcount;
			return;
		}
	} else if (Rtype == REALSXP) {
		if (out_sv->Rtype == REALSXP) {
			dArith_doubleSV_doubles(MULT_OPCODE, sv1, &double0, 1,
						out_sv);
			return;
		}
	}
	error("mult_sv1_zero() only supports input of "
	      "type \"integer\" or \"double\" at the moment");
}

void _Arith_sv1_zero(int opcode, const SparseVec *sv1, SEXPTYPE Rtype2,
		SparseVec *out_sv)
{
	if (out_sv->na_background != sv1->na_background)
		error("SparseArray internal error in "
		      "_Arith_sv1_zero():\n"
		      "    out_sv->na_background != sv1->na_background");
	if (!sv1->na_background && opcode == MULT_OPCODE) {
		mult_sv1_zero(sv1, out_sv);
		return;
	}
	switch (Rtype2) {
	    case INTSXP: {
		int ovflow = 0;
		Arith_SV_ints(opcode, sv1, &int0, 1, out_sv, &ovflow);
		if (ovflow)
			error("SparseArray internal error in "
			      "_Arith_sv1_zero():\n"
			      "    unexpected integer overflow");
		return;
	    }
	    case REALSXP:
		Arith_SV_doubles(opcode, sv1, &double0, 1, out_sv);
		return;
	}
	error("SparseArray internal error in _Arith_sv1_zero():\n"
	      "    zero of type \"%s\" not supported yet",
	      type2char(Rtype2));
}

void _Arith_sv1_na(int opcode, const SparseVec *sv1, SEXPTYPE Rtype2,
		SparseVec *out_sv)
{
	if (!out_sv->na_background)
		error("SparseArray internal error in "
		      "_Arith_sv1_na():\n"
		      "    'out_sv->na_background' is FALSE");
	switch (Rtype2) {
	    case INTSXP: {
		int ovflow = 0;
		Arith_SV_ints(opcode, sv1, &intNA, 1, out_sv, &ovflow);
		if (ovflow)
			error("SparseArray internal error in "
			      "_Arith_sv1_na():\n"
			      "    unexpected integer overflow");
		return;
	    }
	    case REALSXP:
		Arith_SV_doubles(opcode, sv1, &doubleNA, 1, out_sv);
		return;
	}
	error("SparseArray internal error in _Arith_sv1_na():\n"
	      "    NA of type \"%s\" not supported yet",
	      type2char(Rtype2));
}


/****************************************************************************
 * _Arith_zero_sv2()
 * _Arith_na_sv2()
 */

void _Arith_zero_sv2(int opcode, SEXPTYPE Rtype1, const SparseVec *sv2,
		SparseVec *out_sv)
{
	if (out_sv->na_background != sv2->na_background)
		error("SparseArray internal error in "
		      "_Arith_zero_sv2():\n"
		      "    out_sv->na_background != sv2->na_background");
	switch (Rtype1) {
	    case INTSXP: {
		int ovflow = 0;
		Arith_ints_SV(opcode, &int0, 1, sv2, out_sv, &ovflow);
		if (ovflow)
			error("SparseArray internal error in "
			      "_Arith_zero_sv2():\n"
			      "    unexpected integer overflow");
		return;
	    }
	    case REALSXP:
		Arith_doubles_SV(opcode, &double0, 1, sv2, out_sv);
		return;
	}
	error("SparseArray internal error in _Arith_zero_sv2():\n"
	      "    zero of type \"%s\" not supported yet",
	      type2char(Rtype1));
}

void _Arith_na_sv2(int opcode, SEXPTYPE Rtype1, const SparseVec *sv2,
		SparseVec *out_sv)
{
	if (!out_sv->na_background)
		error("SparseArray internal error in "
		      "_Arith_na_sv2():\n"
		      "    'out_sv->na_background' is FALSE");
	switch (Rtype1) {
	    case INTSXP: {
		int ovflow = 0;
		Arith_ints_SV(opcode, &intNA, 1, sv2, out_sv, &ovflow);
		if (ovflow)
			error("SparseArray internal error in "
			      "_Arith_na_sv2():\n"
			      "    unexpected integer overflow");
		return;
	    }
	    case REALSXP:
		Arith_doubles_SV(opcode, &doubleNA, 1, sv2, out_sv);
		return;
	}
	error("SparseArray internal error in _Arith_na_sv2():\n"
	      "    NA of type \"%s\" not supported yet",
	      type2char(Rtype1));
}


/****************************************************************************
 * _Arith_sv1_sv2()
 */

void _Arith_sv1_sv2(int opcode, const SparseVec *sv1, const SparseVec *sv2,
		SparseVec *out_sv, int *ovflow)
{
	SEXPTYPE Rtype1 = get_SV_Rtype(sv1);
	SEXPTYPE Rtype2 = get_SV_Rtype(sv2);
	switch (Rtype1) {
	    case INTSXP:
		if (Rtype2 == INTSXP) {
			if (opcode == DIV_OPCODE || opcode == POW_OPCODE) {
				dArith_intSV_intSV(opcode, sv1, sv2, out_sv);
			} else {
				iArith_intSV_intSV(opcode, sv1, sv2, out_sv,
						   ovflow);
			}
			return;
		}
		if (Rtype2 == REALSXP) {
			dArith_intSV_doubleSV(opcode, sv1, sv2, out_sv);
			return;
		}
		break;
	    case REALSXP:
		if (Rtype2 == INTSXP) {
			dArith_doubleSV_intSV(opcode, sv1, sv2, out_sv);
			return;
		}
		if (Rtype2 == REALSXP) {
			dArith_doubleSV_doubleSV(opcode, sv1, sv2, out_sv);
			return;
		}
		break;
	}
	error("_Arith_sv1_sv2() only supports input of "
	      "type \"integer\" or \"double\" at the moment");
	return;
}

