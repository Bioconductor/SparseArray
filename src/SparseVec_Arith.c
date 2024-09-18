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
static inline int Arith_int(int opcode, int x, int y, int *ovflow)
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

static inline double Arith_double(int opcode, double x, double y)
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

static int Arith_intSV_int(int opcode,
		const SparseVec *sv1, int y,
		int *out_nzvals, int *out_nzoffs, int out_has_NAbg,
		int *ovflow)
{
	int out_background = out_has_NAbg ? intNA : int0;
	const int *nzvals1_p = get_intSV_nzvals_p(sv1);
	if (nzvals1_p == NULL) {  /* lacunar SparseVec */
		int out_val = Arith_int(opcode, int1, y, ovflow);
		if (out_val == out_background)
			return 0;
		out_nzvals[0] = out_val;
		return PROPAGATE_NZOFFS;
	}
	/* regular SparseVec */
	int nzcount1 = get_SV_nzcount(sv1);
	int out_nzcount = 0;
	for (int k = 0; k < nzcount1; k++) {
		int out_val = Arith_int(opcode, nzvals1_p[k], y, ovflow);
		if (out_val == out_background)
			continue;
		out_nzvals[out_nzcount] = out_val;
		out_nzoffs[out_nzcount] = sv1->nzoffs[k];
		out_nzcount++;
	}
	return out_nzcount;
}

static int Arith_int_intSV(int opcode,
		int x, const SparseVec *sv2,
		int *out_nzvals, int *out_nzoffs, int out_has_NAbg,
		int *ovflow)
{
	int out_background = out_has_NAbg ? intNA : int0;
	const int *nzvals2_p = get_intSV_nzvals_p(sv2);
	if (nzvals2_p == NULL) {  /* lacunar SparseVec */
		int out_val = Arith_int(opcode, x, int1, ovflow);
		if (out_val == out_background)
			return 0;
		out_nzvals[0] = out_val;
		return PROPAGATE_NZOFFS;
	}
	/* regular SparseVec */
	int nzcount2 = get_SV_nzcount(sv2);
	int out_nzcount = 0;
	for (int k = 0; k < nzcount2; k++) {
		int out_val = Arith_int(opcode, x, nzvals2_p[k], ovflow);
		if (out_val == out_background)
			continue;
		out_nzvals[out_nzcount] = out_val;
		out_nzoffs[out_nzcount] = sv2->nzoffs[k];
		out_nzcount++;
	}
	return out_nzcount;
}

static void div_intSV_intSV(const SparseVec *sv1, const SparseVec *sv2,
		SparseVec *out_sv)
{
	check_outRtype(out_sv->Rtype, REALSXP, "div_intSV_intSV");
	double *out_nzvals = (double *) out_sv->nzvals;
	int out_nzcount = 0, k1 = 0, k2 = 0, off, x, y;
	while (next_int_int_vals(sv1, sv2, &k1, &k2, &off, &x, &y)) {
		double xx = x == intNA ? doubleNA : (double) x;
		double yy = y == intNA ? doubleNA : (double) y;
		double out_val = xx / yy;
		if (IS_BACKGROUND_VAL(out_val, out_sv->na_background))
			continue;
		out_nzvals[out_nzcount] = out_val;
		out_sv->nzoffs[out_nzcount] = off;
		out_nzcount++;
	}
	out_sv->nzcount = out_nzcount;
	return;
}

static void pow_intSV_intSV(const SparseVec *sv1, const SparseVec *sv2,
		SparseVec *out_sv)
{
	check_outRtype(out_sv->Rtype, REALSXP, "pow_intSV_intSV");
	double *out_nzvals = (double *) out_sv->nzvals;
	int out_nzcount = 0, k1 = 0, k2 = 0, off, x, y;
	while (next_int_int_vals(sv1, sv2, &k1, &k2, &off, &x, &y)) {
		double xx = x == intNA ? doubleNA : (double) x;
		double yy = y == intNA ? doubleNA : (double) y;
		double out_val = Rpow_double(xx, yy);
		if (IS_BACKGROUND_VAL(out_val, out_sv->na_background))
			continue;
		out_nzvals[out_nzcount] = out_val;
		out_sv->nzoffs[out_nzcount] = off;
		out_nzcount++;
	}
	out_sv->nzcount = out_nzcount;
	return;
}

static void Arith_intSV_intSV(int opcode,
		const SparseVec *sv1, const SparseVec *sv2,
		SparseVec *out_sv, int *ovflow)
{
	if (out_sv->len != sv1->len || out_sv->len != sv2->len)
		error("SparseArray internal error in "
		      "Arith_intSV_intSV():\n"
		      "    sv1, sv2, out_sv are incompatible");
	if (opcode == DIV_OPCODE) {
		div_intSV_intSV(sv1, sv2, out_sv);
		return;
	}
	if (opcode == POW_OPCODE) {
		pow_intSV_intSV(sv1, sv2, out_sv);
		return;
	}
	check_outRtype(out_sv->Rtype, INTSXP, "Arith_intSV_intSV");
	int *out_nzvals = (int *) out_sv->nzvals;
	int out_background = out_sv->na_background ? intNA : int0;
	int out_nzcount = 0, k1 = 0, k2 = 0, off, x, y;
	while (next_int_int_vals(sv1, sv2, &k1, &k2, &off, &x, &y)) {
		int out_val = Arith_int(opcode, x, y, ovflow);
		if (out_val == out_background)
			continue;
		out_nzvals[out_nzcount] = out_val;
		out_sv->nzoffs[out_nzcount] = off;
		out_nzcount++;
	}
	out_sv->nzcount = out_nzcount;
	return;
}

static int Arith_intSV_double(int opcode,
		const SparseVec *sv1, double y,
		double *out_nzvals, int *out_nzoffs, int out_has_NAbg)
{
	const int *nzvals1_p = get_intSV_nzvals_p(sv1);
	if (nzvals1_p == NULL) {  /* lacunar SparseVec */
		double out_val = Arith_double(opcode, double1, y);
		if (IS_BACKGROUND_VAL(out_val, out_has_NAbg))
			return 0;
		out_nzvals[0] = out_val;
		return PROPAGATE_NZOFFS;
	}
	/* regular SparseVec */
	int nzcount1 = get_SV_nzcount(sv1);
	int out_nzcount = 0;
	for (int k = 0; k < nzcount1; k++) {
		int x = nzvals1_p[k];
		double xx = x == intNA ? doubleNA : (double) x;
		double out_val = Arith_double(opcode, xx, y);
		if (IS_BACKGROUND_VAL(out_val, out_has_NAbg))
			continue;
		out_nzvals[out_nzcount] = out_val;
		out_nzoffs[out_nzcount] = sv1->nzoffs[k];
		out_nzcount++;
	}
	return out_nzcount;
}

static int Arith_double_intSV(int opcode,
		double x, const SparseVec *sv2,
		double *out_nzvals, int *out_nzoffs, int out_has_NAbg)
{
	const int *nzvals2_p = get_intSV_nzvals_p(sv2);
	if (nzvals2_p == NULL) {  /* lacunar SparseVec */
		double out_val = Arith_double(opcode, x, double1);
		if (IS_BACKGROUND_VAL(out_val, out_has_NAbg))
			return 0;
		out_nzvals[0] = out_val;
		return PROPAGATE_NZOFFS;
	}
	/* regular SparseVec */
	int nzcount2 = get_SV_nzcount(sv2);
	int out_nzcount = 0;
	for (int k = 0; k < nzcount2; k++) {
		int y = nzvals2_p[k];
		double yy = y == intNA ? doubleNA : (double) y;
		double out_val = Arith_double(opcode, x, yy);
		if (IS_BACKGROUND_VAL(out_val, out_has_NAbg))
			continue;
		out_nzvals[out_nzcount] = out_val;
		out_nzoffs[out_nzcount] = sv2->nzoffs[k];
		out_nzcount++;
	}
	return out_nzcount;
}

static void Arith_intSV_doubleSV(int opcode,
		const SparseVec *sv1, const SparseVec *sv2,
		SparseVec *out_sv)
{
	if (out_sv->len != sv1->len || out_sv->len != sv2->len)
		error("SparseArray internal error in "
		      "Arith_intSV_doubleSV():\n"
		      "    sv1, sv2, out_sv are incompatible");
	check_outRtype(out_sv->Rtype, REALSXP, "Arith_intSV_doubleSV");
	double *out_nzvals = (double *) out_sv->nzvals;
	int out_nzcount = 0, k1 = 0, k2 = 0, off, x;
	double y;
	while (next_int_double_vals(sv1, sv2, &k1, &k2, &off, &x, &y)) {
		double xx = x == intNA ? doubleNA : (double) x;
		double out_val = Arith_double(opcode, xx, y);
		if (IS_BACKGROUND_VAL(out_val, out_sv->na_background))
			continue;
		out_nzvals[out_nzcount] = out_val;
		out_sv->nzoffs[out_nzcount] = off;
		out_nzcount++;
	}
	out_sv->nzcount = out_nzcount;
	return;
}

static void Arith_doubleSV_intSV(int opcode,
		const SparseVec *sv1, const SparseVec *sv2,
		SparseVec *out_sv)
{
	if (out_sv->len != sv1->len || out_sv->len != sv2->len)
		error("SparseArray internal error in "
		      "Arith_doubleSV_intSV():\n"
		      "    sv1, sv2, out_sv are incompatible");
	check_outRtype(out_sv->Rtype, REALSXP, "Arith_doubleSV_intSV");
	double *out_nzvals = (double *) out_sv->nzvals;
	int out_nzcount = 0, k1 = 0, k2 = 0, off, y;
	double x;
	while (next_double_int_vals(sv1, sv2, &k1, &k2, &off, &x, &y)) {
		double yy = y == intNA ? doubleNA : (double) y;
		double out_val = Arith_double(opcode, x, yy);
		if (IS_BACKGROUND_VAL(out_val, out_sv->na_background))
			continue;
		out_nzvals[out_nzcount] = out_val;
		out_sv->nzoffs[out_nzcount] = off;
		out_nzcount++;
	}
	out_sv->nzcount = out_nzcount;
	return;
}

static int Arith_doubleSV_double(int opcode,
		const SparseVec *sv1, double y,
		double *out_nzvals, int *out_nzoffs, int out_has_NAbg)
{
	const double *nzvals1_p = get_doubleSV_nzvals_p(sv1);
	if (nzvals1_p == NULL) {  /* lacunar SparseVec */
		double out_val = Arith_double(opcode, double1, y);
		if (IS_BACKGROUND_VAL(out_val, out_has_NAbg))
			return 0;
		out_nzvals[0] = out_val;
		return PROPAGATE_NZOFFS;
	}
	/* regular SparseVec */
	int nzcount1 = get_SV_nzcount(sv1);
	int out_nzcount = 0;
	for (int k = 0; k < nzcount1; k++) {
		double out_val = Arith_double(opcode, nzvals1_p[k], y);
		if (IS_BACKGROUND_VAL(out_val, out_has_NAbg))
			continue;
		out_nzvals[out_nzcount] = out_val;
		out_nzoffs[out_nzcount] = sv1->nzoffs[k];
		out_nzcount++;
	}
	return out_nzcount;
}

static int Arith_double_doubleSV(int opcode,
		double x, const SparseVec *sv2,
		double *out_nzvals, int *out_nzoffs, int out_has_NAbg)
{
	const double *nzvals2_p = get_doubleSV_nzvals_p(sv2);
	if (nzvals2_p == NULL) {  /* lacunar SparseVec */
		double out_val = Arith_double(opcode, x, double1);
		if (IS_BACKGROUND_VAL(out_val, out_has_NAbg))
			return 0;
		out_nzvals[0] = out_val;
		return PROPAGATE_NZOFFS;
	}
	/* regular SparseVec */
	int nzcount2 = get_SV_nzcount(sv2);
	int out_nzcount = 0;
	for (int k = 0; k < nzcount2; k++) {
		double out_val = Arith_double(opcode, x, nzvals2_p[k]);
		if (IS_BACKGROUND_VAL(out_val, out_has_NAbg))
			continue;
		out_nzvals[out_nzcount] = out_val;
		out_nzoffs[out_nzcount] = sv2->nzoffs[k];
		out_nzcount++;
	}
	return out_nzcount;
}

static void Arith_doubleSV_doubleSV(int opcode,
		const SparseVec *sv1, const SparseVec *sv2,
		SparseVec *out_sv)
{
	if (out_sv->len != sv1->len || out_sv->len != sv2->len)
		error("SparseArray internal error in "
		      "Arith_doubleSV_doubleSV():\n"
		      "    sv1, sv2, out_sv are incompatible");
	check_outRtype(out_sv->Rtype, REALSXP, "Arith_doubleSV_doubleSV");
	double *out_nzvals = (double *) out_sv->nzvals;
	int out_nzcount = 0, k1 = 0, k2 = 0, off;
	double x, y;
	while (next_double_double_vals(sv1, sv2, &k1, &k2, &off, &x, &y)) {
		double out_val = Arith_double(opcode, x, y);
		if (IS_BACKGROUND_VAL(out_val, out_sv->na_background))
			continue;
		out_nzvals[out_nzcount] = out_val;
		out_sv->nzoffs[out_nzcount] = off;
		out_nzcount++;
	}
	out_sv->nzcount = out_nzcount;
	return;
}


/****************************************************************************
 * Arith_sv1_int()
 * Arith_int_sv2()
 * Arith_sv1_double()
 * Arith_double_sv2()
 */

static void Arith_sv1_int(int opcode, const SparseVec *sv1, int y,
		SparseVec *out_sv, int *ovflow)
{
	SEXPTYPE Rtype1 = get_SV_Rtype(sv1);
	if (Rtype1 == INTSXP && opcode != DIV_OPCODE && opcode != POW_OPCODE) {
		check_outRtype(out_sv->Rtype, INTSXP, "Arith_sv1_int");
		out_sv->nzcount =
			Arith_intSV_int(opcode, sv1, y,
				(int *) out_sv->nzvals, out_sv->nzoffs,
				out_sv->na_background, ovflow);
		return;
	}
	check_outRtype(out_sv->Rtype, REALSXP, "Arith_sv1_int");
	double yy = y == intNA ? doubleNA : (double) y;
	switch (Rtype1) {
	    case INTSXP:
		out_sv->nzcount =
			Arith_intSV_double(opcode, sv1, yy,
				(double *) out_sv->nzvals, out_sv->nzoffs,
				out_sv->na_background);
		return;
	    case REALSXP:
		out_sv->nzcount =
			Arith_doubleSV_double(opcode, sv1, yy,
				(double *) out_sv->nzvals, out_sv->nzoffs,
				out_sv->na_background);
		return;
	}
	error("SparseArray internal error in Arith_sv1_int():\n"
	      "    'sv1' of type \"%s\" not supported yet",
	      type2char(Rtype1));
}

static void Arith_int_sv2(int opcode, int x, const SparseVec *sv2,
		SparseVec *out_sv, int *ovflow)
{
	SEXPTYPE Rtype2 = get_SV_Rtype(sv2);
	if (Rtype2 == INTSXP && opcode != DIV_OPCODE && opcode != POW_OPCODE) {
		check_outRtype(out_sv->Rtype, INTSXP, "Arith_int_sv2");
		out_sv->nzcount =
			Arith_int_intSV(opcode, x, sv2,
				(int *) out_sv->nzvals, out_sv->nzoffs,
				out_sv->na_background, ovflow);
		return;
	}
	check_outRtype(out_sv->Rtype, REALSXP, "Arith_int_sv2");
	double xx = x == intNA ? doubleNA : (double) x;
	switch (Rtype2) {
	    case INTSXP:
		out_sv->nzcount =
			Arith_double_intSV(opcode, xx, sv2,
				(double *) out_sv->nzvals, out_sv->nzoffs,
				out_sv->na_background);
		return;
	    case REALSXP:
		out_sv->nzcount =
			Arith_double_doubleSV(opcode, xx, sv2,
				(double *) out_sv->nzvals, out_sv->nzoffs,
				out_sv->na_background);
		return;
	}
	error("SparseArray internal error in Arith_int_sv2():\n"
	      "    'sv2' of type \"%s\" not supported yet",
	      type2char(Rtype2));
}

static void Arith_sv1_double(int opcode, const SparseVec *sv1, double y,
		SparseVec *out_sv)
{
	check_outRtype(out_sv->Rtype, REALSXP, "Arith_sv1_double");
	SEXPTYPE Rtype1 = get_SV_Rtype(sv1);
	switch (Rtype1) {
	    case INTSXP:
		out_sv->nzcount =
			Arith_intSV_double(opcode, sv1, y,
				out_sv->nzvals, out_sv->nzoffs,
				out_sv->na_background);
		return;
	    case REALSXP:
		out_sv->nzcount =
			Arith_doubleSV_double(opcode, sv1, y,
				out_sv->nzvals, out_sv->nzoffs,
				out_sv->na_background);
		return;
	}
	error("SparseArray internal error in Arith_sv1_double():\n"
	      "    'sv1' of type \"%s\" not supported yet",
	      type2char(Rtype1));
}

static void Arith_double_sv2(int opcode, double x, const SparseVec *sv2,
		SparseVec *out_sv)
{
	check_outRtype(out_sv->Rtype, REALSXP, "Arith_double_sv2");
	SEXPTYPE Rtype2 = get_SV_Rtype(sv2);
	switch (Rtype2) {
	    case INTSXP:
		out_sv->nzcount =
			Arith_double_intSV(opcode, x, sv2,
				out_sv->nzvals, out_sv->nzoffs,
				out_sv->na_background);
		return;
	    case REALSXP:
		out_sv->nzcount =
			Arith_double_doubleSV(opcode, x, sv2,
				out_sv->nzvals, out_sv->nzoffs,
				out_sv->na_background);
		return;
	}
	error("SparseArray internal error in Arith_double_sv2():\n"
	      "    'sv2' of type \"%s\" not supported yet",
	      type2char(Rtype2));
}


/****************************************************************************
 * _Arith_sv1_scalar()
 * _Arith_scalar_sv2()
 */

/* 'scalar' is assumed to be an atomic vector of length 1.
   This is NOT checked! */
void _Arith_sv1_scalar(int opcode, const SparseVec *sv1, SEXP scalar,
		SparseVec *out_sv, int *ovflow)
{
	if (out_sv->na_background != sv1->na_background)
		error("SparseArray internal error in "
		      "_Arith_sv1_scalar():\n"
		      "    out_sv->na_background != sv1->na_background");
	SEXPTYPE Rtype2 = TYPEOF(scalar);
	switch (Rtype2) {
	    case INTSXP:
		Arith_sv1_int(opcode, sv1, INTEGER(scalar)[0], out_sv, ovflow);
		return;
	    case REALSXP:
		Arith_sv1_double(opcode, sv1, REAL(scalar)[0], out_sv);
		return;
	}
	error("SparseArray internal error in _Arith_sv1_scalar():\n"
	      "    'scalar' of type \"%s\" not supported yet",
	      type2char(Rtype2));
}

/* 'scalar' is assumed to be an atomic vector of length 1.
   This is NOT checked! */
void _Arith_scalar_sv2(int opcode, SEXP scalar, const SparseVec *sv2,
		SparseVec *out_sv, int *ovflow)
{
	if (out_sv->na_background != sv2->na_background)
		error("SparseArray internal error in "
		      "_Arith_scalar_sv2():\n"
		      "    out_sv->na_background != sv2->na_background");
	SEXPTYPE Rtype1 = TYPEOF(scalar);
	switch (Rtype1) {
	    case INTSXP:
		Arith_int_sv2(opcode, INTEGER(scalar)[0], sv2, out_sv, ovflow);
		return;
	    case REALSXP:
		Arith_double_sv2(opcode, REAL(scalar)[0], sv2, out_sv);
		return;
	}
	error("SparseArray internal error in _Arith_scalar_sv2():\n"
	      "    'scalar' of type \"%s\" not supported yet",
	      type2char(Rtype1));
}


/****************************************************************************
 * _Arith_sv1_zero()
 * _Arith_sv1_na()
 */

/* Multiplies the vals in 'sv' with zero. Will return 0 (i.e. no output) if
   the nonzero values in 'sv' are finite (i.e. no NA, NaN, Inf, or -Inf).
   Note that this could simply be achieved by calling:

     _Arith_sv1_scalar(MULT_OPCODE, sv1, 0, ...)

   but mult_sv1_zero() takes a lot of shortcuts so is A LOT more efficient.
   Assumes that 'outRtype' is equal or bigger than the type of the
   nonzero values in 'sv1'. */
static int mult_sv1_zero(const SparseVec *sv1,
		SEXPTYPE outRtype, void *out_nzvals, int *out_nzoffs)
{
	const int *nzvals_p = get_intSV_nzvals_p(sv1);
	if (nzvals_p == NULL)  /* lacunar SparseVec */
		return 0;
	/* regular SparseVec */
	int nzcount = NZCOUNT_IS_NOT_SET;
	SEXPTYPE Rtype = get_SV_Rtype(sv1);
	if (Rtype == INTSXP) {
		int in_nzcount = get_SV_nzcount(sv1);
		if (outRtype == INTSXP) {
			/* We only keep NAs. */
			int *out_nzvals_p = (int *) out_nzvals;
			for (int k = nzcount = 0; k < in_nzcount; k++) {
				int x = nzvals_p[k];
				if (x == NA_INTEGER) {
					out_nzvals_p[nzcount] = NA_INTEGER;
					out_nzoffs[nzcount] = sv1->nzoffs[k];
					nzcount++;
				}
			}
		} else if (outRtype == REALSXP) {
			/* We only keep NAs. */
			double *out_nzvals_p = (double *) out_nzvals;
			for (int k = nzcount = 0; k < in_nzcount; k++) {
				int x = nzvals_p[k];
				if (x == NA_INTEGER) {
					out_nzvals_p[nzcount] = NA_REAL;
					out_nzoffs[nzcount] = sv1->nzoffs[k];
					nzcount++;
				}
			}
		}
	} else if (Rtype == REALSXP) {
		if (outRtype == REALSXP) {
			nzcount = Arith_doubleSV_double(MULT_OPCODE,
					sv1, 0.0,
					(double *) out_nzvals, out_nzoffs,
					sv1->na_background);
		}
	}
	if (nzcount == NZCOUNT_IS_NOT_SET)
		error("mult_sv1_zero() only supports input of "
		      "type \"integer\" or \"double\" at the moment");
	return nzcount;
}

void _Arith_sv1_zero(int opcode, const SparseVec *sv1, SEXPTYPE Rtype2,
		SparseVec *out_sv)
{
	if (out_sv->na_background != sv1->na_background)
		error("SparseArray internal error in "
		      "_Arith_sv1_zero():\n"
		      "    out_sv->na_background != sv1->na_background");
	if (!sv1->na_background && opcode == MULT_OPCODE) {
		out_sv->nzcount =
			mult_sv1_zero(sv1, out_sv->Rtype,
				      out_sv->nzvals, out_sv->nzoffs);
		return;
	}
	switch (Rtype2) {
	    case INTSXP: {
		int ovflow = 0;
		Arith_sv1_int(opcode, sv1, int0, out_sv, &ovflow);
		if (ovflow)
			error("SparseArray internal error in "
			      "_Arith_sv1_zero():\n"
			      "    unexpected integer overflow");
		return;
	    }
	    case REALSXP:
		Arith_sv1_double(opcode, sv1, double0, out_sv);
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
		Arith_sv1_int(opcode, sv1, intNA, out_sv, &ovflow);
		if (ovflow)
			error("SparseArray internal error in "
			      "_Arith_sv1_na():\n"
			      "    unexpected integer overflow");
		return;
	    }
	    case REALSXP:
		Arith_sv1_double(opcode, sv1, doubleNA, out_sv);
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
		Arith_int_sv2(opcode, int0, sv2, out_sv, &ovflow);
		if (ovflow)
			error("SparseArray internal error in "
			      "_Arith_zero_sv2():\n"
			      "    unexpected integer overflow");
		return;
	    }
	    case REALSXP:
		Arith_double_sv2(opcode, double0, sv2, out_sv);
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
		Arith_int_sv2(opcode, intNA, sv2, out_sv, &ovflow);
		if (ovflow)
			error("SparseArray internal error in "
			      "_Arith_na_sv2():\n"
			      "    unexpected integer overflow");
		return;
	    }
	    case REALSXP:
		Arith_double_sv2(opcode, doubleNA, sv2, out_sv);
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
	if (Rtype1 == INTSXP) {
		if (Rtype2 == INTSXP) {
			Arith_intSV_intSV(opcode, sv1, sv2, out_sv, ovflow);
			return;
		}
		if (Rtype2 == REALSXP) {
			Arith_intSV_doubleSV(opcode, sv1, sv2, out_sv);
			return;
		}
	} else if (Rtype1 == REALSXP) {
		if (Rtype2 == INTSXP) {
			Arith_doubleSV_intSV(opcode, sv1, sv2, out_sv);
			return;
		}
		if (Rtype2 == REALSXP) {
			Arith_doubleSV_doubleSV(opcode, sv1, sv2, out_sv);
			return;
		}
	}
	error("_Arith_sv1_sv2() only supports input of "
	      "type \"integer\" or \"double\" at the moment");
	return;
}

