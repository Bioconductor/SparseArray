/****************************************************************************
 *                   'Arith' operations on sparse vectors                   *
 ****************************************************************************/
#include "SparseVec_Arith.h"

//#include "Rvector_utils.h"
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
	if ((x < 0 && y == R_PosInf) || (x == R_NegInf && y != trunc(y)))
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

static int Arith_intSV_int(int opcode,
		const SparseVec *sv1, int y,
		int *out_nzoffs, int *out_nzvals, int *ovflow)
{
	const int *nzvals1 = get_intSV_nzvals(sv1);
	int nzcount1 = get_SV_nzcount(sv1);
	int out_nzcount = 0;
	for (int k = 0; k < nzcount1; k++) {
		int v = Arith_int(opcode, nzvals1[k], y, ovflow);
		if (v != 0) {
			out_nzoffs[out_nzcount] = sv1->nzoffs[k];
			out_nzvals[out_nzcount] = v;
			out_nzcount++;
		}
	}
	return out_nzcount;
}

static int Arith_intSV_intSV(int opcode,
		const SparseVec *sv1, const SparseVec *sv2,
		int *out_nzoffs, int *out_nzvals, int *ovflow)
{
	int k1, k2, off, x, y;

	int out_nzcount = 0;
	k1 = k2 = 0;
	while (next_nzvals_int_int(sv1, sv2,
				   &k1, &k2, &off, &x, &y))
	{
		int v = Arith_int(opcode, x, y, ovflow);
		if (v != 0) {
			out_nzoffs[out_nzcount] = off;
			out_nzvals[out_nzcount] = v;
			out_nzcount++;
		}
	}
	return out_nzcount;
}

static int Arith_intSV_double(int opcode,
                const SparseVec *sv1, double y,
		int *out_nzoffs, double *out_nzvals)
{
	const int *nzvals1 = get_intSV_nzvals(sv1);
	int nzcount1 = get_SV_nzcount(sv1);
	int out_nzcount = 0;
	for (int k = 0; k < nzcount1; k++) {
		double v;
		int x = nzvals1[k];
		if (x == NA_INTEGER) {
			v = NA_REAL;
		} else {
			v = Arith_double(opcode, (double) x, y);
		}
		if (v != 0) {
			out_nzoffs[out_nzcount] = sv1->nzoffs[k];
			out_nzvals[out_nzcount] = v;
			out_nzcount++;
		}
	}
	return out_nzcount;
}

static int Arith_intSV_doubleSV(int opcode,
		const SparseVec *sv1, const SparseVec *sv2,
		int *out_nzoffs, double *out_nzvals)
{
	int k1, k2, off, x;
	double y;

	int out_nzcount = 0;
	k1 = k2 = 0;
	while (next_nzvals_int_double(sv1, sv2,
				      &k1, &k2, &off, &x, &y))
	{
		double v;
		if (x == NA_INTEGER) {
			v = NA_REAL;
		} else {
			v = Arith_double(opcode, (double) x, y);
		}
		if (v != 0.0) {
			out_nzoffs[out_nzcount] = off;
			out_nzvals[out_nzcount] = v;
			out_nzcount++;
		}
	}
	return out_nzcount;
}

static int Arith_doubleSV_intSV(int opcode,
		const SparseVec *sv1, const SparseVec *sv2,
		int *out_nzoffs, double *out_nzvals)
{
	int k1, k2, off, y;
	double x;

	int out_nzcount = 0;
	k1 = k2 = 0;
	while (next_nzvals_double_int(sv1, sv2,
				      &k1, &k2, &off, &x, &y))
	{
		double v;
		if (y == NA_INTEGER) {
			v = NA_REAL;
		} else {
			v = Arith_double(opcode, x, (double) y);
		}
		if (v != 0.0) {
			out_nzoffs[out_nzcount] = off;
			out_nzvals[out_nzcount] = v;
			out_nzcount++;
		}
	}
	return out_nzcount;
}

static int Arith_doubleSV_double(int opcode,
		const SparseVec *sv1, double y,
		int *out_nzoffs, double *out_nzvals)
{
	const double *nzvals1 = get_doubleSV_nzvals(sv1);
	int nzcount1 = get_SV_nzcount(sv1);
	int out_nzcount = 0;
	for (int k = 0; k < nzcount1; k++) {
		double v = Arith_double(opcode, nzvals1[k], y);
		if (v != 0.0) {
			out_nzoffs[out_nzcount] = sv1->nzoffs[k];
			out_nzvals[out_nzcount] = v;
			out_nzcount++;
		}
	}
	return out_nzcount;
}

static int Arith_doubleSV_doubleSV(int opcode,
		const SparseVec *sv1, const SparseVec *sv2,
		int *out_nzoffs, double *out_nzvals)
{
	int k1, k2, off;
	double x, y;

	int out_nzcount = 0;
	k1 = k2 = 0;
	while (next_nzvals_double_double(sv1, sv2,
					 &k1, &k2, &off, &x, &y))
	{
		double v = Arith_double(opcode, x, y);
		if (v != 0.0) {
			out_nzoffs[out_nzcount] = off;
			out_nzvals[out_nzcount] = v;
			out_nzcount++;
		}
	}
	return out_nzcount;
}

/* 'scalar' is assumed to be an atomic vector of length 1.
   This is NOT checked! */
int _Arith_sv1_scalar(int opcode, const SparseVec *sv1, SEXP scalar,
		SEXPTYPE expected_outRtype,
		int *out_nzoffs, void *out_nzvals, int *ovflow)
{
	SEXPTYPE effective_outRtype = REALSXP;
	int nzcount = -1;
	SEXPTYPE Rtype1 = get_SV_Rtype(sv1);
	if (Rtype1 == INTSXP) {
		if (TYPEOF(scalar) == INTSXP) {
			effective_outRtype = INTSXP;
			nzcount = Arith_intSV_int(opcode,
				sv1, INTEGER(scalar)[0],
				out_nzoffs, (int *) out_nzvals, ovflow);
		} else if (TYPEOF(scalar) == REALSXP) {
			nzcount = Arith_intSV_double(opcode,
				sv1, REAL(scalar)[0],
				out_nzoffs, (double *) out_nzvals);
		}
	} else if (Rtype1 == REALSXP) {
		if (TYPEOF(scalar) == REALSXP) {
			nzcount = Arith_doubleSV_double(opcode,
				sv1, REAL(scalar)[0],
				out_nzoffs, (double *) out_nzvals);
		}
	}
	if (nzcount == -1)
		error("_Arith_sv1_scalar() only supports input of "
		      "type \"integer\" or \"double\" at the moment");
	if (expected_outRtype != effective_outRtype)
		error("SparseArray internal error in "
		      "_Arith_sv1_scalar():\n"
		      "    expected_outRtype != effective_outRtype");
	return nzcount;
}

/* Multiplies the vals in 'sv' with zero. Will return 0 (i.e. no output) if
   the nonzero values in 'sv' are finite (i.e. no NA, NaN, Inf, or -Inf).
   Note that this could simply be achieved by calling:

     _Arith_sv1_scalar(MULT_OPCODE, sv, 0, ...)

   but _mult_SV_zero() takes a lot of shortcuts so is A LOT more efficient.
   Assumes that 'outRtype' is equal or bigger than the type of the
   nonzero values in 'sv'. */
int _mult_SV_zero(const SparseVec *sv,
		SEXPTYPE outRtype, int *out_nzoffs, void *out_nzvals)
{
	int nzcount = -1;
	SEXPTYPE Rtype = get_SV_Rtype(sv);
	if (Rtype == INTSXP) {
		const int *nzvals = get_intSV_nzvals(sv);
		int in_nzcount = get_SV_nzcount(sv);
		if (outRtype == INTSXP) {
			/* We only keep NAs. */
			int *out_nzvals_p = (int *) out_nzvals;
			for (int k = nzcount = 0; k < in_nzcount; k++) {
				int x = nzvals[k];
				if (x == NA_INTEGER) {
					out_nzoffs[nzcount] = sv->nzoffs[k];
					out_nzvals_p[nzcount] = NA_INTEGER;
					nzcount++;
				}
			}
		} else if (outRtype == REALSXP) {
			/* We only keep NAs. */
			double *out_nzvals_p = (double *) out_nzvals;
			for (int k = nzcount = 0; k < in_nzcount; k++) {
				int x = nzvals[k];
				if (x == NA_INTEGER) {
					out_nzoffs[nzcount] = sv->nzoffs[k];
					out_nzvals_p[nzcount] = NA_REAL;
					nzcount++;
				}
			}
		}
	} else if (Rtype == REALSXP) {
		if (outRtype == REALSXP) {
			nzcount = Arith_doubleSV_double(MULT_OPCODE,
					sv, 0.0,
					out_nzoffs, (double *) out_nzvals);
		}
	}
	if (nzcount == -1)
		error("_mult_SV_zero() only supports input "
		      "of type \"integer\" or \"double\" at the moment");
	return nzcount;
}

int _Arith_sv1_sv2(int opcode, const SparseVec *sv1, const SparseVec *sv2,
		SEXPTYPE expected_outRtype,
		int *out_nzoffs, void *out_nzvals, int *ovflow)
{
	SEXPTYPE effective_outRtype = REALSXP;
	int nzcount = -1;
	SEXPTYPE Rtype1 = get_SV_Rtype(sv1);
	SEXPTYPE Rtype2 = get_SV_Rtype(sv2);
	if (Rtype1 == INTSXP) {
		if (Rtype2 == INTSXP) {
			effective_outRtype = INTSXP;
			nzcount = Arith_intSV_intSV(opcode, sv1, sv2,
					out_nzoffs, (int *) out_nzvals,
					ovflow);
		} else if (Rtype2 == REALSXP) {
			nzcount = Arith_intSV_doubleSV(opcode, sv1, sv2,
					out_nzoffs, (double *) out_nzvals);
		}
	} else if (Rtype1 == REALSXP) {
		if (Rtype2 == INTSXP) {
			nzcount = Arith_doubleSV_intSV(opcode, sv1, sv2,
					out_nzoffs, (double *) out_nzvals);
		} else if (Rtype2 == REALSXP) {
			nzcount = Arith_doubleSV_doubleSV(opcode, sv1, sv2,
					out_nzoffs, (double *) out_nzvals);
		}
	}
	if (nzcount == -1)
		error("_Arith_sv1_sv2() only supports input of "
		      "type \"integer\" or \"double\" at the moment");
	if (expected_outRtype != effective_outRtype)
		error("SparseArray internal error in "
		      "_Arith_sv1_sv2():\n"
		      "    expected_outRtype != effective_outRtype");
	return nzcount;
}

