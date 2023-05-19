/****************************************************************************
 *                   Summarization of an R atomic vector                    *
 ****************************************************************************/
#include "Rvector_summarization.h"

#include <math.h>  /* for sqrt() */


/****************************************************************************
 * _get_summarize_opcode()
 */

int _get_summarize_opcode(SEXP op, SEXPTYPE Rtype)
{
	const char *s;

	if (!IS_CHARACTER(op) || LENGTH(op) != 1)
		error("'op' must be a single string");
	op = STRING_ELT(op, 0);
	if (op == NA_STRING)
		error("'op' cannot be NA");
	s = CHAR(op);
	if (Rtype != LGLSXP && Rtype != INTSXP && Rtype != REALSXP)
		error("%s() does not support SparseArray objects "
		      "of type() \"%s\"", s, type2char(Rtype));
	if (strcmp(s, "min") == 0)
		return MIN_OPCODE;
	if (strcmp(s, "max") == 0)
		return MAX_OPCODE;
	if (strcmp(s, "range") == 0)
		return RANGE_OPCODE;
	if (strcmp(s, "sum") == 0)
		return SUM_OPCODE;
	if (strcmp(s, "prod") == 0)
		return PROD_OPCODE;
	if (strcmp(s, "mean") == 0)
		return MEAN_OPCODE;
	if (strcmp(s, "sum_centered_X2") == 0)
		return SUM_CENTERED_X2_OPCODE;
	if (strcmp(s, "sum_X_X2") == 0)
		return SUM_X_X2_OPCODE;
	if (strcmp(s, "var1") == 0)
		return VAR1_OPCODE;
	if (strcmp(s, "var2") == 0)
		return VAR2_OPCODE;
	if (strcmp(s, "sd1") == 0)
		return SD1_OPCODE;
	if (strcmp(s, "sd2") == 0)
		return SD2_OPCODE;
	if (Rtype != LGLSXP && Rtype != INTSXP)
		error("%s() does not support SparseArray objects "
		      "of type() \"%s\"", s, type2char(Rtype));
	if (strcmp(s, "any") == 0)
		return ANY_OPCODE;
	if (strcmp(s, "all") == 0)
		return ALL_OPCODE;
	error("'op' must be one of: \"any\", \"all\", \"min\", \"max\", "
	      "\"range\", \"sum\", \"prod\",\n"
	      "                       \"mean\", "
	      "\"sum_centered_X2\", \"sum_X_X2\", \"var1\", \"var2\",\n"
	      "                       \"sd1\", \"sd2\"");
	return 0;
}


/****************************************************************************
 * _make_SummarizeOp()
 * _init_SummarizeResult()
 */

SummarizeOp _make_SummarizeOp(int opcode, SEXPTYPE in_Rtype,
			      int na_rm, double center)
{
	SummarizeOp summarize_op;

	summarize_op.opcode = opcode;
	summarize_op.in_Rtype = in_Rtype;
	summarize_op.na_rm = na_rm;
	summarize_op.center = center;
	return summarize_op;
}

void _init_SummarizeResult(const SummarizeOp *summarize_op,
			   SummarizeResult *res)
{
	res->totalcount = res->nzcount = res->nacount = 0;
	res->postprocess_one_zero = 0;
	res->outbuf_is_set = 1;
	res->warn = 0;
	switch (summarize_op->opcode) {
	    case ANY_OPCODE:
		res->out_Rtype = LGLSXP;
		res->outbuf.one_int[0] = 0;
		return;
	    case ALL_OPCODE:
		res->postprocess_one_zero = 1;
		res->out_Rtype = LGLSXP;
		res->outbuf.one_int[0] = 1;
		return;
	    case SUM_OPCODE: case MEAN_OPCODE:
		res->out_Rtype = REALSXP;
		res->outbuf.one_double[0] = 0.0;
		return;
	    case PROD_OPCODE:
		res->postprocess_one_zero = 1;
		res->out_Rtype = REALSXP;
		res->outbuf.one_double[0] = 1.0;
		return;
	    case SUM_CENTERED_X2_OPCODE: case VAR1_OPCODE: case SD1_OPCODE:
		res->out_Rtype = REALSXP;
		res->outbuf.one_double[0] = 0.0;
		return;
	    case SUM_X_X2_OPCODE: case VAR2_OPCODE: case SD2_OPCODE:
		res->out_Rtype = REALSXP;
		res->outbuf.two_doubles[0] = res->outbuf.two_doubles[1] = 0.0;
		return;
	}
	/* From now on, 'summarize_op->opcode' can only be MIN_OPCODE,
	   MAX_OPCODE, or RANGE_OPCODE. */
	res->postprocess_one_zero = 1;
	if (summarize_op->in_Rtype == INTSXP) {
		res->outbuf_is_set = 0;
		res->out_Rtype = INTSXP;
		return;
	}
	res->out_Rtype = REALSXP;
	switch (summarize_op->opcode) {
	    case MIN_OPCODE:
		res->outbuf.one_double[0] = R_PosInf;
		return;
	    case MAX_OPCODE:
		res->outbuf.one_double[0] = R_NegInf;
		return;
	    case RANGE_OPCODE:
		res->outbuf.two_doubles[0] = R_PosInf;
		res->outbuf.two_doubles[1] = R_NegInf;
		return;
	}
	return;
}


/****************************************************************************
 * Low-level summarization functions
 *
 * They all return 1 when they bail out early (on reaching a break condition),
 * and 0 otherwise.
 */

/* No need to read 'outbuf[0]' because we can assume that it's set to
   FALSE or NA (it cannot be TRUE). */
static inline int any_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		int outbuf[1])
{
	int set_outbuf_to_NA;

	set_outbuf_to_NA = 0;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			set_outbuf_to_NA = 1;
			continue;
		}
		if (*x != 0) {
			outbuf[0] = 1;
			return 1;  // bail out early
		}
	}
	if (set_outbuf_to_NA)
		outbuf[0] = NA_INTEGER;
	return 0;
}

/* No need to read 'outbuf[0]' because we can assume that it's set to
   TRUE or NA (it cannot be FALSE). */
static inline int all_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		int outbuf[1])
{
	int set_outbuf_to_NA;

	set_outbuf_to_NA = 0;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			set_outbuf_to_NA = 1;
			continue;
		}
		if (*x == 0) {
			outbuf[0] = 0;
			return 1;  // bail out early
		}
	}
	if (set_outbuf_to_NA)
		outbuf[0] = NA_INTEGER;
	return 0;
}

static inline int min_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		int *outbuf_is_set, int outbuf[1])
{
	int out0;

	out0 = outbuf[0];
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			outbuf[0] = NA_INTEGER;
			*outbuf_is_set = 1;
			return 1;  // bail out early
		}
		if (!*outbuf_is_set || *x < out0) {
			out0 = *x;
			*outbuf_is_set = 1;
		}
	}
	outbuf[0] = out0;
	return 0;
}

static inline int min_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[1])
{
	double out0, xx;
	int out0_is_not_NaN;

	out0 = outbuf[0];
	out0_is_not_NaN = !R_IsNaN(out0);
	for (int i = 0; i < n; i++, x++) {
		xx = *x;
		if (ISNAN(xx)) {  // True for *both* NA and NaN
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			if (R_IsNA(xx)) {
				outbuf[0] = NA_REAL;
				return 1;  // bail out early
			}
			out0 = xx;
			out0_is_not_NaN = 0;
			continue;
		}
		if (out0_is_not_NaN && xx < out0)
			out0 = xx;
	}
	outbuf[0] = out0;
	return 0;
}

static inline int max_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		int *outbuf_is_set, int outbuf[1])
{
	int out0;

	out0 = outbuf[0];
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			outbuf[0] = NA_INTEGER;
			*outbuf_is_set = 1;
			return 1;  // bail out early
		}
		if (!*outbuf_is_set || *x > out0) {
			out0 = *x;
			*outbuf_is_set = 1;
		}
	}
	outbuf[0] = out0;
	return 0;
}

static inline int max_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[1])
{
	double out0, xx;
	int out0_is_not_NaN;

	out0 = outbuf[0];
	out0_is_not_NaN = !R_IsNaN(out0);
	for (int i = 0; i < n; i++, x++) {
		xx = *x;
		if (ISNAN(xx)) {  // True for *both* NA and NaN
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			if (R_IsNA(xx)) {
				outbuf[0] = NA_REAL;
				return 1;  // bail out early
			}
			out0 = xx;
			out0_is_not_NaN = 0;
			continue;
		}
		if (out0_is_not_NaN && xx > out0)
			out0 = xx;
	}
	outbuf[0] = out0;
	return 0;
}

static inline int range_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		int *outbuf_is_set, int outbuf[2])
{
	int out0, out1;

	out0 = outbuf[0];
	out1 = outbuf[1];
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			outbuf[0] = outbuf[1] = NA_INTEGER;
			*outbuf_is_set = 1;
			return 1;  // bail out early
		}
		if (*outbuf_is_set) {
			if (*x < out0)
				out0 = *x;
			if (*x > out1)
				out1 = *x;
		} else {
			out0 = out1 = *x;
			*outbuf_is_set = 1;
		}
	}
	outbuf[0] = out0;
	outbuf[1] = out1;
	return 0;
}

static inline int range_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[2])
{
	double out0, out1, xx;
	int out0_is_not_NaN;

	out0 = outbuf[0];
	out1 = outbuf[1];
	out0_is_not_NaN = !R_IsNaN(out0);
	for (int i = 0; i < n; i++, x++) {
		xx = *x;
		if (ISNAN(xx)) {  // True for *both* NA and NaN
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			if (R_IsNA(xx)) {
				outbuf[0] = outbuf[1] = NA_REAL;
				return 1;  // bail out early
			}
			out0 = out1 = xx;
			out0_is_not_NaN = 0;
			continue;
		}
		if (out0_is_not_NaN) {
			if (xx < out0)
				out0 = xx;
			if (xx > out1)
				out1 = xx;
		}
	}
	outbuf[0] = out0;
	outbuf[1] = out1;
	return 0;
}

static inline int sum_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[1])
{
	double out0;

	out0 = outbuf[0];
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			outbuf[0] = NA_REAL;
			return 1;  // bail out early
		}
		out0 += (double) *x;
	}
	outbuf[0] = out0;
	return 0;
}

static inline int sum_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[1])
{
	double out0, xx;
	int out0_is_not_NaN;

	out0 = outbuf[0];
	out0_is_not_NaN = !R_IsNaN(out0);
	for (int i = 0; i < n; i++, x++) {
		xx = *x;
		if (ISNAN(xx)) {  // True for *both* NA and NaN
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			if (R_IsNA(xx)) {
				outbuf[0] = NA_REAL;
				return 1;  // bail out early
			}
			out0 = xx;
			out0_is_not_NaN = 0;
			continue;
		}
		if (out0_is_not_NaN)
			out0 += xx;
	}
	outbuf[0] = out0;
	return 0;
}

static inline int prod_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[1])
{
	double out0;

	out0 = outbuf[0];
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			outbuf[0] = NA_REAL;
			return 1;  // bail out early
		}
		out0 *= (double) *x;
	}
	outbuf[0] = out0;
	return 0;
}

static inline int prod_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[1])
{
	double out0, xx;
	int out0_is_not_NaN;

	out0 = outbuf[0];
	out0_is_not_NaN = !R_IsNaN(out0);
	for (int i = 0; i < n; i++, x++) {
		xx = *x;
		if (ISNAN(xx)) {  // True for *both* NA and NaN
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			if (R_IsNA(xx)) {
				outbuf[0] = NA_REAL;
				return 1;  // bail out early
			}
			out0 = xx;
			out0_is_not_NaN = 0;
			continue;
		}
		if (out0_is_not_NaN)
			out0 *= xx;
	}
	outbuf[0] = out0;
	return 0;
}

static inline int sum_centered_X2_ints(const int *x, int n,
		int na_rm, double center, R_xlen_t *nacount,
		double outbuf[1])
{
	double out0, y;

	out0 = outbuf[0];
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			outbuf[0] = NA_REAL;
			return 1;  // bail out early
		}
		y = (double) *x - center;
		out0 += y * y;
	}
	outbuf[0] = out0;
	return 0;
}

static inline int sum_centered_X2_doubles(const double *x, int n,
		int na_rm, double center, R_xlen_t *nacount,
		double outbuf[1])
{
	double out0, xx, y;
	int out0_is_not_NaN;

	out0 = outbuf[0];
	out0_is_not_NaN = !R_IsNaN(out0);
	for (int i = 0; i < n; i++, x++) {
		xx = *x;
		if (ISNAN(xx)) {  // True for *both* NA and NaN
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			if (R_IsNA(xx)) {
				outbuf[0] = NA_REAL;
				return 1;  // bail out early
			}
			out0 = xx;
			out0_is_not_NaN = 0;
			continue;
		}
		if (out0_is_not_NaN) {
			y = xx - center;
			out0 += y * y;
		}
	}
	outbuf[0] = out0;
	return 0;
}

static inline int sum_X_X2_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[2])
{
	double out0, out1, xx;

	out0 = outbuf[0];
	out1 = outbuf[1];
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			outbuf[0] = outbuf[1] = NA_REAL;
			return 1;  // bail out early
		}
		xx = (double) *x;
		out0 += xx;
		out1 += xx * xx;
	}
	outbuf[0] = out0;
	outbuf[1] = out1;
	return 0;
}

static inline int sum_X_X2_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[2])
{
	double out0, out1, xx;
	int out0_is_not_NaN;

	out0 = outbuf[0];
	out1 = outbuf[1];
	out0_is_not_NaN = !R_IsNaN(out0);
	for (int i = 0; i < n; i++, x++) {
		xx = *x;
		if (ISNAN(xx)) {  // True for *both* NA and NaN
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			if (R_IsNA(xx)) {
				outbuf[0] = outbuf[1] = NA_REAL;
				return 1;  // bail out early
			}
			out0 = out1 = xx;
			out0_is_not_NaN = 0;
			continue;
		}
		if (out0_is_not_NaN) {
			out0 += xx;
			out1 += xx * xx;
		}
	}
	outbuf[0] = out0;
	outbuf[1] = out1;
	return 0;
}


/****************************************************************************
 * _summarize_Rvector()
 */

static int summarize_ints(const int *x, int x_len,
		int opcode, int na_rm, double center, SummarizeResult *res)
{
	switch (opcode) {
	    case MIN_OPCODE:
		return min_ints(x, x_len, na_rm, &(res->nacount),
				&(res->outbuf_is_set), res->outbuf.one_int);
	    case MAX_OPCODE:
		return max_ints(x, x_len, na_rm, &(res->nacount),
				&(res->outbuf_is_set), res->outbuf.one_int);
	    case RANGE_OPCODE:
		return range_ints(x, x_len, na_rm, &(res->nacount),
				&(res->outbuf_is_set), res->outbuf.two_ints);
	    case SUM_OPCODE: case MEAN_OPCODE:
		return sum_ints(x, x_len, na_rm, &(res->nacount),
				res->outbuf.one_double);
	    case PROD_OPCODE:
		return prod_ints(x, x_len, na_rm, &(res->nacount),
				 res->outbuf.one_double);
	    case ANY_OPCODE:
		return any_ints(x, x_len, na_rm, &(res->nacount),
				res->outbuf.one_int);
	    case ALL_OPCODE:
		return all_ints(x, x_len, na_rm, &(res->nacount),
				res->outbuf.one_int);
	    case SUM_CENTERED_X2_OPCODE: case VAR1_OPCODE: case SD1_OPCODE:
		return sum_centered_X2_ints(x, x_len, na_rm, center,
				&(res->nacount), res->outbuf.one_double);
	    case SUM_X_X2_OPCODE: case VAR2_OPCODE: case SD2_OPCODE:
		return sum_X_X2_ints(x, x_len, na_rm, &(res->nacount),
				res->outbuf.two_doubles);
	}
	error("SparseArray internal error in summarize_ints():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static int summarize_doubles(const double *x, int x_len,
		int opcode, int na_rm, double center, SummarizeResult *res)
{
	switch (opcode) {
	    case MIN_OPCODE:
		return min_doubles(x, x_len, na_rm, &(res->nacount),
				res->outbuf.one_double);
	    case MAX_OPCODE:
		return max_doubles(x, x_len, na_rm, &(res->nacount),
				res->outbuf.one_double);
	    case RANGE_OPCODE:
		return range_doubles(x, x_len, na_rm, &(res->nacount),
				res->outbuf.two_doubles);
	    case SUM_OPCODE: case MEAN_OPCODE:
		return sum_doubles(x, x_len, na_rm, &(res->nacount),
				res->outbuf.one_double);
	    case PROD_OPCODE:
		return prod_doubles(x, x_len, na_rm, &(res->nacount),
				res->outbuf.one_double);
	    case SUM_CENTERED_X2_OPCODE: case VAR1_OPCODE: case SD1_OPCODE:
		return sum_centered_X2_doubles(x, x_len, na_rm, center,
				&(res->nacount), res->outbuf.one_double);
	    case SUM_X_X2_OPCODE: case VAR2_OPCODE: case SD2_OPCODE:
		return sum_X_X2_doubles(x, x_len, na_rm, &(res->nacount),
				res->outbuf.two_doubles);
	}
	error("SparseArray internal error in summarize_doubles():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

int _summarize_Rvector(SEXP x, const SummarizeOp *summarize_op,
		       SummarizeResult *res)
{
	SEXPTYPE x_Rtype;
	int x_len, bailout;

	x_Rtype = TYPEOF(x);
	if (x_Rtype != summarize_op->in_Rtype)
		error("SparseArray internal error in _summarize_Rvector():\n"
		      "    x_Rtype != summarize_op->in_Rtype");
	x_len = LENGTH(x);
	res->totalcount += x_len;
	switch (x_Rtype) {
	    case LGLSXP: case INTSXP:
		bailout = summarize_ints(INTEGER(x), x_len,
				summarize_op->opcode, summarize_op->na_rm,
				summarize_op->center, res);
		break;
	    case REALSXP:
		bailout = summarize_doubles(REAL(x), x_len,
				summarize_op->opcode, summarize_op->na_rm,
				summarize_op->center, res);
		break;
	    default:
		error("SparseArray internal error in _summarize_Rvector():\n"
		      "    input type \"%s\" is not supported",
		      type2char(x_Rtype));
	}
	if (bailout)
		res->postprocess_one_zero = 0;
	return bailout;
}


/****************************************************************************
 * _postprocess_SummarizeResult()
 */

/* Does NOT increase 'res->totalcount' by 1. */
static int summarize_one_zero(const SummarizeOp *summarize_op,
			      SummarizeResult *res)
{
	switch (summarize_op->in_Rtype) {
	    case LGLSXP: case INTSXP: {
		int zero = 0;
		return summarize_ints(&zero, 1,
				summarize_op->opcode, summarize_op->na_rm,
				summarize_op->center, res);
	    }
	    case REALSXP: {
		double zero = 0.0;
		return summarize_doubles(&zero, 1,
				summarize_op->opcode, summarize_op->na_rm,
				summarize_op->center, res);
	    }
	}
	error("SparseArray internal error in summarize_one_zero():\n"
	      "    input type \"%s\" is not supported",
	      type2char(summarize_op->in_Rtype));
	return 0;  /* will never reach this */
}


void _postprocess_SummarizeResult(const SummarizeOp *summarize_op,
				  SummarizeResult *res)
{
	int opcode;
	R_xlen_t zerocount, totalcount;

	opcode = summarize_op->opcode;
	zerocount = res->totalcount - res->nzcount;
	totalcount = res->totalcount;
	if (summarize_op->na_rm)
		totalcount -= res->nacount;

	if (res->postprocess_one_zero && zerocount != 0)
		summarize_one_zero(summarize_op, res);

	if (!res->outbuf_is_set) {
		if (res->out_Rtype == INTSXP && (opcode == MIN_OPCODE ||
						 opcode == MAX_OPCODE ||
						 opcode == RANGE_OPCODE))
		{
			/* Will happen if the virtual vector we're summarizing
			   has length 0 (i.e. 'res->totalcount == 0'), or if
			   it contains only NAs (i.e. 'res->nacount ==
			   res->totalcount') and 'summarize_op->na_rm' is True.
			   This is a case where we intentional deviate from
			   base::min(), base::max(), and base::range(). */
			if (opcode == RANGE_OPCODE) {
				res->outbuf.two_ints[0] =
				res->outbuf.two_ints[1] = NA_INTEGER;
			} else {
				res->outbuf.one_int[0] = NA_INTEGER;
			}
			res->warn = res->outbuf_is_set = 1;
		} else {
			error("SparseArray internal error in "
			      "_postprocess_SummarizeResult():\n"
			      "    'res->outbuf_is_set' is False");
		}
	}

	switch (opcode) {
	    case MEAN_OPCODE:
		res->outbuf.one_double[0] /= (double) totalcount;
		return;

	    case SUM_CENTERED_X2_OPCODE: case VAR1_OPCODE: case SD1_OPCODE:
		double center = summarize_op->center;
		res->outbuf.one_double[0] += center * center * zerocount;
		if (opcode == SUM_CENTERED_X2_OPCODE)
			return;
		if (totalcount <= 1) {
			res->outbuf.one_double[0] = NA_REAL;
			return;
		}
		res->outbuf.one_double[0] /= (totalcount - 1.0);
		if (opcode == VAR1_OPCODE)
			return;
		res->outbuf.one_double[0] = sqrt(res->outbuf.one_double[0]);
		return;

	    case VAR2_OPCODE: case SD2_OPCODE:
		if (totalcount <= 1) {
			res->outbuf.one_double[0] = NA_REAL;
			return;
		}
		double sum_X  = res->outbuf.two_doubles[0];
		double sum_X2 = res->outbuf.two_doubles[1];
		double var2   = (sum_X2 - sum_X * sum_X / totalcount) /
			        (totalcount - 1.0);
		res->outbuf.one_double[0] = var2;
		if (opcode == VAR2_OPCODE)
			return;
		res->outbuf.one_double[0] = sqrt(res->outbuf.one_double[0]);
		return;
	}
	return;
}


/****************************************************************************
 * _make_SEXP_from_summarize_result()
 */

/* ScalarLogical() is kind of broken e.g. it seems to have the 3 possible
   results precomputed. This means that if we set attributes on the result
   of, say, ScalarLogical(1) (like we do in _make_SEXP_from_summarize_result()
   below), then any subsequent call to ScalarLogical(1) will return a TRUE
   with the same attributes on it. In other words, we've permanently altered
   the behavior of ScalarLogical(1) for the current session!
   Note that this will alter base R functionalities as well e.g.:

     > any(matrix(c(TRUE, TRUE, FALSE, NA)))
     [1] TRUE
     attr(,"nacount")
     [1] 0

   ScalarLogical2() fixes that. Caller must PROTECT() the result. */
static SEXP ScalarLogical2(int i)
{
	SEXP ans;

	ans = duplicate(PROTECT(ScalarLogical(i)));
	UNPROTECT(1);
	return ans;
}

static SEXP sum_X_X2_as_SEXP(double sum_X, double sum_X2, SEXPTYPE in_Rtype)
{
	SEXP ans;

	/* Either 'sum_X' and 'sum_X2' are both set to NA or NaN or none
	   of them is. Furthermore, in the former case, they should both
	   be set to the same kind of NA i.e. both are set to NA or both
	   are set NaN. */
	if (in_Rtype == INTSXP) {
		/* Note that when 'in_Rtype' is INTSXP, the only kind of
		   NAs that can end up in 'sum_X' is NA_REAL. No NaNs. */
		if (ISNAN(sum_X)) {
			ans = PROTECT(NEW_INTEGER(2));
			INTEGER(ans)[0] = INTEGER(ans)[1] = NA_INTEGER;
			UNPROTECT(1);
			return ans;
		}
		if (sum_X  <= INT_MAX && sum_X  >= -INT_MAX &&
		    sum_X2 <= INT_MAX && sum_X2 >= -INT_MAX)
		{
			ans = PROTECT(NEW_INTEGER(2));
			/* We round 'sum_X' and 'sum_X2' to the
			   nearest integer. */
			INTEGER(ans)[0] = (int) (sum_X  + 0.5);
			INTEGER(ans)[1] = (int) (sum_X2 + 0.5);
			UNPROTECT(1);
			return ans;
		}
	}
	ans = PROTECT(NEW_NUMERIC(2));
	if (ISNAN(sum_X)) {  // True for *both* NA and NaN
		/* 'sum_X' is NA_REAL or NaN. */
		REAL(ans)[0] = REAL(ans)[1] = sum_X;
	} else {
		REAL(ans)[0] = sum_X;
		REAL(ans)[1] = sum_X2;
	}
	UNPROTECT(1);
	return ans;
}

/* Returns an integer or numeric vector of length 1 or 2.
   Looks at 'res->outbuf_is_set' only when 'opcode' is MIN/MAX/RANGE_OPCODE
   and 'Rtype' is INTSXP. */
static SEXP res2nakedSEXP(const SummarizeResult *res,
			  int opcode, SEXPTYPE in_Rtype)
{
	SEXP ans;

	if (opcode == ANY_OPCODE || opcode == ALL_OPCODE)
		return ScalarLogical2(res->outbuf.one_int[0]);

	if (opcode == MIN_OPCODE || opcode == MAX_OPCODE) {
		if (in_Rtype == REALSXP)
			return ScalarReal(res->outbuf.one_double[0]);
		if (in_Rtype == INTSXP || in_Rtype == LGLSXP)
			return ScalarInteger(res->outbuf.one_int[0]);
		error("SparseArray internal error in res2nakedSEXP()\n"
		      "    type \"%s\" not supported for min() or max()",
		      type2char(in_Rtype));
	}

	if (opcode == RANGE_OPCODE) {
		if (in_Rtype == REALSXP) {
			ans = PROTECT(NEW_NUMERIC(2));
			REAL(ans)[0] = res->outbuf.two_doubles[0];
			REAL(ans)[1] = res->outbuf.two_doubles[1];
		} else if (in_Rtype == INTSXP || in_Rtype == LGLSXP) {
			ans = PROTECT(NEW_INTEGER(2));
			INTEGER(ans)[0] = res->outbuf.two_ints[0];
			INTEGER(ans)[1] = res->outbuf.two_ints[1];
		} else {
			error("SparseArray internal error in res2nakedSEXP()\n"
			      "    type \"%s\" not supported for range()",
			      type2char(in_Rtype));
		}
		UNPROTECT(1);
		return ans;
	}

	if (opcode == SUM_X_X2_OPCODE)
		return sum_X_X2_as_SEXP(res->outbuf.two_doubles[0],
					res->outbuf.two_doubles[1], in_Rtype);

	if (in_Rtype == INTSXP && (opcode == SUM_OPCODE ||
				   opcode == PROD_OPCODE))
	{
		if (ISNAN(res->outbuf.one_double[0]))
			return ScalarInteger(NA_INTEGER);
		if (res->outbuf.one_double[0] < -INT_MAX ||
		    res->outbuf.one_double[0] >  INT_MAX)
			return ScalarReal(res->outbuf.one_double[0]);
		/* Round 'res->outbuf.one_double[0]' to the nearest integer. */
		return ScalarInteger((int) (res->outbuf.one_double[0] + 0.5));
	}

	return ScalarReal(res->outbuf.one_double[0]);
}

/* Returns an integer or numeric vector of length 1 or 2.
   If 'na_rm' is TRUE, then the "nacount" attribute is set on the
   returned vector. */
SEXP _make_SEXP_from_summarize_result(const SummarizeOp *summarize_op,
		const SummarizeResult *res)
{
	SEXP ans, ans_attrib;

	ans = res2nakedSEXP(res, summarize_op->opcode, summarize_op->in_Rtype);
	if (!summarize_op->na_rm)
		return ans;
	PROTECT(ans);
	if (res->nacount > INT_MAX)
		ans_attrib = ScalarReal((double) res->nacount);
	else
		ans_attrib = ScalarInteger((int) res->nacount);
	PROTECT(ans_attrib);
	setAttrib(ans, install("nacount"), ans_attrib);
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * _count_Rvector_NAs() and _Rvector_has_any_NA()
 */

static int count_NA_int_elts(const int *x, int n)
{
	int count;

	for (int i = count = 0; i < n; i++, x++)
		if (*x == NA_INTEGER)
			count++;
	return count;
}
static int any_NA_int_elt(const int *x, int n)
{
	for (int i = 0; i < n; i++, x++)
		if (*x == NA_INTEGER)
			return 1;
	return 0;
}

static int count_NA_double_elts(const double *x, int n)
{
	int count;

	for (int i = count = 0; i < n; i++, x++)
		if (ISNAN(*x))  // True for *both* NA and NaN
			count++;
	return count;
}
static int any_NA_double_elt(const double *x, int n)
{
	for (int i = 0; i < n; i++, x++)
		if (ISNAN(*x))  // True for *both* NA and NaN
			return 1;
	return 0;
}

#define RCOMPLEX_IS_NA(x) (ISNAN((x)->r) || ISNAN((x)->i))
static int count_NA_Rcomplex_elts(const Rcomplex *x, int n)
{
	int count;

	for (int i = count = 0; i < n; i++, x++)
		if (RCOMPLEX_IS_NA(x))
			count++;
	return count;
}
static int any_NA_Rcomplex_elt(const Rcomplex *x, int n)
{
	for (int i = 0; i < n; i++, x++)
		if (RCOMPLEX_IS_NA(x))
			return 1;
	return 0;
}

static int count_NA_character_elts(SEXP x)
{
	int n, count;

	n = LENGTH(x);
	for (int i = count = 0; i < n; i++)
		if (STRING_ELT(x, i) == NA_STRING)
			count++;
	return count;
}
static int any_NA_character_elt(SEXP x)
{
	int n;

	n = LENGTH(x);
	for (int i = 0; i < n; i++)
		if (STRING_ELT(x, i) == NA_STRING)
			return 1;
	return 0;
}

static inline int is_single_NA(SEXP x)
{
	switch (TYPEOF(x)) {
	    case LGLSXP: case INTSXP:
		if (LENGTH(x) == 1 && INTEGER(x)[0] == NA_INTEGER)
			return 1;
	    break;
	    case REALSXP:
		if (LENGTH(x) == 1 && ISNAN(REAL(x)[0]))
			return 1;
	    break;
	    case CPLXSXP:
		if (LENGTH(x) == 1 && RCOMPLEX_IS_NA(COMPLEX(x)))
			return 1;
	    break;
	    case STRSXP:
		if (LENGTH(x) == 1 && STRING_ELT(x, 0) == NA_STRING)
			return 1;
	    break;
	}
	return 0;
}

static int count_NA_list_elts(SEXP x)
{
	int n, count;

	n = LENGTH(x);
	for (int i = count = 0; i < n; i++)
		if (is_single_NA(VECTOR_ELT(x, i)))
			count++;
	return count;
}
static int any_NA_list_elt(SEXP x)
{
	int n;

	n = LENGTH(x);
	for (int i = 0; i < n; i++) {
		if (is_single_NA(VECTOR_ELT(x, i)))
			return 1;
	}
	return 0;
}

int _count_Rvector_NAs(SEXP Rvector)
{
	SEXPTYPE Rtype;
	int n;

	Rtype = TYPEOF(Rvector);
	n = LENGTH(Rvector);
	switch (Rtype) {
	    case LGLSXP: case INTSXP:
			  return count_NA_int_elts(INTEGER(Rvector), n);
	    case REALSXP: return count_NA_double_elts(REAL(Rvector), n);
	    case CPLXSXP: return count_NA_Rcomplex_elts(COMPLEX(Rvector), n);
	    case RAWSXP:  return 0;
	    case STRSXP:  return count_NA_character_elts(Rvector);
	    case VECSXP:  return count_NA_list_elts(Rvector);
	}
	error("SparseArray internal error in _count_Rvector_NAs():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return -1;  /* will never reach this */
}

int _Rvector_has_any_NA(SEXP Rvector)
{
	SEXPTYPE Rtype;
	int n;

	Rtype = TYPEOF(Rvector);
	n = LENGTH(Rvector);
	switch (Rtype) {
	    case LGLSXP: case INTSXP:
			  return any_NA_int_elt(INTEGER(Rvector), n);
	    case REALSXP: return any_NA_double_elt(REAL(Rvector), n);
	    case CPLXSXP: return any_NA_Rcomplex_elt(COMPLEX(Rvector), n);
	    case RAWSXP:  return 0;
	    case STRSXP:  return any_NA_character_elt(Rvector);
	    case VECSXP:  return any_NA_list_elt(Rvector);
	}
	error("SparseArray internal error in _Rvector_has_any_NA():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return -1;  /* will never reach this */
}

