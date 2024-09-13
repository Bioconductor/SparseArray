/****************************************************************************
 ****************************************************************************
 **									   **
 **                  Summarization of an atomic R vector                   **
 **									   **
 ****************************************************************************
 ****************************************************************************/
#include "Rvector_summarization.h"

#include "Rvector_utils.h"

#include <math.h>  /* for sqrt() */


/****************************************************************************
 * _get_summarize_opcode()
 */

int _get_summarize_opcode(SEXP op, SEXPTYPE Rtype)
{
	if (!(IS_CHARACTER(op) && LENGTH(op) == 1))
		error("'op' must be a single string");
	op = STRING_ELT(op, 0);
	if (op == NA_STRING)
		error("'op' cannot be NA");
	const char *s = CHAR(op);

	if (Rtype != LGLSXP && Rtype != INTSXP && Rtype != REALSXP &&
	    Rtype != CPLXSXP && Rtype != STRSXP)
		error("%s() does not support SparseArray objects "
		      "of type() \"%s\"", s, type2char(Rtype));

	if (strcmp(s, "anyNA") == 0)
		return ANYNA_OPCODE;
	if (strcmp(s, "countNAs") == 0)
		return COUNTNAS_OPCODE;
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
	if (strcmp(s, "centered_X2_sum") == 0)
		return CENTERED_X2_SUM_OPCODE;
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
	error("'op' must be one of: "
	      "\"anyNA\", \"countNAs\", \"any\", \"all\",\n"
	      "                       \"min\", \"max\", "
	      "\"range\", \"sum\", \"prod\", \"mean\",\n"
	      "                       \"centered_X2_sum\", \"sum_X_X2\",\n"
	      "                       \"var1\", \"var2\", \"sd1\", \"sd2\"");
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
	res->in_length = res->in_nzcount = res->in_nacount = 0;
	res->outbuf_status = OUTBUF_IS_SET;
	res->postprocess_one_zero = 0;
	res->warn = 0;
	switch (summarize_op->opcode) {
	    case ANYNA_OPCODE: case ANY_OPCODE:
		res->out_Rtype = LGLSXP;
		res->outbuf.one_int[0] = 0;
		return;
	    case COUNTNAS_OPCODE:
		res->out_Rtype = REALSXP;
		res->outbuf.one_double[0] = 0.0;
		return;
	    case ALL_OPCODE:
		res->out_Rtype = LGLSXP;
		res->outbuf.one_int[0] = 1;
		res->postprocess_one_zero = 1;
		return;
	    case SUM_OPCODE: case MEAN_OPCODE:
		res->out_Rtype = REALSXP;
		res->outbuf.one_double[0] = 0.0;
		return;
	    case PROD_OPCODE:
		res->out_Rtype = REALSXP;
		res->outbuf.one_double[0] = 1.0;
		res->postprocess_one_zero = 1;
		return;
	    case CENTERED_X2_SUM_OPCODE: case VAR1_OPCODE: case SD1_OPCODE:
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
	if (summarize_op->in_Rtype == LGLSXP ||
	    summarize_op->in_Rtype == INTSXP)
	{
		res->out_Rtype = INTSXP;
		res->outbuf_status = OUTBUF_IS_NOT_SET;
		return;
	}
	if (summarize_op->in_Rtype == REALSXP) {
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
	}
	error("SparseArray internal error in _init_SummarizeResult():\n"
	      "    operation not supported on SparseArray objects "
	      "of type() \"%s\"", type2char(summarize_op->in_Rtype));
	return;
}


/****************************************************************************
 * Low-level summarization utilities - part I
 *
 * Support R summarization functions that use Interface 1: FUN(x)
 *
 * For all of them, 'outbuf' is initialized by _init_SummarizeResult() above.
 * They all return the new "outbuf status".
 */

static inline int anyNA_ints(const int *x, int n, int outbuf[1])
{
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			/* Bail out early. */
			outbuf[0] = 1;
			return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
		}
	}
	return OUTBUF_IS_SET;
}

static inline int anyNA_doubles(const double *x, int n, int outbuf[1])
{
	for (int i = 0; i < n; i++, x++) {
		if (ISNAN(*x)) {  // True for *both* NA and NaN
			/* Bail out early. */
			outbuf[0] = 1;
			return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
		}
	}
	return OUTBUF_IS_SET;
}

static inline int anyNA_Rcomplexes(const Rcomplex *x, int n, int outbuf[1])
{
	for (int i = 0; i < n; i++, x++) {
		if (RCOMPLEX_IS_NA_OR_NaN(x)) {
			/* Bail out early. */
			outbuf[0] = 1;
			return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
		}
	}
	return OUTBUF_IS_SET;
}

static inline int anyNA_Rstrings(SEXP x, int outbuf[1])
{
	int n = LENGTH(x);
	for (int i = 0; i < n; i++) {
		if (STRING_ELT(x, i) == NA_STRING) {
			/* Bail out early. */
			outbuf[0] = 1;
			return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
		}
	}
	return OUTBUF_IS_SET;
}

static inline int countNAs_ints(const int *x, int n, double outbuf[1])
{
	double out0 = outbuf[0];
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER)
			out0++;
	}
	outbuf[0] = out0;
	return OUTBUF_IS_SET;
}

static inline int countNAs_doubles(const double *x, int n, double outbuf[1])
{
	double out0 = outbuf[0];
	for (int i = 0; i < n; i++, x++) {
		if (ISNAN(*x))
			out0++;
	}
	outbuf[0] = out0;
	return OUTBUF_IS_SET;
}

static inline int countNAs_Rcomplexes(const Rcomplex *x, int n,
				      double outbuf[1])
{
	double out0 = outbuf[0];
	for (int i = 0; i < n; i++, x++) {
		if (RCOMPLEX_IS_NA_OR_NaN(x))
			out0++;
	}
	outbuf[0] = out0;
	return OUTBUF_IS_SET;
}

static inline int countNAs_Rstrings(SEXP x, double outbuf[1])
{
	double out0 = outbuf[0];
	int n = LENGTH(x);
	for (int i = 0; i < n; i++) {
		if (STRING_ELT(x, i) == NA_STRING)
			out0++;
	}
	outbuf[0] = out0;
	return OUTBUF_IS_SET;
}


/****************************************************************************
 * Low-level summarization utilities - part II
 *
 * Support R summarization functions that use Interface 2 & Interface 3:
 * - Interface 2: FUN(x, na.rm)
 * - Interface 3: FUN(x, na.rm, center)
 *
 * They all return the new "outbuf status".
 */

/* 'outbuf' initialized by _init_SummarizeResult() above.
   No need to read 'outbuf[0]' because we can assume that it's set to
   FALSE or NA (it cannot be TRUE). */
static inline int any_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		int outbuf[1])
{
	int set_outbuf_to_NA = 0;
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
			/* Bail out early. */
			outbuf[0] = 1;
			return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
		}
	}
	if (set_outbuf_to_NA)
		outbuf[0] = NA_INTEGER;
	return OUTBUF_IS_SET;
}

/* 'outbuf' initialized by _init_SummarizeResult() above.
   No need to read 'outbuf[0]' because we can assume that it's set to
   TRUE or NA (it cannot be FALSE). */
static inline int all_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		int outbuf[1])
{
	int set_outbuf_to_NA = 0;
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
			/* Bail out early. */
			outbuf[0] = 0;
			return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
		}
	}
	if (set_outbuf_to_NA)
		outbuf[0] = NA_INTEGER;
	return OUTBUF_IS_SET;
}

/* 'outbuf' NOT initialized by _init_SummarizeResult() above. */
static inline int min_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		int outbuf[1], int outbuf_status)
{
	int out0 = outbuf[0];
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			/* Bail out early. */
			outbuf[0] = NA_INTEGER;
			return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
		}
		if (outbuf_status == OUTBUF_IS_NOT_SET || *x < out0) {
			out0 = *x;
			outbuf_status = OUTBUF_IS_SET;
		}
	}
	outbuf[0] = out0;
	return outbuf_status;
}

/* 'outbuf' initialized by _init_SummarizeResult() above. */
static inline int min_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[1])
{
	double out0 = outbuf[0];
	int out0_is_not_NaN = !R_IsNaN(out0);
	for (int i = 0; i < n; i++, x++) {
		double xx = *x;
		if (ISNAN(xx)) {  // True for *both* NA and NaN
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			if (R_IsNA(xx)) {
				/* Bail out early. */
				outbuf[0] = NA_REAL;
				return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
			}
			out0 = xx;
			out0_is_not_NaN = 0;
			continue;
		}
		if (out0_is_not_NaN && xx < out0)
			out0 = xx;
	}
	outbuf[0] = out0;
	return OUTBUF_IS_SET;
}

/* 'outbuf' NOT initialized by _init_SummarizeResult() above. */
static inline int max_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		int outbuf[1], int outbuf_status)
{
	int out0 = outbuf[0];
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			/* Bail out early. */
			outbuf[0] = NA_INTEGER;
			return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
		}
		if (outbuf_status == OUTBUF_IS_NOT_SET || *x > out0) {
			out0 = *x;
			outbuf_status = OUTBUF_IS_SET;
		}
	}
	outbuf[0] = out0;
	return outbuf_status;
}

/* 'outbuf' initialized by _init_SummarizeResult() above. */
static inline int max_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[1])
{
	double out0 = outbuf[0];
	int out0_is_not_NaN = !R_IsNaN(out0);
	for (int i = 0; i < n; i++, x++) {
		double xx = *x;
		if (ISNAN(xx)) {  // True for *both* NA and NaN
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			if (R_IsNA(xx)) {
				/* Bail out early. */
				outbuf[0] = NA_REAL;
				return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
			}
			out0 = xx;
			out0_is_not_NaN = 0;
			continue;
		}
		if (out0_is_not_NaN && xx > out0)
			out0 = xx;
	}
	outbuf[0] = out0;
	return OUTBUF_IS_SET;
}

/* 'outbuf' NOT initialized by _init_SummarizeResult() above. */
static inline int range_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		int outbuf[2], int outbuf_status)
{
	int out0 = outbuf[0];
	int out1 = outbuf[1];
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			/* Bail out early. */
			outbuf[0] = outbuf[1] = NA_INTEGER;
			return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
		}
		if (outbuf_status == OUTBUF_IS_NOT_SET) {
			out0 = out1 = *x;
			outbuf_status = OUTBUF_IS_SET;
		} else {
			if (*x < out0)
				out0 = *x;
			if (*x > out1)
				out1 = *x;
		}
	}
	outbuf[0] = out0;
	outbuf[1] = out1;
	return outbuf_status;
}

/* 'outbuf' initialized by _init_SummarizeResult() above. */
static inline int range_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[2])
{
	double out0 = outbuf[0];
	double out1 = outbuf[1];
	int out0_is_not_NaN = !R_IsNaN(out0);
	for (int i = 0; i < n; i++, x++) {
		double xx = *x;
		if (ISNAN(xx)) {  // True for *both* NA and NaN
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			if (R_IsNA(xx)) {
				/* Bail out early. */
				outbuf[0] = outbuf[1] = NA_REAL;
				return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
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
	return OUTBUF_IS_SET;
}

/* 'outbuf' initialized by _init_SummarizeResult() above. */
static inline int sum_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[1])
{
	double out0 = outbuf[0];
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			/* Bail out early. */
			outbuf[0] = NA_REAL;
			return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
		}
		out0 += (double) *x;
	}
	outbuf[0] = out0;
	return OUTBUF_IS_SET;
}

/* 'outbuf' initialized by _init_SummarizeResult() above. */
static inline int sum_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[1])
{
	double out0 = outbuf[0];
	int out0_is_not_NaN = !R_IsNaN(out0);
	for (int i = 0; i < n; i++, x++) {
		double xx = *x;
		if (ISNAN(xx)) {  // True for *both* NA and NaN
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			if (R_IsNA(xx)) {
				/* Bail out early. */
				outbuf[0] = NA_REAL;
				return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
			}
			out0 = xx;
			out0_is_not_NaN = 0;
			continue;
		}
		if (out0_is_not_NaN)
			out0 += xx;
	}
	outbuf[0] = out0;
	return OUTBUF_IS_SET;
}

/* 'outbuf' initialized by _init_SummarizeResult() above. */
static inline int prod_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[1])
{
	double out0 = outbuf[0];
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			/* Bail out early. */
			outbuf[0] = NA_REAL;
			return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
		}
		out0 *= (double) *x;
	}
	outbuf[0] = out0;
	return OUTBUF_IS_SET;
}

/* 'outbuf' initialized by _init_SummarizeResult() above. */
static inline int prod_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[1])
{
	double out0 = outbuf[0];
	int out0_is_not_NaN = !R_IsNaN(out0);
	for (int i = 0; i < n; i++, x++) {
		double xx = *x;
		if (ISNAN(xx)) {  // True for *both* NA and NaN
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			if (R_IsNA(xx)) {
				/* Bail out early. */
				outbuf[0] = NA_REAL;
				return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
			}
			out0 = xx;
			out0_is_not_NaN = 0;
			continue;
		}
		if (out0_is_not_NaN)
			out0 *= xx;
	}
	outbuf[0] = out0;
	return OUTBUF_IS_SET;
}

/* 'outbuf' initialized by _init_SummarizeResult() above. */
static inline int centered_X2_sum_ints(const int *x, int n,
		int na_rm, double center, R_xlen_t *nacount,
		double outbuf[1])
{
	double out0 = outbuf[0];
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			/* Bail out early. */
			outbuf[0] = NA_REAL;
			return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
		}
		double delta = (double) *x - center;
		out0 += delta * delta;
	}
	outbuf[0] = out0;
	return OUTBUF_IS_SET;
}

/* 'outbuf' initialized by _init_SummarizeResult() above. */
static inline int centered_X2_sum_doubles(const double *x, int n,
		int na_rm, double center, R_xlen_t *nacount,
		double outbuf[1])
{
	double out0 = outbuf[0];
	int out0_is_not_NaN = !R_IsNaN(out0);
	for (int i = 0; i < n; i++, x++) {
		double xx = *x;
		if (ISNAN(xx)) {  // True for *both* NA and NaN
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			if (R_IsNA(xx)) {
				/* Bail out early. */
				outbuf[0] = NA_REAL;
				return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
			}
			out0 = xx;
			out0_is_not_NaN = 0;
			continue;
		}
		if (out0_is_not_NaN) {
			double delta = xx - center;
			out0 += delta * delta;
		}
	}
	outbuf[0] = out0;
	return OUTBUF_IS_SET;
}

/* 'outbuf' initialized by _init_SummarizeResult() above. */
static inline int sum_X_X2_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[2])
{
	double out0 = outbuf[0];
	double out1 = outbuf[1];
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			/* Bail out early. */
			outbuf[0] = outbuf[1] = NA_REAL;
			return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
		}
		double xx = (double) *x;
		out0 += xx;
		out1 += xx * xx;
	}
	outbuf[0] = out0;
	outbuf[1] = out1;
	return OUTBUF_IS_SET;
}

/* 'outbuf' initialized by _init_SummarizeResult() above. */
static inline int sum_X_X2_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[2])
{
	double out0 = outbuf[0];
	double out1 = outbuf[1];
	int out0_is_not_NaN = !R_IsNaN(out0);
	for (int i = 0; i < n; i++, x++) {
		double xx = *x;
		if (ISNAN(xx)) {  // True for *both* NA and NaN
			if (na_rm) {
				(*nacount)++;
				continue;
			}
			if (R_IsNA(xx)) {
				/* Bail out early. */
				outbuf[0] = outbuf[1] = NA_REAL;
				return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
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
	return OUTBUF_IS_SET;
}


/****************************************************************************
 * _summarize_ones()
 * _summarize_Rvector()
 */

static int summarize_ones(SEXPTYPE x_Rtype, int x_len,
		int opcode, double center, SummarizeResult *res)
{
	if (x_len == 0)
		return res->outbuf_status;
	switch (opcode) {
	    case ANYNA_OPCODE:
	    case COUNTNAS_OPCODE:
	    case ALL_OPCODE:
		return OUTBUF_IS_SET;
	    case ANY_OPCODE:
		res->outbuf.one_int[0] = 1;
		return OUTBUF_IS_SET_WITH_BREAKING_VALUE;
	    case MIN_OPCODE:
		if (x_Rtype == INTSXP || x_Rtype == LGLSXP) {
			if (res->outbuf_status == OUTBUF_IS_NOT_SET ||
			    res->outbuf.one_int[0] > int1)
			{
				res->outbuf.one_int[0] = int1;
			}
		} else {
			if (res->outbuf_status == OUTBUF_IS_NOT_SET ||
			    res->outbuf.one_double[0] > double1)
			{
				res->outbuf.one_double[0] = double1;
			}
		}
		return OUTBUF_IS_SET;
	    case MAX_OPCODE:
		if (x_Rtype == INTSXP || x_Rtype == LGLSXP) {
			if (res->outbuf_status == OUTBUF_IS_NOT_SET ||
			    res->outbuf.one_int[0] < int1)
			{
				res->outbuf.one_int[0] = int1;
			}
		} else {
			if (res->outbuf_status == OUTBUF_IS_NOT_SET ||
			    res->outbuf.one_double[0] < double1)
			{
				res->outbuf.one_double[0] = double1;
			}
		}
		return OUTBUF_IS_SET;
	    case RANGE_OPCODE:
		if (x_Rtype == INTSXP || x_Rtype == LGLSXP) {
			if (res->outbuf_status == OUTBUF_IS_NOT_SET) {
				res->outbuf.two_ints[0] =
				res->outbuf.two_ints[1] = int1;
			} else {
				if (res->outbuf.two_ints[0] > int1)
					res->outbuf.two_ints[0] = int1;
				if (res->outbuf.two_ints[1] < int1)
					res->outbuf.two_ints[1] = int1;
			}
		} else {
			if (res->outbuf_status == OUTBUF_IS_NOT_SET) {
				res->outbuf.two_doubles[0] =
				res->outbuf.two_doubles[1] = double1;
			} else {
				if (res->outbuf.two_doubles[0] > double1)
					res->outbuf.two_doubles[0] = double1;
				if (res->outbuf.two_doubles[1] < double1)
					res->outbuf.two_doubles[1] = double1;
			}
		}
		return OUTBUF_IS_SET;
	    case SUM_OPCODE: case MEAN_OPCODE:
		res->outbuf.one_double[0] += (double) x_len;
	    case PROD_OPCODE:
		return OUTBUF_IS_SET;
	    case CENTERED_X2_SUM_OPCODE: case VAR1_OPCODE: case SD1_OPCODE: {
		double delta = 1.0 - center;
		res->outbuf.one_double[0] += delta * delta * x_len;
		return OUTBUF_IS_SET;
	    }
	    case SUM_X_X2_OPCODE: case VAR2_OPCODE: case SD2_OPCODE:
		res->outbuf.two_doubles[0] += 1.0;
		res->outbuf.two_doubles[1] += 1.0;
		return OUTBUF_IS_SET;
	}
	error("SparseArray internal error in summarize_ones():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static int summarize_ints(const int *x, int x_len,
		int opcode, int na_rm, double center, SummarizeResult *res)
{
	R_xlen_t *nacount_p = &(res->in_nacount);
	switch (opcode) {
	    case ANYNA_OPCODE:
		return anyNA_ints(x, x_len, res->outbuf.one_int);
	    case COUNTNAS_OPCODE:
		return countNAs_ints(x, x_len, res->outbuf.one_double);
	    case ANY_OPCODE:
		return any_ints(x, x_len, na_rm, nacount_p,
				res->outbuf.one_int);
	    case ALL_OPCODE:
		return all_ints(x, x_len, na_rm, nacount_p,
				res->outbuf.one_int);
	    case MIN_OPCODE:
		return min_ints(x, x_len, na_rm, nacount_p,
				res->outbuf.one_int, res->outbuf_status);
	    case MAX_OPCODE:
		return max_ints(x, x_len, na_rm, nacount_p,
				res->outbuf.one_int, res->outbuf_status);
	    case RANGE_OPCODE:
		return range_ints(x, x_len, na_rm, nacount_p,
				res->outbuf.two_ints, res->outbuf_status);
	    case SUM_OPCODE: case MEAN_OPCODE:
		return sum_ints(x, x_len, na_rm, nacount_p,
				res->outbuf.one_double);
	    case PROD_OPCODE:
		return prod_ints(x, x_len, na_rm, nacount_p,
				res->outbuf.one_double);
	    case CENTERED_X2_SUM_OPCODE: case VAR1_OPCODE: case SD1_OPCODE:
		return centered_X2_sum_ints(x, x_len, na_rm,
				center, nacount_p,
				res->outbuf.one_double);
	    case SUM_X_X2_OPCODE: case VAR2_OPCODE: case SD2_OPCODE:
		return sum_X_X2_ints(x, x_len, na_rm, nacount_p,
				res->outbuf.two_doubles);
	}
	error("SparseArray internal error in summarize_ints():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static int summarize_doubles(const double *x, int x_len,
		int opcode, int na_rm, double center, SummarizeResult *res)
{
	R_xlen_t *nacount_p = &(res->in_nacount);
	switch (opcode) {
	    case ANYNA_OPCODE:
		return anyNA_doubles(x, x_len, res->outbuf.one_int);
	    case COUNTNAS_OPCODE:
		return countNAs_doubles(x, x_len, res->outbuf.one_double);
	    case MIN_OPCODE:
		return min_doubles(x, x_len, na_rm, nacount_p,
				res->outbuf.one_double);
	    case MAX_OPCODE:
		return max_doubles(x, x_len, na_rm, nacount_p,
				res->outbuf.one_double);
	    case RANGE_OPCODE:
		return range_doubles(x, x_len, na_rm, nacount_p,
				res->outbuf.two_doubles);
	    case SUM_OPCODE: case MEAN_OPCODE:
		return sum_doubles(x, x_len, na_rm, nacount_p,
				res->outbuf.one_double);
	    case PROD_OPCODE:
		return prod_doubles(x, x_len, na_rm, nacount_p,
				res->outbuf.one_double);
	    case CENTERED_X2_SUM_OPCODE: case VAR1_OPCODE: case SD1_OPCODE:
		return centered_X2_sum_doubles(x, x_len, na_rm,
				center, nacount_p,
				res->outbuf.one_double);
	    case SUM_X_X2_OPCODE: case VAR2_OPCODE: case SD2_OPCODE:
		return sum_X_X2_doubles(x, x_len, na_rm, nacount_p,
				res->outbuf.two_doubles);
	}
	error("SparseArray internal error in summarize_doubles():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

/* Arguments 'na_rm' and 'center' ignored at the moment. */
static int summarize_Rcomplexes(const Rcomplex *x, int x_len,
		int opcode, int na_rm, double center, SummarizeResult *res)
{
	switch (opcode) {
	    case ANYNA_OPCODE:
		return anyNA_Rcomplexes(x, x_len, res->outbuf.one_int);
	    case COUNTNAS_OPCODE:
		return countNAs_Rcomplexes(x, x_len, res->outbuf.one_double);
	}
	error("SparseArray internal error in summarize_Rcomplexes():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

/* Arguments 'na_rm' and 'center' ignored at the moment. */
static int summarize_Rstrings(SEXP x,
		int opcode, int na_rm, double center, SummarizeResult *res)
{
	switch (opcode) {
	    case ANYNA_OPCODE:
		return anyNA_Rstrings(x, res->outbuf.one_int);
	    case COUNTNAS_OPCODE:
		return countNAs_Rstrings(x, res->outbuf.one_double);
	}
	error("SparseArray internal error in summarize_Rstrings():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

/* Summarizes a fictive vector of ones. */
void _summarize_ones(int x_len, const SummarizeOp *summarize_op,
		     SummarizeResult *res)
{
	if (res->outbuf_status == OUTBUF_IS_SET_WITH_BREAKING_VALUE)
		error("SparseArray internal error in _summarize_ones():\n"
		      "    outbuf already set with breaking value");
	res->in_length += x_len;
	int new_status = summarize_ones(summarize_op->in_Rtype, x_len,
					summarize_op->opcode,
					summarize_op->center, res);
	res->outbuf_status = new_status;
	if (new_status == OUTBUF_IS_SET_WITH_BREAKING_VALUE)
		res->postprocess_one_zero = 0;
	return;
}

void _summarize_Rvector(SEXP x, const SummarizeOp *summarize_op,
			SummarizeResult *res)
{
	if (res->outbuf_status == OUTBUF_IS_SET_WITH_BREAKING_VALUE)
		error("SparseArray internal error in _summarize_Rvector():\n"
		      "    outbuf already set with breaking value");
	SEXPTYPE x_Rtype = TYPEOF(x);
	if (x_Rtype != summarize_op->in_Rtype)
		error("SparseArray internal error in _summarize_Rvector():\n"
		      "    x_Rtype != summarize_op->in_Rtype");
	int x_len = LENGTH(x);
	res->in_length += x_len;
	int new_status;
	switch (x_Rtype) {
	    case INTSXP: case LGLSXP:
		new_status = summarize_ints(INTEGER(x), x_len,
				summarize_op->opcode, summarize_op->na_rm,
				summarize_op->center, res);
		break;
	    case REALSXP:
		new_status = summarize_doubles(REAL(x), x_len,
				summarize_op->opcode, summarize_op->na_rm,
				summarize_op->center, res);
		break;
	    case CPLXSXP:
		new_status = summarize_Rcomplexes(COMPLEX(x), x_len,
				summarize_op->opcode, summarize_op->na_rm,
				summarize_op->center, res);
		break;
	    case STRSXP:
		new_status = summarize_Rstrings(x,
				summarize_op->opcode, summarize_op->na_rm,
				summarize_op->center, res);
		break;
	    default:
		error("SparseArray internal error in _summarize_Rvector():\n"
		      "    input type \"%s\" is not supported",
		      type2char(x_Rtype));
	}
	res->outbuf_status = new_status;
	if (new_status == OUTBUF_IS_SET_WITH_BREAKING_VALUE)
		res->postprocess_one_zero = 0;
	return;
}


/****************************************************************************
 * _postprocess_SummarizeResult()
 */

/* Does NOT increase 'res->in_length' by 1. */
static void summarize_one_zero(const SummarizeOp *summarize_op,
			       SummarizeResult *res)
{
	if (res->outbuf_status == OUTBUF_IS_SET_WITH_BREAKING_VALUE)
		error("SparseArray internal error in summarize_one_zero():\n"
		      "    outbuf already set with breaking value");
	int new_status;
	switch (summarize_op->in_Rtype) {
	    case INTSXP: case LGLSXP:
		new_status = summarize_ints(&int0, 1,
				summarize_op->opcode, summarize_op->na_rm,
				summarize_op->center, res);
		break;
	    case REALSXP:
		new_status = summarize_doubles(&double0, 1,
				summarize_op->opcode, summarize_op->na_rm,
				summarize_op->center, res);
		break;
	    default:
		error("SparseArray internal error in summarize_one_zero():\n"
		      "    input type \"%s\" is not supported",
		      type2char(summarize_op->in_Rtype));
	}
	res->outbuf_status = new_status;
	return;
}

/* Does NOT increase 'res->in_length' by 1. */
static void summarize_one_NA(const SummarizeOp *summarize_op,
			     SummarizeResult *res)
{
	if (res->outbuf_status == OUTBUF_IS_SET_WITH_BREAKING_VALUE)
		error("SparseArray internal error in summarize_one_NA():\n"
		      "    outbuf already set with breaking value");
	if (summarize_op->na_rm)
		return;
	int new_status;
	switch (summarize_op->in_Rtype) {
	    case INTSXP: case LGLSXP: {
		new_status = summarize_ints(&intNA, 1,
				summarize_op->opcode, 0,
				summarize_op->center, res);
		break;
	    }
	    case REALSXP: {
		new_status = summarize_doubles(&doubleNA, 1,
				summarize_op->opcode, 0,
				summarize_op->center, res);
		break;
	    }
	    case CPLXSXP: {
		new_status = summarize_Rcomplexes(&RcomplexNA, 1,
				summarize_op->opcode, 0,
				summarize_op->center, res);
		break;
	    }
	    case STRSXP: {
		SEXP val = PROTECT(ScalarString(NA_STRING));
		new_status = summarize_Rstrings(val,
				summarize_op->opcode, 0,
				summarize_op->center, res);
		UNPROTECT(1);
		break;
	    }
	    default:
		error("SparseArray internal error in summarize_one_NA():\n"
		      "    input type \"%s\" is not supported",
		      type2char(summarize_op->in_Rtype));
	}
	res->outbuf_status = new_status;
	return;
}

void _postprocess_SummarizeResult(SummarizeResult *res, int na_background,
				  const SummarizeOp *summarize_op)
{
	/* There's nothing to do if a break condition was reached. */
	if (res->outbuf_status == OUTBUF_IS_SET_WITH_BREAKING_VALUE)
		return;

	int opcode = summarize_op->opcode;
	R_xlen_t zerocount = res->in_length - res->in_nzcount;
	if (opcode == COUNTNAS_OPCODE) {
		if (na_background)
			res->outbuf.one_double[0] += zerocount;
		return;
	}
	R_xlen_t effective_len = res->in_length;

	if (summarize_op->na_rm) {
		if (na_background)
			effective_len = res->in_nzcount;
		effective_len -= res->in_nacount;
	}

	if (zerocount != 0) {
		if (na_background) {
			summarize_one_NA(summarize_op, res);
		} else if (res->postprocess_one_zero) {
			summarize_one_zero(summarize_op, res);
		}
	}

	if (res->outbuf_status == OUTBUF_IS_NOT_SET) {
		if ((opcode == MIN_OPCODE || opcode == MAX_OPCODE ||
		     opcode == RANGE_OPCODE) &&
		    (res->out_Rtype == LGLSXP || res->out_Rtype == INTSXP))
		{
			/* Will happen if the virtual vector we're summarizing
			   has length 0 (i.e. 'res->in_length == 0'), or if
			   it contains only NAs (i.e. 'res->in_nacount ==
			   res->in_length') and 'summarize_op->na_rm' is True.
			   This is a case where we intentionally deviate from
			   base::min(), base::max(), and base::range(). */
			if (opcode == RANGE_OPCODE) {
				res->outbuf.two_ints[0] =
				res->outbuf.two_ints[1] = NA_INTEGER;
			} else {
				res->outbuf.one_int[0] = NA_INTEGER;
			}
			res->warn = 1;
			res->outbuf_status = OUTBUF_IS_SET;
			return;
		}
		error("SparseArray internal error in "
		      "_postprocess_SummarizeResult():\n"
		      "    outbuf is not set");
	}

	switch (opcode) {
	    /* For some strange reasons the Apple clang compiler is not
	       happy if we don't use curly brackets to wrap the blocks
	       following the case labels.
	       See https://github.com/Bioconductor/SparseArray/issues/3 */
	    case MEAN_OPCODE: {
		res->outbuf.one_double[0] /= (double) effective_len;
		return;
	    }
	    case CENTERED_X2_SUM_OPCODE: case VAR1_OPCODE: case SD1_OPCODE: {
		double center = summarize_op->center;
		if (!na_background)
			res->outbuf.one_double[0] +=
				center * center * zerocount;
		if (opcode == CENTERED_X2_SUM_OPCODE)
			return;
		if (effective_len <= 1) {
			res->outbuf.one_double[0] = NA_REAL;
			return;
		}
		res->outbuf.one_double[0] /= (effective_len - 1.0);
		if (opcode == VAR1_OPCODE)
			return;
		res->outbuf.one_double[0] = sqrt(res->outbuf.one_double[0]);
		return;
	    }
	    case VAR2_OPCODE: case SD2_OPCODE: {
		if (effective_len <= 1) {
			res->outbuf.one_double[0] = NA_REAL;
			return;
		}
		double sum_X  = res->outbuf.two_doubles[0];
		double sum_X2 = res->outbuf.two_doubles[1];
		double var2   = (sum_X2 - sum_X * sum_X / effective_len) /
			        (effective_len - 1.0);
		res->outbuf.one_double[0] = var2;
		if (opcode == VAR2_OPCODE)
			return;
		res->outbuf.one_double[0] = sqrt(res->outbuf.one_double[0]);
		return;
	    }
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
	SEXP ans = duplicate(PROTECT(ScalarLogical(i)));
	UNPROTECT(1);
	return ans;
}

/* Round 'x' to nearest int. */
#define	BACK_TO_INT(x) ((int) ((x) >= 0 ? (x) + 0.5 : (x) - 0.5))

static SEXP sum_X_X2_as_SEXP(double sum_X, double sum_X2, SEXPTYPE in_Rtype)
{
	/* Either 'sum_X' and 'sum_X2' are both set to NA or NaN or none
	   of them is. Furthermore, in the former case, they should both
	   be set to the same kind of NA i.e. both are set to NA or both
	   are set NaN. */
	if (in_Rtype == INTSXP) {
		/* Note that when 'in_Rtype' is INTSXP, the only kind of
		   NAs that can end up in 'sum_X' is NA_REAL. No NaNs. */
		if (ISNAN(sum_X)) {
			SEXP ans = PROTECT(NEW_INTEGER(2));
			INTEGER(ans)[0] = INTEGER(ans)[1] = NA_INTEGER;
			UNPROTECT(1);
			return ans;
		}
		if (sum_X  <= INT_MAX && sum_X  >= -INT_MAX &&
		    sum_X2 <= INT_MAX && sum_X2 >= -INT_MAX)
		{
			SEXP ans = PROTECT(NEW_INTEGER(2));
			INTEGER(ans)[0] = BACK_TO_INT(sum_X);
			INTEGER(ans)[1] = BACK_TO_INT(sum_X2);
			UNPROTECT(1);
			return ans;
		}
	}
	SEXP ans = PROTECT(NEW_NUMERIC(2));
	if (ISNAN(sum_X)) {
		/* 'sum_X' is NA_REAL or NaN. */
		REAL(ans)[0] = REAL(ans)[1] = sum_X;
	} else {
		REAL(ans)[0] = sum_X;
		REAL(ans)[1] = sum_X2;
	}
	UNPROTECT(1);
	return ans;
}

/* Returns an integer or numeric vector of length 1 or 2. */
static SEXP res2nakedSEXP(const SummarizeResult *res,
			  int opcode, SEXPTYPE in_Rtype)
{
	if (opcode == ANYNA_OPCODE ||
	    opcode == ANY_OPCODE || opcode == ALL_OPCODE)
	{
		return ScalarLogical2(res->outbuf.one_int[0]);
	}

	if (opcode == COUNTNAS_OPCODE) {
		double out0 = res->outbuf.one_double[0];
		if (out0 > INT_MAX)
			return ScalarReal(out0);
		/* Round 'out0' to the nearest integer. */
		return ScalarInteger(BACK_TO_INT(out0));
	}

	if ((opcode == MIN_OPCODE || opcode == MAX_OPCODE) &&
	    in_Rtype != REALSXP)
	{
		return ScalarInteger(res->outbuf.one_int[0]);
	}

	if (opcode == RANGE_OPCODE) {
		SEXP ans;
		if (in_Rtype == REALSXP) {
			ans = PROTECT(NEW_NUMERIC(2));
			REAL(ans)[0] = res->outbuf.two_doubles[0];
			REAL(ans)[1] = res->outbuf.two_doubles[1];
		} else {
			ans = PROTECT(NEW_INTEGER(2));
			INTEGER(ans)[0] = res->outbuf.two_ints[0];
			INTEGER(ans)[1] = res->outbuf.two_ints[1];
		}
		UNPROTECT(1);
		return ans;
	}

	if ((opcode == SUM_OPCODE || opcode == PROD_OPCODE) &&
	    (in_Rtype == LGLSXP || in_Rtype == INTSXP))
	{
		double out0 = res->outbuf.one_double[0];
		if (ISNAN(out0))
			return ScalarInteger(NA_INTEGER);
		if (out0 < -INT_MAX || out0 > INT_MAX)
			return ScalarReal(out0);
		/* Round 'out0' to the nearest integer. */
		return ScalarInteger(BACK_TO_INT(out0));
	}

	if (opcode == SUM_X_X2_OPCODE)
		return sum_X_X2_as_SEXP(res->outbuf.two_doubles[0],
					res->outbuf.two_doubles[1], in_Rtype);

	return ScalarReal(res->outbuf.one_double[0]);
}

/* Returns an integer or numeric vector of length 1 or 2. */
SEXP _make_SEXP_from_summarize_result(const SummarizeOp *summarize_op,
		const SummarizeResult *res)
{
	return res2nakedSEXP(res, summarize_op->opcode, summarize_op->in_Rtype);
/*
	// If 'na_rm' is TRUE, then we set the "nacount" attribute on the
	// returned vector.
	SEXP ans, ans_attrib;

	ans = res2nakedSEXP(res, summarize_op->opcode, summarize_op->in_Rtype);
	if (!summarize_op->na_rm)
		return ans;
	PROTECT(ans);
	if (res->in_nacount > INT_MAX)
		ans_attrib = ScalarReal((double) res->in_nacount);
	else
		ans_attrib = ScalarInteger((int) res->in_nacount);
	PROTECT(ans_attrib);
	setAttrib(ans, install("nacount"), ans_attrib);
	UNPROTECT(2);
	return ans;
*/
}


/****************************************************************************
 * any_NA_list_elt() / count_NA_list_elts()
 *
 * NOT USED but we'll need something like this when we implement the
 * anyNA_list() and countNAs_list() utilities in order to make ops
 * ANYNA_OPCODE and COUNTNAS_OPCODE work on type "list".
 */

static inline int is_single_NA(SEXP x)
{
	switch (TYPEOF(x)) {
	    case INTSXP: case LGLSXP:
		if (LENGTH(x) == 1 && INTEGER(x)[0] == NA_INTEGER)
			return 1;
	    break;
	    case REALSXP:
		if (LENGTH(x) == 1 && ISNAN(REAL(x)[0]))
			return 1;
	    break;
	    case CPLXSXP:
		if (LENGTH(x) == 1 && RCOMPLEX_IS_NA_OR_NaN(COMPLEX(x)))
			return 1;
	    break;
	    case STRSXP:
		if (LENGTH(x) == 1 && STRING_ELT(x, 0) == NA_STRING)
			return 1;
	    break;
	}
	return 0;
}

static int any_NA_list_elt(SEXP x)
{
	int n = LENGTH(x);
	for (int i = 0; i < n; i++) {
		if (is_single_NA(VECTOR_ELT(x, i)))
			return 1;
	}
	return 0;
}

static int count_NA_list_elts(SEXP x)
{
	int n = LENGTH(x);
	int count = 0;
	for (int i = 0; i < n; i++)
		if (is_single_NA(VECTOR_ELT(x, i)))
			count++;
	return count;
}

