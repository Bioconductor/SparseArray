/****************************************************************************
 *                   Summarization of an R atomic vector                    *
 ****************************************************************************/
#include "Rvector_summarization.h"


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
	if (strcmp(s, "sum_shifted_X2") == 0)
		return SUM_SHIFTED_X2_OPCODE;
	if (strcmp(s, "sum_X_X2") == 0)
		return SUM_X_X2_OPCODE;
	if (Rtype == REALSXP)
		error("%s() does not support SparseArray objects "
		      "of type() \"%s\"", s, type2char(Rtype));
	if (strcmp(s, "any") == 0)
		return ANY_OPCODE;
	if (strcmp(s, "all") == 0)
		return ALL_OPCODE;
	error("'op' must be one of: \"min\", \"max\", \"range\", "
	      "\"sum\", \"prod\", \"any\", \"all\",\n"
	      "                       \"sum_shifted_X2\", \"sum_X_X2\"");
	return 0;
}


/****************************************************************************
 * Low-level summarization functions
 *
 * They all return 1 when they bail out early (on reaching a break condition),
 * and 0 otherwise.
 */

static inline int min_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		int *outbuf_is_set, int outbuf[1])
{
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				outbuf[0] = *x;
				*outbuf_is_set = 1;
				return 1;  // bail out early
			}
			(*nacount)++;
			continue;
		}
		if (!*outbuf_is_set || *x < outbuf[0]) {
			outbuf[0] = *x;
			*outbuf_is_set = 1;
		}
	}
	return 0;
}

static inline int min_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[1])
{
	for (int i = 0; i < n; i++, x++) {
		if (ISNAN(*x)) {  // True for *both* NA and NaN
			if (!na_rm) {
				outbuf[0] = *x;
				if (R_IsNA(*x))
					return 1;  // bail out early
			}
			(*nacount)++;
			continue;
		}
		if (!R_IsNaN(outbuf[0]) && *x < outbuf[0])
			outbuf[0] = *x;
	}
	return 0;
}

static inline int max_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		int *outbuf_is_set, int outbuf[1])
{
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				outbuf[0] = *x;
				*outbuf_is_set = 1;
				return 1;  // bail out early
			}
			(*nacount)++;
			continue;
		}
		if (!*outbuf_is_set || *x > outbuf[0]) {
			outbuf[0] = *x;
			*outbuf_is_set = 1;
		}
	}
	return 0;
}

static inline int max_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[1])
{
	for (int i = 0; i < n; i++, x++) {
		if (ISNAN(*x)) {  // True for *both* NA and NaN
			if (!na_rm) {
				outbuf[0] = *x;
				if (R_IsNA(*x))
					return 1;  // bail out early
			}
			(*nacount)++;
			continue;
		}
		if (!R_IsNaN(outbuf[0]) && *x > outbuf[0])
			outbuf[0] = *x;
	}
	return 0;
}

static inline int range_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		int *outbuf_is_set, int outbuf[2])
{
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				outbuf[0] = outbuf[1] = *x;
				*outbuf_is_set = 1;
				return 1;  // bail out early
			}
			(*nacount)++;
			continue;
		}
		if (!*outbuf_is_set) {
			outbuf[0] = outbuf[1] = *x;
			*outbuf_is_set = 1;
		} else {
			if (*x < outbuf[0])
				outbuf[0] = *x;
			if (*x > outbuf[1])
				outbuf[1] = *x;
		}
	}
	return 0;
}

static inline int range_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[2])
{
	for (int i = 0; i < n; i++, x++) {
		if (ISNAN(*x)) {  // True for *both* NA and NaN
			if (!na_rm) {
				outbuf[0] = outbuf[1] = *x;
				if (R_IsNA(*x))
					return 1;  // bail out early
			}
			(*nacount)++;
			continue;
		}
		if (!R_IsNaN(outbuf[0])) {
			if (*x < outbuf[0])
				outbuf[0] = *x;
			if (*x > outbuf[1])
				outbuf[1] = *x;
		}
	}
	return 0;
}

static inline int sum_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[1])
{
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				outbuf[0] = NA_REAL;
				return 1;  // bail out early
			}
			(*nacount)++;
			continue;
		}
		outbuf[0] += (double) *x;
	}
	return 0;
}

static inline int sum_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[1])
{
	for (int i = 0; i < n; i++, x++) {
		if (ISNAN(*x)) {  // True for *both* NA and NaN
			if (!na_rm) {
				outbuf[0] = *x;
				if (R_IsNA(*x))
					return 1;  // bail out early
			}
			(*nacount)++;
			continue;
		}
		if (!R_IsNaN(outbuf[0]))
			outbuf[0] += *x;
	}
	return 0;
}

static inline int prod_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[1])
{
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				outbuf[0] = NA_REAL;
				return 1;  // bail out early
			}
			(*nacount)++;
			continue;
		}
		outbuf[0] *= (double) *x;
	}
	return 0;
}

static inline int prod_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[1])
{
	for (int i = 0; i < n; i++, x++) {
		if (ISNAN(*x)) {  // True for *both* NA and NaN
			if (!na_rm) {
				outbuf[0] = *x;
				if (R_IsNA(*x))
					return 1;  // bail out early
			}
			(*nacount)++;
			continue;
		}
		if (!R_IsNaN(outbuf[0]))
			outbuf[0] *= *x;
	}
	return 0;
}

static inline int any_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		int outbuf[1])
{
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm)
				outbuf[0] = NA_INTEGER;
			(*nacount)++;
			continue;
		}
		if (*x != 0) {
			outbuf[0] = 1;
			return 1;  // bail out early
		}
	}
	return 0;
}

static inline int all_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		int outbuf[1])
{
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm)
				outbuf[0] = NA_INTEGER;
			(*nacount)++;
			continue;
		}
		if (*x == 0) {
			outbuf[0] = 0;
			return 1;  // bail out early
		}
	}
	return 0;
}

static inline int sum_shifted_X2_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[2])
{
	double y;

	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				outbuf[0] = NA_REAL;
				return 1;  // bail out early
			}
			(*nacount)++;
			continue;
		}
		y = (double) *x - outbuf[1];
		outbuf[0] += y * y;
	}
	return 0;
}

static inline int sum_shifted_X2_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[2])
{
	double y;

	for (int i = 0; i < n; i++, x++) {
		if (ISNAN(*x)) {  // True for *both* NA and NaN
			if (!na_rm) {
				outbuf[0] = *x;
				if (R_IsNA(*x))
					return 1;  // bail out early
			}
			(*nacount)++;
			continue;
		}
		if (!R_IsNaN(outbuf[0])) {
			y = *x - outbuf[1];
			outbuf[0] += y * y;
		}
	}
	return 0;
}

static inline int sum_X_X2_ints(const int *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[2])
{
	double xx;

	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				outbuf[0] = outbuf[1] = NA_REAL;
				return 1;  // bail out early
			}
			(*nacount)++;
			continue;
		}
		xx = (double) *x;
		outbuf[0] += xx;
		outbuf[1] += xx * xx;
	}
	return 0;
}

static inline int sum_X_X2_doubles(const double *x, int n,
		int na_rm, R_xlen_t *nacount,
		double outbuf[2])
{
	double xx;

	for (int i = 0; i < n; i++, x++) {
		xx = *x;
		if (ISNAN(xx)) {  // True for *both* NA and NaN
			if (!na_rm) {
				outbuf[0] = outbuf[1] = xx;
				if (R_IsNA(xx))
					return 1;  // bail out early
			}
			(*nacount)++;
			continue;
		}
		if (!R_IsNaN(outbuf[0])) {
			outbuf[0] += xx;
			outbuf[1] += xx * xx;
		}
	}
	return 0;
}

/* Returns an integer or numeric vector of length 1 or 2.
   Looks at 'res->outbuf_is_set' only when 'opcode' is MIN/MAX/RANGE_OPCODE
   and 'Rtype' is INTSXP. */
static SEXP res2nakedSEXP(const SummarizeResult *res,
			  int opcode, SEXPTYPE in_Rtype)
{
	SEXP ans;

	if (opcode == ANY_OPCODE || opcode == ALL_OPCODE)
		return ScalarLogical(res->outbuf.one_int[0]);

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

	if (opcode == SUM_X_X2_OPCODE) {
		/* Either 'res->outbuf.two_doubles[0]'
		   and 'res->outbuf.two_doubles[1]' are both set
		   to NA or NaN or none of them is. Furthermore, in the former
		   case, they should both be set to the same kind of NA i.e.
		   both are set to NA or both are set NaN. */
		double init0 = res->outbuf.two_doubles[0],
		       init1 = res->outbuf.two_doubles[1];
		if (in_Rtype == INTSXP) {
			/* Note that when 'in_Rtype' is INTSXP, the only
			   kind of NA that can end up in 'init0' is NA_REAL.
			   No NaN. */
			if (ISNAN(init0)) {
				ans = PROTECT(NEW_INTEGER(2));
				INTEGER(ans)[0] = INTEGER(ans)[1] = NA_INTEGER;
				UNPROTECT(1);
				return ans;
			}
			if (init0 <= INT_MAX && init0 >= -INT_MAX &&
			    init1 <= INT_MAX && init1 >= -INT_MAX)
			{
				ans = PROTECT(NEW_INTEGER(2));
				/* We round 'init0' and 'init1' to the
				   nearest integer. */
				INTEGER(ans)[0] = (int) (init0 + 0.5);
				INTEGER(ans)[1] = (int) (init1 + 0.5);
				UNPROTECT(1);
				return ans;
			}
		}
		ans = PROTECT(NEW_NUMERIC(2));
		if (ISNAN(init0)) {  // True for *both* NA and NaN
			/* init0' is NA_REAL or NaN. */
			REAL(ans)[0] = REAL(ans)[1] = init0;
		} else {
			REAL(ans)[0] = init0;
			REAL(ans)[1] = init1;
		}
		UNPROTECT(1);
		return ans;
	}

	/* 'opcode' is either SUM_OPCODE, PROD_OPCODE, or
	   SUM_SHIFTED_X2_OPCODE. */
	if (in_Rtype == REALSXP || opcode == SUM_SHIFTED_X2_OPCODE)
		return ScalarReal(res->outbuf.one_double[0]);

	if (ISNAN(res->outbuf.one_double[0]))  // True for *both* NA and NaN
		return ScalarInteger(NA_INTEGER);
	if (res->outbuf.one_double[0] <= INT_MAX &&
	    res->outbuf.one_double[0] >= -INT_MAX)
		/* Round 'res->outbuf.one_double[0]' to the nearest integer. */
		return ScalarInteger((int) (res->outbuf.one_double[0] + 0.5));
	return ScalarReal(res->outbuf.one_double[0]);
}


/****************************************************************************
 * _make_SummarizeOp()
 * _init_SummarizeResult()
 */

SummarizeOp _make_SummarizeOp(int opcode, SEXPTYPE in_Rtype,
			      int na_rm, double shift)
{
	SummarizeOp summarize_op;

	summarize_op.opcode = opcode;
	summarize_op.in_Rtype = in_Rtype;
	summarize_op.na_rm = na_rm;
	summarize_op.shift = shift;
	return summarize_op;
}

void _init_SummarizeResult(const SummarizeOp *summarize_op,
			   SummarizeResult *res)
{
	res->totalcount = res->nzcount = res->nacount = 0;
	res->outbuf_is_set = 1;
	res->warn = 0;
	switch (summarize_op->opcode) {
	    case ANY_OPCODE:
		res->out_Rtype = LGLSXP;
		res->outbuf.one_int[0] = 0;
		return;
	    case ALL_OPCODE:
		res->out_Rtype = LGLSXP;
		res->outbuf.one_int[0] = 1;
		return;
	    case SUM_OPCODE:
		res->out_Rtype = REALSXP;
		res->outbuf.one_double[0] = 0.0;
		return;
	    case PROD_OPCODE:
		res->out_Rtype = REALSXP;
		res->outbuf.one_double[0] = 1.0;
		return;
	    case SUM_SHIFTED_X2_OPCODE:
		res->out_Rtype = REALSXP;
		res->outbuf.two_doubles[0] = 0.0;
		res->outbuf.two_doubles[1] = summarize_op->shift;
		return;
	    case SUM_X_X2_OPCODE:
		res->out_Rtype = REALSXP;
		res->outbuf.two_doubles[0] = res->outbuf.two_doubles[1] = 0.0;
		return;
	}
	/* From now on, 'summarize_op->opcode' can only be MIN_OPCODE,
	   MAX_OPCODE, or RANGE_OPCODE. */
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
 * _summarize_Rvector()
 */

static int summarize_ints(const int *x, int x_len,
		int opcode, int na_rm, SummarizeResult *res)
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
	    case SUM_OPCODE:
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
	    case SUM_SHIFTED_X2_OPCODE:
		return sum_shifted_X2_ints(x, x_len, na_rm, &(res->nacount),
				res->outbuf.two_doubles);
	    case SUM_X_X2_OPCODE:
		return sum_X_X2_ints(x, x_len, na_rm, &(res->nacount),
				res->outbuf.two_doubles);
	}
	error("SparseArray internal error in summarize_ints():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

static int summarize_doubles(const double *x, int x_len,
		int opcode, int na_rm, SummarizeResult *res)
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
	    case SUM_OPCODE:
		return sum_doubles(x, x_len, na_rm, &(res->nacount),
				res->outbuf.one_double);
	    case PROD_OPCODE:
		return prod_doubles(x, x_len, na_rm, &(res->nacount),
				res->outbuf.one_double);
	    case SUM_SHIFTED_X2_OPCODE:
		return sum_shifted_X2_doubles(x, x_len, na_rm, &(res->nacount),
				res->outbuf.two_doubles);
	    case SUM_X_X2_OPCODE:
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
	int x_len, bailout;

	if (TYPEOF(x) != summarize_op->in_Rtype)
		error("SparseArray internal error in "
		      "_summarize_Rvector():\n"
		      "    TYPEOF(x) != summarize_op->in_Rtype");
	x_len = LENGTH(x);
	if (summarize_op->in_Rtype == INTSXP) {
		bailout = summarize_ints(INTEGER(x), x_len,
				summarize_op->opcode, summarize_op->na_rm,
				res);
	} else {
		bailout = summarize_doubles(REAL(x), x_len,
				summarize_op->opcode, summarize_op->na_rm,
				res);
	}
	res->totalcount += x_len;
	return bailout;
}


/****************************************************************************
 * _postprocess_SummarizeResult()
 */

/* Does NOT increase 'res->totalcount' by 1. */
static int summarize_one_zero(const SummarizeOp *summarize_op,
			      SummarizeResult *res)
{
	int bailout;

	if (summarize_op->in_Rtype == INTSXP) {
		int zero = 0;
		bailout = summarize_ints(&zero, 1,
				summarize_op->opcode, summarize_op->na_rm, res);
	} else {
		double zero = 0.0;
		bailout = summarize_doubles(&zero, 1,
				summarize_op->opcode, summarize_op->na_rm, res);
	}
	return bailout;
}


void _postprocess_SummarizeResult(const SummarizeOp *summarize_op,
				  SummarizeResult *res)
{
	int opcode;

	opcode = summarize_op->opcode;
	if (res->nzcount < res->totalcount && opcode != SUM_SHIFTED_X2_OPCODE)
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
			      "_postprocess_SummarizeResult():\n",
			      "    'res->outbuf_is_set' is False");
		}
	}
	return;
}


/****************************************************************************
 * _make_SEXP_from_summarize_result()
 */

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
	error("SparseArray internal error in "
	      "_count_Rvector_NAs():\n"
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
	error("SparseArray internal error in "
	      "_Rvector_has_any_NA():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return -1;  /* will never reach this */
}

