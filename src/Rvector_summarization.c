/****************************************************************************
 *                     Summarization an R atomic vector                     *
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
		      "of type \"%s\"", s, type2char(Rtype));
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
		      "of type \"%s\"", s, type2char(Rtype));
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
 * Callback functions used for summarize operations
 *
 * All these callback functions return the "new status":
 *    0 = 'init' has not been set yet
 *    1 = 'init' has been set
 *    2 = 'init' has been set and we don't need to continue (break condition)
 *
 * IMPORTANT NOTE: Most of them ignore the supplied 'status', only
 * min/max/range_ints() don't. This is because 'init' should have been set
 * by select_summarize_FUN() before the callback function gets called (see
 * below in this file). So it doesn't matter whether the supplied 'status'
 * is 0 or 1, and it doesn't matter if the callback function returns a new
 * status of 0 or 1.
 * Note that the init2nakedSEXP() function below in this file will ignore
 * the final status anyways, except when 'opcode' is MIN/MAX/RANGE_OPCODE
 * and 'Rtype' is INTSXP.
 * So for the callback functions that ignore the supplied 'status', the only
 * thing that matters is whether the returned 'status' is 2 or not, so the
 * caller knows whether to bail out or not (break condition).
 */

#define DOUBLE_IS_NA(x) (R_IsNA(x) || R_IsNaN(x))

/* Does NOT ignore 'status'. */
static inline int min_ints(void *init, const int *x, int n,
		int na_rm, R_xlen_t *na_rm_count, int status)
{
	int *int_init;

	int_init = (int *) init;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				int_init[0] = *x;
				return 2;
			}
			(*na_rm_count)++;
			continue;
		}
		if (status == 0 || *x < int_init[0]) {
			int_init[0] = *x;
			status = 1;
		}
	}
	return status;
}

/* Ignores 'status'. */
static inline int min_doubles(void *init, const double *x, int n,
		int na_rm, R_xlen_t *na_rm_count, int status)
{
	double *double_init;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (DOUBLE_IS_NA(*x)) {
			if (!na_rm) {
				double_init[0] = *x;
				if (R_IsNA(*x))
					return 2;
			}
			(*na_rm_count)++;
			continue;
		}
		if (!R_IsNaN(double_init[0]) && *x < double_init[0])
			double_init[0] = *x;
	}
	return status;
}

/* Does NOT ignore 'status'. */
static inline int max_ints(void *init, const int *x, int n,
		int na_rm, R_xlen_t *na_rm_count, int status)
{
	int *int_init;

	int_init = (int *) init;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				int_init[0] = *x;
				return 2;
			}
			(*na_rm_count)++;
			continue;
		}
		if (status == 0 || *x > int_init[0]) {
			int_init[0] = *x;
			status = 1;
		}
	}
	return status;
}

/* Ignores 'status'. */
static inline int max_doubles(void *init, const double *x, int n,
		int na_rm, R_xlen_t *na_rm_count, int status)
{
	double *double_init;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (DOUBLE_IS_NA(*x)) {
			if (!na_rm) {
				double_init[0] = *x;
				if (R_IsNA(*x))
					return 2;
			}
			(*na_rm_count)++;
			continue;
		}
		if (!R_IsNaN(double_init[0]) && *x > double_init[0])
			double_init[0] = *x;
	}
	return status;
}

/* Does NOT ignore 'status'. */
static inline int range_ints(void *init, const int *x, int n,
		int na_rm, R_xlen_t *na_rm_count, int status)
{
	int *int_init;

	int_init = (int *) init;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				int_init[0] = int_init[1] = *x;
				return 2;
			}
			(*na_rm_count)++;
			continue;
		}
		if (status == 0) {
			int_init[0] = int_init[1] = *x;
			status = 1;
		} else {
			if (*x < int_init[0])
				int_init[0] = *x;
			if (*x > int_init[1])
				int_init[1] = *x;
		}
	}
	return status;
}

/* Ignores 'status'. */
static inline int range_doubles(void *init, const double *x, int n,
		int na_rm, R_xlen_t *na_rm_count, int status)
{
	double *double_init;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (DOUBLE_IS_NA(*x)) {
			if (!na_rm) {
				double_init[0] = double_init[1] = *x;
				if (R_IsNA(*x))
					return 2;
			}
			(*na_rm_count)++;
			continue;
		}
		if (!R_IsNaN(double_init[0])) {
			if (*x < double_init[0])
				double_init[0] = *x;
			if (*x > double_init[1])
				double_init[1] = *x;
		}
	}
	return status;
}

/* Ignores 'status'. */
static inline int sum_ints(void *init, const int *x, int n,
		int na_rm, R_xlen_t *na_rm_count, int status)
{
	double *double_init;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				double_init[0] = NA_REAL;
				return 2;
			}
			(*na_rm_count)++;
			continue;
		}
		double_init[0] += (double) *x;
	}
	return status;
}

/* Ignores 'status'. */
static inline int sum_doubles(void *init, const double *x, int n,
		int na_rm, R_xlen_t *na_rm_count, int status)
{
	double *double_init;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (DOUBLE_IS_NA(*x)) {
			if (!na_rm) {
				double_init[0] = *x;
				if (R_IsNA(*x))
					return 2;
			}
			(*na_rm_count)++;
			continue;
		}
		if (!R_IsNaN(double_init[0]))
			double_init[0] += *x;
	}
	return status;
}

/* Ignores 'status'. */
static inline int prod_ints(void *init, const int *x, int n,
		int na_rm, R_xlen_t *na_rm_count, int status)
{
	double *double_init;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				double_init[0] = NA_REAL;
				return 2;
			}
			(*na_rm_count)++;
			continue;
		}
		double_init[0] *= (double) *x;
	}
	return status;
}

/* Ignores 'status'. */
static inline int prod_doubles(void *init, const double *x, int n,
		int na_rm, R_xlen_t *na_rm_count, int status)
{
	double *double_init;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (DOUBLE_IS_NA(*x)) {
			if (!na_rm) {
				double_init[0] = *x;
				if (R_IsNA(*x))
					return 2;
			}
			(*na_rm_count)++;
			continue;
		}
		if (!R_IsNaN(double_init[0]))
			double_init[0] *= *x;
	}
	return status;
}

/* Ignores 'status'. */
static inline int any_ints(void *init, const int *x, int n,
		int na_rm, R_xlen_t *na_rm_count, int status)
{
	int *int_init;

	int_init = (int *) init;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm)
				int_init[0] = *x;
			(*na_rm_count)++;
			continue;
		}
		if (*x != 0) {
			int_init[0] = *x;
			return 2;
		}
	}
	return status;
}

/* Ignores 'status'. */
static inline int all_ints(void *init, const int *x, int n,
		int na_rm, R_xlen_t *na_rm_count, int status)
{
	int *int_init;

	int_init = (int *) init;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm)
				int_init[0] = *x;
			(*na_rm_count)++;
			continue;
		}
		if (*x == 0) {
			int_init[0] = *x;
			return 2;
		}
	}
	return status;
}

/* Ignores 'status'. */
static inline int sum_shifted_X2_ints(void *init, const int *x, int n,
		int na_rm, R_xlen_t *na_rm_count, int status)
{
	double *double_init, y;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				double_init[0] = NA_REAL;
				return 2;
			}
			(*na_rm_count)++;
			continue;
		}
		y = (double) *x - double_init[1];
		double_init[0] += y * y;
	}
	return status;
}

/* Ignores 'status'. */
static inline int sum_shifted_X2_doubles(void *init, const double *x, int n,
		int na_rm, R_xlen_t *na_rm_count, int status)
{
	double *double_init, y;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (DOUBLE_IS_NA(*x)) {
			if (!na_rm) {
				double_init[0] = *x;
				if (R_IsNA(*x))
					return 2;
			}
			(*na_rm_count)++;
			continue;
		}
		if (!R_IsNaN(double_init[0])) {
			y = *x - double_init[1];
			double_init[0] += y * y;
		}
	}
	return status;
}

/* Ignores 'status'. */
static inline int sum_X_X2_ints(void *init, const int *x, int n,
		int na_rm, R_xlen_t *na_rm_count, int status)
{
	double *double_init, xx;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				double_init[0] = double_init[1] = NA_REAL;
				return 2;
			}
			(*na_rm_count)++;
			continue;
		}
		xx = (double) *x;
		double_init[0] += xx;
		double_init[1] += xx * xx;
	}
	return status;
}

/* Ignores 'status'. */
static inline int sum_X_X2_doubles(void *init, const double *x, int n,
		int na_rm, R_xlen_t *na_rm_count, int status)
{
	double *double_init, xx;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		xx = *x;
		if (DOUBLE_IS_NA(xx)) {
			if (!na_rm) {
				double_init[0] = double_init[1] = xx;
				if (R_IsNA(xx))
					return 2;
			}
			(*na_rm_count)++;
			continue;
		}
		if (!R_IsNaN(double_init[0])) {
			double_init[0] += xx;
			double_init[1] += xx * xx;
		}
	}
	return status;
}

/* Only one of '*summarize_ints_FUN' or '*summarize_doubles_FUN' will be set
   to a non-NULL value. The other one will be set to NULL. */
static void select_summarize_FUN(int opcode, SEXPTYPE Rtype, double shift,
		SummarizeInts_FUNType *summarize_ints_FUN,
		SummarizeDoubles_FUNType *summarize_doubles_FUN,
		void *init)
{
	int *int_init = (int *) init;
	double *double_init = (double *) init;

	*summarize_ints_FUN = NULL;
	*summarize_doubles_FUN = NULL;
	if (opcode == ANY_OPCODE) {
		*summarize_ints_FUN = any_ints;
		int_init[0] = 0;
		return;
	}
	if (opcode == ALL_OPCODE) {
		*summarize_ints_FUN = all_ints;
		int_init[0] = 1;
		return;
	}
	if (opcode == SUM_OPCODE) {
		if (Rtype == REALSXP) {
			*summarize_doubles_FUN = sum_doubles;
		} else {
			*summarize_ints_FUN = sum_ints;
		}
		double_init[0] = 0.0;
		return;
	}
	if (opcode == PROD_OPCODE) {
		if (Rtype == REALSXP) {
			*summarize_doubles_FUN = prod_doubles;
		} else {
			*summarize_ints_FUN = prod_ints;
		}
		double_init[0] = 1.0;
		return;
	}
	if (opcode == SUM_SHIFTED_X2_OPCODE) {
		if (Rtype == REALSXP) {
			*summarize_doubles_FUN = sum_shifted_X2_doubles;
		} else {
			*summarize_ints_FUN = sum_shifted_X2_ints;
		}
		double_init[0] = 0.0;
		double_init[1] = shift;
		return;
	}
	if (opcode == SUM_X_X2_OPCODE) {
		if (Rtype == REALSXP) {
			*summarize_doubles_FUN = sum_X_X2_doubles;
		} else {
			*summarize_ints_FUN = sum_X_X2_ints;
		}
		double_init[0] = double_init[1] = 0.0;
		return;
	}
	if (Rtype == REALSXP) {
		switch (opcode) {
		    case MIN_OPCODE:
			*summarize_doubles_FUN = min_doubles;
			double_init[0] = R_PosInf;
			break;
		    case MAX_OPCODE:
			*summarize_doubles_FUN = max_doubles;
			double_init[0] = R_NegInf;
			break;
		    case RANGE_OPCODE:
			*summarize_doubles_FUN = range_doubles;
			double_init[0] = R_PosInf;
			double_init[1] = R_NegInf;
			break;
		}
		return;
	}
	/* NO initial value! */
	switch (opcode) {
	    case MIN_OPCODE:   *summarize_ints_FUN = min_ints;   break;
	    case MAX_OPCODE:   *summarize_ints_FUN = max_ints;   break;
	    case RANGE_OPCODE: *summarize_ints_FUN = range_ints; break;
	}
	return;
}

/* Returns an integer or numeric vector of length 1 or 2. */
static SEXP init2nakedSEXP(int opcode, SEXPTYPE Rtype, void *init, int status)
{
	int *int_init = (int *) init;
	double *double_init = (double *) init;
	SEXP ans;

	if (opcode == ANY_OPCODE || opcode == ALL_OPCODE)
		return ScalarLogical(int_init[0]);

	if (opcode == MIN_OPCODE || opcode == MAX_OPCODE) {
		if (Rtype == REALSXP)
			return ScalarReal(double_init[0]);
		if (status == 0) {
			return ScalarReal(opcode == MIN_OPCODE ? R_PosInf
							       : R_NegInf);
		} else {
			return ScalarInteger(int_init[0]);
		}
	}

	if (opcode == RANGE_OPCODE) {
		if (Rtype == REALSXP) {
			ans = PROTECT(NEW_NUMERIC(2));
			REAL(ans)[0] = double_init[0];
			REAL(ans)[1] = double_init[1];
		} else {
			if (status == 0) {
				ans = PROTECT(NEW_NUMERIC(2));
				REAL(ans)[0] = R_PosInf;
				REAL(ans)[1] = R_NegInf;
			} else {
				ans = PROTECT(NEW_INTEGER(2));
				INTEGER(ans)[0] = int_init[0];
				INTEGER(ans)[1] = int_init[1];
			}
		}
		UNPROTECT(1);
		return ans;
	}

	if (opcode == SUM_X_X2_OPCODE) {
		/* Either 'double_init[0]' and 'double_init[1]' are both set
		   to NA or NaN or none of them is. Furthermore, in the former
		   case, they should both be set to the same kind of NA i.e.
		   both are set to NA or both are set NaN. */
		double init0 = double_init[0], init1 = double_init[1];
		if (Rtype == INTSXP) {
			/* When 'Rtype' is INTSXP, the only kind of NA that
			   can end up in 'init0' is NA_REAL. No NaN. */
			if (R_IsNA(init0)) {
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
		if (DOUBLE_IS_NA(init0)) {
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
	if (Rtype == REALSXP || opcode == SUM_SHIFTED_X2_OPCODE)
		return ScalarReal(double_init[0]);

	if (R_IsNA(double_init[0]))
		return ScalarInteger(NA_INTEGER);
	if (double_init[0] <= INT_MAX && double_init[0] >= -INT_MAX)
		/* Round 'double_init[0]' to the nearest integer. */
		return ScalarInteger((int) (double_init[0] + 0.5));
	return ScalarReal(double_init[0]);
}


/****************************************************************************
 * _init_SummarizeOp()
 * _apply_summarize_op()
 * _make_SEXP_from_summarize_result()
 */

SummarizeOp _init_SummarizeOp(int opcode, SEXPTYPE Rtype,
		int na_rm, double shift, void *init)
{
	SummarizeOp summarize_op;
	SummarizeInts_FUNType summarize_ints_FUN;
	SummarizeDoubles_FUNType summarize_doubles_FUN;

	select_summarize_FUN(opcode, Rtype, shift,
		&summarize_ints_FUN, &summarize_doubles_FUN, init);

	summarize_op.opcode = opcode;
	summarize_op.Rtype = Rtype;
	summarize_op.na_rm = na_rm;
	summarize_op.shift = shift;
	summarize_op.summarize_ints_FUN = summarize_ints_FUN;
	summarize_op.summarize_doubles_FUN = summarize_doubles_FUN;
	return summarize_op;
}

int _apply_summarize_op(const SummarizeOp *summarize_op,
		void *init, const void *x, int n, R_xlen_t *na_rm_count,
		int status)
{
	if (summarize_op->Rtype == INTSXP) {
		status = summarize_op->summarize_ints_FUN(
					init, (const int *) x, n,
					summarize_op->na_rm, na_rm_count,
					status);
	} else {
		status = summarize_op->summarize_doubles_FUN(
					init, (const double *) x, n,
					summarize_op->na_rm, na_rm_count,
					status);
	}
	return status;
}

/* Returns an integer or numeric vector of length 1 or 2.
   If 'na_rm' is TRUE, then the "na_rm_count" attribute is set on the
   returned vector. */
SEXP _make_SEXP_from_summarize_result(int opcode, SEXPTYPE Rtype,
		void *init, int na_rm, R_xlen_t na_rm_count, int status)
{
	SEXP ans, ans_attrib;

	ans = init2nakedSEXP(opcode, Rtype, init, status);
	if (!na_rm)
		return ans;
	PROTECT(ans);
	if (na_rm_count > INT_MAX)
		ans_attrib = ScalarReal((double) na_rm_count);
	else
		ans_attrib = ScalarInteger((int) na_rm_count);
	PROTECT(ans_attrib);
	setAttrib(ans, install("na_rm_count"), ans_attrib);
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
		if (DOUBLE_IS_NA(*x))
			count++;
	return count;
}
static int any_NA_double_elt(const double *x, int n)
{
	for (int i = 0; i < n; i++, x++)
		if (DOUBLE_IS_NA(*x))
			return 1;
	return 0;
}

#define RCOMPLEX_IS_NA(x) (DOUBLE_IS_NA((x)->r) || DOUBLE_IS_NA((x)->i))
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
	SEXPTYPE Rtype;

	Rtype = TYPEOF(x);
	switch (Rtype) {
	    case LGLSXP: case INTSXP:
		if (LENGTH(x) == 1 && INTEGER(x)[0] == NA_INTEGER)
			return 1;
	    break;
	    case REALSXP:
		if (LENGTH(x) == 1 && DOUBLE_IS_NA(REAL(x)[0]))
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

