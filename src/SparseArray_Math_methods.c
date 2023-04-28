/****************************************************************************
 *                   Math methods for SparseArray objects                   *
 ****************************************************************************/
#include "SparseArray_Math_methods.h"

#include "leaf_vector_utils.h"

#include <math.h>   /* for fabs(), sqrt(), floor(), ceil(), trunc(),
		       log1p(), expm1(), sin(), asin(), tan(), atan(),
		       sinh(), asinh(), tanh(), atanh() */
#include <Rmath.h>  /* for sign(), sinpi(), tanpi(), fround(), fprec() */

static int NaNs_produced_flag;

static inline void set_NaNs_produced_flag(int flag)
{
	NaNs_produced_flag = flag;
	return;
}

static inline int get_NaNs_produced_flag()
{
	return NaNs_produced_flag;
}


/****************************************************************************
 * select_double2double_FUN()
 */

#define	ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(funname)(double x)	\
{								\
	double v = funname(x);					\
	if (ISNAN(v) && !ISNAN(x))				\
		set_NaNs_produced_flag(1);			\
	return v;						\
}

static inline double Rabs_double ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(fabs)
static inline double Rsign_double ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(sign)
static inline double Rsqrt_double ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(sqrt)
static inline double Rfloor_double ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(floor)
static inline double Rceiling_double ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(ceil)
static inline double Rtrunc_double ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(trunc)
static inline double Rlog1p_double ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(log1p)
static inline double Rexpm1_double ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(expm1)
static inline double Rsin_double ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(sin)
static inline double Rsinpi_double ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(sinpi)
static inline double Rasin_double ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(asin)
static inline double Rtan_double ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(tan)
static inline double Rtanpi_double ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(tanpi)
static inline double Ratan_double ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(atan)
static inline double Rsinh_double ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(sinh)
static inline double Rasinh_double ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(asinh)
static inline double Rtanh_double ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(tanh)
static inline double Ratanh_double ARGS_AND_BODY_OF_MATH_DOUBLE_FUN(atanh)

static double digits0;
static inline double Rround_double(double x) { return fround(x, digits0); }
static inline double Rsignif_double(double x) { return fprec(x, digits0); }

static double (*select_double2double_FUN(const char *op))(double)
{
	/* 'Math' group */
	if (strcmp(op, "abs") == 0)
		return Rabs_double;
	if (strcmp(op, "sign") == 0)
		return Rsign_double;
	if (strcmp(op, "sqrt") == 0)
		return Rsqrt_double;
	if (strcmp(op, "floor") == 0)
		return Rfloor_double;
	if (strcmp(op, "ceiling") == 0)
		return Rceiling_double;
	if (strcmp(op, "trunc") == 0)
		return Rtrunc_double;
	if (strcmp(op, "log1p") == 0)
		return Rlog1p_double;
	if (strcmp(op, "expm1") == 0)
		return Rexpm1_double;
	if (strcmp(op, "sin") == 0)
		return Rsin_double;
	if (strcmp(op, "sinpi") == 0)
		return Rsinpi_double;
	if (strcmp(op, "asin") == 0)
		return Rasin_double;
	if (strcmp(op, "tan") == 0)
		return Rtan_double;
	if (strcmp(op, "tanpi") == 0)
		return Rtanpi_double;
	if (strcmp(op, "atan") == 0)
		return Ratan_double;
	if (strcmp(op, "sinh") == 0)
		return Rsinh_double;
	if (strcmp(op, "asinh") == 0)
		return Rasinh_double;
	if (strcmp(op, "tanh") == 0)
		return Rtanh_double;
	if (strcmp(op, "atanh") == 0)
		return Ratanh_double;

	/* 'Math2' group */
	if (strcmp(op, "round") == 0)
		return Rround_double;
	if (strcmp(op, "signif") == 0)
		return Rsignif_double;

	error("SparseArray internal error in select_double2double_FUN():\n"
	      "    unsupported 'Math' or 'Math2' function: \"%s\"", op);
	return NULL;  /* will never reach this */
}


/****************************************************************************
 * C_Math_SVT() and C_Math_SVT2()
 */

static SEXP REC_Math_SVT(SEXP SVT, const int *dims, int ndim,
			 apply_2double_FUNS *funs,
			 int *offs_buf, double *vals_buf)
{
	int ans_len, is_empty, i;
	SEXP ans, ans_elt, subSVT;

	if (SVT == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		return _lv_apply_to_REALSXP(SVT, funs, offs_buf, vals_buf);
	}

	/* 'SVT' is a list. */
	ans_len = dims[ndim - 1];
	ans = PROTECT(NEW_LIST(ans_len));
	is_empty = 1;
	for (i = 0; i < ans_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		ans_elt = REC_Math_SVT(subSVT, dims, ndim - 1,
				       funs, offs_buf, vals_buf);
		if (ans_elt != R_NilValue) {
			PROTECT(ans_elt);
			SET_VECTOR_ELT(ans, i, ans_elt);
			UNPROTECT(1);
			is_empty = 0;
		}
	}
	UNPROTECT(1);
	return is_empty ? R_NilValue : ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_Math_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP op)
{
	SEXPTYPE x_Rtype;
	apply_2double_FUNS funs;
	int *offs_buf;
	double *vals_buf;
	SEXP ans;

	x_Rtype = _get_Rtype_from_Rstring(x_type);
	if (x_Rtype == 0)
		error("SparseArray internal error in C_Math_SVT():\n"
		      "    invalid 'x_type' value");

	if (!IS_CHARACTER(op) || LENGTH(op) != 1)
		error("SparseArray internal error in C_Math_SVT():\n"
		      "    'op' must be a single string");
	op = STRING_ELT(op, 0);
	if (op == NA_STRING)
		error("SparseArray internal error in C_Math_SVT():\n"
		      "    'op' cannot be NA");

	funs.Rbyte2double_FUN = NULL;
	funs.int2double_FUN = NULL;
	funs.double2double_FUN = select_double2double_FUN(CHAR(op));
	funs.Rcomplex2double_FUN = NULL;

	offs_buf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	vals_buf = (double *) R_alloc(INTEGER(x_dim)[0], sizeof(double));
	set_NaNs_produced_flag(0);
	ans = REC_Math_SVT(x_SVT, INTEGER(x_dim), LENGTH(x_dim),
			   &funs, offs_buf, vals_buf);
	if (ans != R_NilValue && get_NaNs_produced_flag()) {
		PROTECT(ans);
		warning("NaNs produced");
		UNPROTECT(1);
	}
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_Math2_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP op, SEXP digits)
{
	SEXPTYPE x_Rtype;
	apply_2double_FUNS funs;
	int *offs_buf;
	double *vals_buf;

	x_Rtype = _get_Rtype_from_Rstring(x_type);
	if (x_Rtype == 0)
		error("SparseArray internal error in C_Math_SVT():\n"
		      "    invalid 'x_type' value");

	if (!IS_CHARACTER(op) || LENGTH(op) != 1)
		error("SparseArray internal error in C_Math_SVT():\n"
		      "    'op' must be a single string");
	op = STRING_ELT(op, 0);
	if (op == NA_STRING)
		error("SparseArray internal error in C_Math_SVT():\n"
		      "    'op' cannot be NA");

	funs.Rbyte2double_FUN = NULL;
	funs.int2double_FUN = NULL;
	funs.double2double_FUN = select_double2double_FUN(CHAR(op));
	funs.Rcomplex2double_FUN = NULL;
	digits0 = REAL(digits)[0];

	offs_buf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	vals_buf = (double *) R_alloc(INTEGER(x_dim)[0], sizeof(double));
	return REC_Math_SVT(x_SVT, INTEGER(x_dim), LENGTH(x_dim),
			    &funs, offs_buf, vals_buf);
}

