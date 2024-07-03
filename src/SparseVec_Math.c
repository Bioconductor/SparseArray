/****************************************************************************
 *                   'Math' operations on sparse vectors                    *
 ****************************************************************************/
#include "SparseVec_Math.h"

#include "Rvector_utils.h"

#include <math.h>   /* for fabs(), sqrt(), floor(), ceil(), trunc(),
		       log1p(), expm1(), sin(), asin(), tan(), atan(),
		       sinh(), asinh(), tanh(), atanh() */
#include <Rmath.h>  /* for sign(), sinpi(), tanpi(), fround(), fprec() */

#include <string.h> /* for strcmp() */


static int NaNs_produced_flag;

static inline void set_NaNs_produced_flag(int flag)
{
	NaNs_produced_flag = flag;
	return;
}

static inline int get_NaNs_produced_flag(void)
{
	return NaNs_produced_flag;
}


/****************************************************************************
 * _get_MathFUN()
 */

#define	FUNDEF_Rmath_double(funname)		\
	(double x)				\
{						\
	double v = funname(x);			\
	if (ISNAN(v) && !ISNAN(x))		\
		set_NaNs_produced_flag(1);	\
	return v;				\
}

static inline double Rabs_double	FUNDEF_Rmath_double(fabs)
static inline double Rsign_double	FUNDEF_Rmath_double(sign)
static inline double Rsqrt_double	FUNDEF_Rmath_double(sqrt)
static inline double Rfloor_double	FUNDEF_Rmath_double(floor)
static inline double Rceiling_double	FUNDEF_Rmath_double(ceil)
static inline double Rtrunc_double	FUNDEF_Rmath_double(trunc)
static inline double Rlog1p_double	FUNDEF_Rmath_double(log1p)
static inline double Rexpm1_double	FUNDEF_Rmath_double(expm1)
static inline double Rsin_double	FUNDEF_Rmath_double(sin)
static inline double Rsinpi_double	FUNDEF_Rmath_double(sinpi)
static inline double Rasin_double	FUNDEF_Rmath_double(asin)
static inline double Rtan_double	FUNDEF_Rmath_double(tan)
static inline double Rtanpi_double	FUNDEF_Rmath_double(tanpi)
static inline double Ratan_double	FUNDEF_Rmath_double(atan)
static inline double Rsinh_double	FUNDEF_Rmath_double(sinh)
static inline double Rasinh_double	FUNDEF_Rmath_double(asinh)
static inline double Rtanh_double	FUNDEF_Rmath_double(tanh)
static inline double Ratanh_double	FUNDEF_Rmath_double(atanh)

static double digits0;  /* yes, double! don't ask me why */
static inline double Rround_double(double x) { return fround(x, digits0); }
static inline double Rsignif_double(double x) { return fprec(x, digits0); }

MathFUN _get_MathFUN(const char *op)
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

	error("SparseArray internal error in _get_MathFUN():\n"
	      "    unsupported 'Math' or 'Math2' function: \"%s\"", op);
	return NULL;  /* will never reach this */
}


/****************************************************************************
 * _Math_doubleSV()
 */

int _Math_doubleSV(MathFUN fun, const SparseVec *sv, double digits,
		double *out_nzvals, int *out_nzoffs, int *newNaNs)
{
	set_NaNs_produced_flag(0);
	digits0 = digits;
	if (sv->nzvals == R_NilValue) {  /* lacunar SparseVec */
		double v = fun(1.0);
		if (v == double0)
			return 0;
		out_nzvals[0] = v;
		return PROPAGATE_NZOFFS;
	}
	/* regular SparseVec */
	const double *nzvals_p = get_doubleSV_nzvals_p(sv);
	int nzcount = get_SV_nzcount(sv);
	int out_nzcount = 0;
	for (int k = 0; k < nzcount; k++) {
		double v = fun(nzvals_p[k]);
		if (v != double0) {
			out_nzvals[out_nzcount] = v;
			out_nzoffs[out_nzcount] = sv->nzoffs[k];
			out_nzcount++;
		}
	}
	if (get_NaNs_produced_flag())
		*newNaNs = 1;
	return out_nzcount;
}

