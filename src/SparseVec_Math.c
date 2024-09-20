/****************************************************************************
 *                   'Math' operations on sparse vectors                    *
 ****************************************************************************/
#include "SparseVec_Math.h"

#include <math.h>   /* for fabs(), sqrt(), floor(), ceil(), trunc(),
		       log(), log10(), log2(), log1p(), exp(), expm1(),
		       sin(), asin(), sinh(), asinh(),
		       cos(), acos(), cosh(), acosh(),
		       tan(), atan(), tanh(), atanh() */
#include <Rmath.h>  /* for sign(), sinpi(), cospi(), tanpi(),
		       gammafn(), lgammafn(), digamma(), trigamma(),
		       fround(), fprec() */

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

static inline double Rlog_double	FUNDEF_Rmath_double(log)
static inline double Rlog10_double	FUNDEF_Rmath_double(log10)
static inline double Rlog2_double	FUNDEF_Rmath_double(log2)
static inline double Rlog1p_double	FUNDEF_Rmath_double(log1p)
static inline double Rexp_double	FUNDEF_Rmath_double(exp)
static inline double Rexpm1_double	FUNDEF_Rmath_double(expm1)

static inline double Rsin_double	FUNDEF_Rmath_double(sin)
static inline double Rasin_double	FUNDEF_Rmath_double(asin)
static inline double Rsinh_double	FUNDEF_Rmath_double(sinh)
static inline double Rasinh_double	FUNDEF_Rmath_double(asinh)
static inline double Rsinpi_double	FUNDEF_Rmath_double(sinpi)

static inline double Rcos_double	FUNDEF_Rmath_double(cos)
static inline double Racos_double	FUNDEF_Rmath_double(acos)
static inline double Rcosh_double	FUNDEF_Rmath_double(cosh)
static inline double Racosh_double	FUNDEF_Rmath_double(acosh)
static inline double Rcospi_double	FUNDEF_Rmath_double(cospi)

static inline double Rtan_double	FUNDEF_Rmath_double(tan)
static inline double Ratan_double	FUNDEF_Rmath_double(atan)
static inline double Rtanh_double	FUNDEF_Rmath_double(tanh)
static inline double Ratanh_double	FUNDEF_Rmath_double(atanh)
static inline double Rtanpi_double	FUNDEF_Rmath_double(tanpi)

static inline double Rgamma_double	FUNDEF_Rmath_double(gammafn)
static inline double Rlgamma_double	FUNDEF_Rmath_double(lgammafn)
static inline double Rdigamma_double	FUNDEF_Rmath_double(digamma)
static inline double Rtrigamma_double	FUNDEF_Rmath_double(trigamma)

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

	if (strcmp(op, "log") == 0)        /* does NOT propagate zeros */
		return Rlog_double;
	if (strcmp(op, "log10") == 0)      /* does NOT propagate zeros */
		return Rlog10_double;
	if (strcmp(op, "log2") == 0)       /* does NOT propagate zeros */
		return Rlog2_double;
	if (strcmp(op, "log1p") == 0)
		return Rlog1p_double;
	if (strcmp(op, "exp") == 0)        /* does NOT propagate zeros */
		return Rexp_double;
	if (strcmp(op, "expm1") == 0)
		return Rexpm1_double;

	if (strcmp(op, "sin") == 0)
		return Rsin_double;
	if (strcmp(op, "asin") == 0)
		return Rasin_double;
	if (strcmp(op, "sinh") == 0)
		return Rsinh_double;
	if (strcmp(op, "asinh") == 0)
		return Rasinh_double;
	if (strcmp(op, "sinpi") == 0)
		return Rsinpi_double;

	if (strcmp(op, "cos") == 0)        /* does NOT propagate zeros */
		return Rcos_double;
	if (strcmp(op, "acos") == 0)       /* does NOT propagate zeros */
		return Racos_double;
	if (strcmp(op, "cosh") == 0)       /* does NOT propagate zeros */
		return Rcosh_double;
	if (strcmp(op, "acosh") == 0)      /* does NOT propagate zeros */
		return Racosh_double;
	if (strcmp(op, "cospi") == 0)      /* does NOT propagate zeros */
		return Rcospi_double;

	if (strcmp(op, "tan") == 0)
		return Rtan_double;
	if (strcmp(op, "atan") == 0)
		return Ratan_double;
	if (strcmp(op, "tanh") == 0)
		return Rtanh_double;
	if (strcmp(op, "atanh") == 0)
		return Ratanh_double;
	if (strcmp(op, "tanpi") == 0)
		return Rtanpi_double;

	if (strcmp(op, "gamma") == 0)     /* does NOT propagate zeros */
		return Rgamma_double;
	if (strcmp(op, "lgamma") == 0)    /* does NOT propagate zeros */
		return Rlgamma_double;
	if (strcmp(op, "digamma") == 0)   /* does NOT propagate zeros */
		return Rdigamma_double;
	if (strcmp(op, "trigamma") == 0)  /* does NOT propagate zeros */
		return Rtrigamma_double;

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

void _Math_doubleSV(MathFUN fun, const SparseVec *sv, double digits,
		SparseVec *out_sv, int *newNaNs)
{
        if (out_sv->len != sv->len)
		error("SparseArray internal error in "
		      "_Math_doubleSV():\n"
		      "    'sv' and 'out_sv' are incompatible");
	set_NaNs_produced_flag(0);
	digits0 = digits;
	double *out_nzvals = (double *) out_sv->nzvals;
	out_sv->nzcount = 0;
	const double *nzvals_p = get_doubleSV_nzvals_p(sv);
	if (nzvals_p == NULL) {  /* lacunar SparseVec */
		double out_val = fun(1.0);
		if (IS_BACKGROUND_VAL(out_val, out_sv->na_background))
			return;
		out_nzvals[0] = out_val;
		out_sv->nzcount = PROPAGATE_NZOFFS;
		return;
	}
	/* regular SparseVec */
	int nzcount = get_SV_nzcount(sv);
	for (int k = 0; k < nzcount; k++) {
		double out_val = fun(nzvals_p[k]);
		if (IS_BACKGROUND_VAL(out_val, out_sv->na_background))
			continue;
                APPEND_TO_NZVALS_NZOFFS(out_val, sv->nzoffs[k],
				out_nzvals, out_sv->nzoffs, out_sv->nzcount);
	}
	if (get_NaNs_produced_flag())
		*newNaNs = 1;
	return;
}

