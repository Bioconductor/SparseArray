#include "coerceVector2.h"


/****************************************************************************
 * _CoercionWarning()
 *
 * Unfortunately CoercionWarning() is already implemented in
 * R/src/main/coerce.c but it cannot be called from package code.
 * Yeah it's a shame that we need to reimplement it here!
 */

/* Coercion warnings will be OR'ed : */
#define	WARN_NA     1
#define	WARN_INT_NA 2
#define	WARN_IMAG   4
#define	WARN_RAW    8

void _CoercionWarning(int warn)
{
	if (warn & WARN_NA)
		warning("NAs introduced by coercion");
	if (warn & WARN_INT_NA)
		warning("NAs introduced by coercion to integer range");
	if (warn & WARN_IMAG)
		warning("imaginary parts discarded in coercion");
	if (warn & WARN_RAW)
		warning("out-of-range values treated as 0 in coercion to raw");
	return;
}

#define COERCE_ERROR_FMT "coercion from type \"%s\" to type \"%s\" " \
			 "is not supported"
#define COERCE_ERROR(from, to) \
	error(COERCE_ERROR_FMT, type2char(from), type2char(to))


/****************************************************************************
 * coerceToLogical()
 *
 * Unfortunately coerceToLogical() is already implemented in
 * R/src/main/coerce.c but it cannot be called from package code.
 * Yeah it's a shame that we need to reimplement the function here!
 */

static inline int LogicalFromInteger(int x)
{
	return x == NA_INTEGER ? NA_LOGICAL : x != 0;
}

static inline int LogicalFromReal(double x)
{
	return ISNAN(x) ? NA_LOGICAL : x != 0.0;
}

static inline int LogicalFromComplex(Rcomplex x)
{
	return ISNAN(x.r) || ISNAN(x.i) ? NA_LOGICAL : x.r != 0 || x.i != 0;
}

static inline int LogicalFromString(SEXP x)
{
    if (x != R_NaString) {
        if (StringTrue(CHAR(x))) return 1;
        if (StringFalse(CHAR(x))) return 0;
    }
    return NA_LOGICAL;
}

static SEXP coerceToLogical(SEXP v)
{
	R_xlen_t ans_len, i;
	SEXP ans;
	int *pa;

	ans_len = XLENGTH(v);
	ans = PROTECT(NEW_LOGICAL(ans_len));
	pa = LOGICAL(ans);
	switch (TYPEOF(v)) {
	    case INTSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = LogicalFromInteger(INTEGER(v)[i]);
		break;
	    case REALSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = LogicalFromReal(REAL(v)[i]);
		break;
	    case CPLXSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = LogicalFromComplex(COMPLEX(v)[i]);
		break;
	    case STRSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = LogicalFromString(STRING_ELT(v, i));
		break;
	    case RAWSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = RAW(v)[i] != 0;
		break;
	    default:
		COERCE_ERROR(TYPEOF(v), LGLSXP);
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * coerceToInteger()
 *
 * Unfortunately coerceToInteger() is already implemented in
 * R/src/main/coerce.c but it cannot be called from package code.
 * Also the original coerceToInteger() does not support the 'warn' argument.
 * Sad!
 */

static inline int IntegerFromLogical(int x, int *warn)
{
	return x == NA_LOGICAL ? NA_INTEGER : x;
}

static inline int IntegerFromReal(double x, int *warn)
{
	if (ISNAN(x))
		return NA_INTEGER;
	if (x >= INT_MAX + 1.0 || x <= INT_MIN ) {
		*warn |= WARN_INT_NA;
		return NA_INTEGER;
	}
	return (int) x;
}

static inline int IntegerFromComplex(Rcomplex x, int *warn)
{
	if (ISNAN(x.r) || ISNAN(x.i))
		return NA_INTEGER;
	if (x.r > INT_MAX + 1.0 || x.r <= INT_MIN ) {
		*warn |= WARN_INT_NA;
		return NA_INTEGER;
	}
	if (x.i != 0.0)
		*warn |= WARN_IMAG;
	return (int) x.r;
}

static inline int IntegerFromString(SEXP x, int *warn)
{
    double xdouble;
    char *endp;
    if (x != R_NaString && !isBlankString(CHAR(x))) { /* ASCII */
        xdouble = R_strtod(CHAR(x), &endp); /* ASCII */
        if (isBlankString(endp)) {
            // behave the same as IntegerFromReal() etc:
            if (xdouble >= INT_MAX+1. || xdouble <= INT_MIN ) {
                *warn |= WARN_INT_NA;
                return NA_INTEGER;
            }
            else
                return (int) xdouble;
        }
        else *warn |= WARN_NA;
    }
    return NA_INTEGER;
}

static SEXP coerceToInteger(SEXP v, int *warn)
{
	R_xlen_t ans_len, i;
	SEXP ans;
	int *pa;

	ans_len = XLENGTH(v);
	ans = PROTECT(NEW_INTEGER(ans_len));
	pa = INTEGER(ans);
	switch (TYPEOF(v)) {
	    case LGLSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = IntegerFromLogical(LOGICAL(v)[i], warn);
		break;
	    case REALSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = IntegerFromReal(REAL(v)[i], warn);
		break;
	    case CPLXSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = IntegerFromComplex(COMPLEX(v)[i], warn);
		break;
	    case STRSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = IntegerFromString(STRING_ELT(v, i), warn);
		break;
	    case RAWSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = (int) RAW(v)[i];
		break;
	    default:
		COERCE_ERROR(TYPEOF(v), INTSXP);
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * coerceToReal()
 *
 * Unfortunately coerceToReal() is already implemented in
 * R/src/main/coerce.c but it cannot be called from package code.
 * Also the original coerceToReal() does not support the 'warn' argument.
 * A real shame!
 */

static inline double RealFromLogical(int x, int *warn)
{
	return x == NA_LOGICAL ? NA_REAL : x;
}

static inline double RealFromInteger(int x, int *warn)
{
	return x == NA_INTEGER ? NA_REAL : x;
}

static inline double RealFromComplex(Rcomplex x, int *warn)
{
	if (ISNAN(x.r) || ISNAN(x.i))
		return NA_REAL;
	if (x.i != 0.0)
		*warn |= WARN_IMAG;
	return x.r;
}

static inline double RealFromString(SEXP x, int *warn)
{
    double xdouble;
    char *endp;
    if (x != R_NaString && !isBlankString(CHAR(x))) { /* ASCII */
        xdouble = R_strtod(CHAR(x), &endp); /* ASCII */
        if (isBlankString(endp))
            return xdouble;
        else
            *warn |= WARN_NA;
    }
    return NA_REAL;
}

static SEXP coerceToReal(SEXP v, int *warn)
{
	R_xlen_t ans_len, i;
	SEXP ans;
	double *pa;

	ans_len = XLENGTH(v);
	ans = PROTECT(NEW_NUMERIC(ans_len));
	pa = REAL(ans);
	switch (TYPEOF(v)) {
	    case LGLSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = RealFromLogical(LOGICAL(v)[i], warn);
		break;
	    case INTSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = RealFromInteger(INTEGER(v)[i], warn);
		break;
	    case CPLXSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = RealFromComplex(COMPLEX(v)[i], warn);
		break;
	    case STRSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = RealFromString(STRING_ELT(v, i), warn);
		break;
	    case RAWSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = (double) RAW(v)[i];
		break;
	    default:
		COERCE_ERROR(TYPEOF(v), REALSXP);
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * coerceToComplex()
 *
 * Really sad!
 */

static inline Rcomplex ComplexFromLogical(int x, int *warn)
{
	Rcomplex z;

	if (x == NA_LOGICAL) {
		z.r = NA_REAL;
		z.i = NA_REAL;
	} else {
		z.r = x;
		z.i = 0.0;
	}
	return z;
}

static inline Rcomplex ComplexFromInteger(int x, int *warn)
{
	Rcomplex z;

	if (x == NA_INTEGER) {
		z.r = NA_REAL;
		z.i = NA_REAL;
	} else {
		z.r = x;
		z.i = 0.0;
	}
	return z;
}

static inline Rcomplex ComplexFromString(SEXP x, int *warn)
{
    double xr, xi;
    Rcomplex z;
    const char *xx = CHAR(x); /* ASCII */
    char *endp;

    z.r = z.i = NA_REAL;
    if (x != R_NaString && !isBlankString(xx)) {
        xr = R_strtod(xx, &endp);
        if (isBlankString(endp)) {
            z.r = xr;
            z.i = 0.0;
        }
        else if (*endp == '+' || *endp == '-') {
            xi = R_strtod(endp, &endp);
            if (*endp++ == 'i' && isBlankString(endp)) {
                z.r = xr;
                z.i = xi;
            }
            else *warn |= WARN_NA;
        }
        else *warn |= WARN_NA;
    }
    return z;
}

static SEXP coerceToComplex(SEXP v, int *warn)
{
	R_xlen_t ans_len, i;
	SEXP ans;
	Rcomplex *pa;

	ans_len = XLENGTH(v);
	ans = PROTECT(NEW_COMPLEX(ans_len));
	pa = COMPLEX(ans);
	switch (TYPEOF(v)) {
	    case LGLSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = ComplexFromLogical(LOGICAL(v)[i], warn);
		break;
	    case INTSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = ComplexFromInteger(INTEGER(v)[i], warn);
		break;
	    case REALSXP:
		for (i = 0; i < ans_len; i++) {
		    pa->r = REAL(v)[i];
		    pa->i = 0.0;
		    pa++;
		}
		break;
	    case STRSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = ComplexFromString(STRING_ELT(v, i), warn);
		break;
	    case RAWSXP:
		for (i = 0; i < ans_len; i++) {
		    pa->r = (double) RAW(v)[i];
		    pa->i = 0.0;
		    pa++;
		}
		break;
	    default:
		COERCE_ERROR(TYPEOF(v), CPLXSXP);
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * coerceToRaw()
 *
 * Very very sad!
 */

static inline Rbyte RawFromLogical(int x, int *warn)
{
	if (x == NA_LOGICAL) {
		*warn |= WARN_RAW;
		return 0;
	}
	return (Rbyte) x;
}

static inline Rbyte RawFromInteger(int x, int *warn)
{
	if (x == NA_INTEGER || x < 0 || x > 255) {
		*warn |= WARN_RAW;
		return 0;
	}
	return (Rbyte) x;
}

static inline Rbyte RawFromReal(double x, int *warn)
{
	int tmp;

	if (ISNAN(x) || x <= -1.0 || x >= 256.0) {
		*warn |= WARN_RAW;
		return 0;
	}
	tmp = (int) x;
	return (Rbyte) tmp;
}

static inline Rbyte RawFromComplex(Rcomplex x, int *warn)
{
	int tmp;

        if (ISNAN(x.r) || ISNAN(x.i) || x.r <= -1.0 || x.r >= 256.0) {
		*warn |= WARN_RAW;
		return 0;
	}
	tmp = (int) x.r;
	if (x.i != 0.0)
		*warn |= WARN_IMAG;
	return (Rbyte) tmp;
}

static inline Rbyte RawFromString(SEXP x, int *warn)
{
    double xdouble;
    char *endp;
    if (x != R_NaString && !isBlankString(CHAR(x))) { /* ASCII */
        xdouble = R_strtod(CHAR(x), &endp); /* ASCII */
        if (isBlankString(endp)) {
            int tmp = (int) xdouble;
            if (tmp < 0 || tmp > 255) {
                *warn |= WARN_RAW;
                return 0;
            }
            return (Rbyte) tmp;
        }
    }
    *warn |= WARN_RAW;
    return 0;
}

static SEXP coerceToRaw(SEXP v, int *warn)
{
	R_xlen_t ans_len, i;
	SEXP ans;
	Rbyte *pa;

	ans_len = XLENGTH(v);
	ans = PROTECT(NEW_RAW(ans_len));
	pa = RAW(ans);
	switch (TYPEOF(v)) {
	    case LGLSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = RawFromLogical(LOGICAL(v)[i], warn);
		break;
	    case INTSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = RawFromInteger(INTEGER(v)[i], warn);
		break;
	    case REALSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = RawFromReal(REAL(v)[i], warn);
		break;
	    case CPLXSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = RawFromComplex(COMPLEX(v)[i], warn);
		break;
	    case STRSXP:
		for (i = 0; i < ans_len; i++)
		    *(pa++) = RawFromString(STRING_ELT(v, i), warn);
		break;
	    default:
		COERCE_ERROR(TYPEOF(v), RAWSXP);
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * _coerceVector2()
 *
 * Like coerceVector() in R/src/main/coerce.c but with the 'warn' argument.
 */

SEXP _coerceVector2(SEXP v, SEXPTYPE type, int *warn)
{
	if (TYPEOF(v) == VECSXP) {
		/* Coercion from list to atomic vector.
		   Not a common use case and typically quite inefficient.
		   Supported for completeness only.
		   We just use coerceVector() to handle this case, which will
		   call coerceVectorList() from R/src/main/coerce.c.
		   IMPORTANT NOTE: Will possibly issue tons of warnings e.g.
		   one per list element in 'v' in the almost-worst-case
		   scenario! */
		return coerceVector(v, type);
	}
	switch (type) {
	    case LGLSXP:  return coerceToLogical(v);
	    case INTSXP:  return coerceToInteger(v, warn);
	    case REALSXP: return coerceToReal(v, warn);
	    case CPLXSXP: return coerceToComplex(v, warn);
	    case RAWSXP:  return coerceToRaw(v, warn);
	}
	/* We use coerceVector() to handle the remaining cases (type==STRSXP
	   and type==VECSXP). This is ok because these cases never issue
	   warnings (see coerceToString() and coerceToVectorList() in
	   R/src/main/coerce.c). */
	return coerceVector(v, type);
}

