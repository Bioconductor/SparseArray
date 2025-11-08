#ifndef _SPARSEVEC_H_
#define _SPARSEVEC_H_

#include <Rdefines.h>

#include "Rvector_utils.h"

#include <limits.h>  /* for INT_MAX */


/****************************************************************************
 * SparseVec struct
 */

/* Set 'nzvals' to R_NilValue to represent a lacunar SparseVec. */
typedef struct sparse_vec_t {
	SEXPTYPE Rtype;     /* type of the values in 'nzvals' */
	void *nzvals;       /* NULL or array of nonzero values */
	int *nzoffs;        /* array of offsets for the nonzero values */
	int nzcount;        /* nb of nonzero values */
	int len;            /* vector length (= nzcount + nb of zeros) */
	int na_background;  /* background value is NA instead of zero */
} SparseVec;


/****************************************************************************
 * Some low-level convenience macros and inline functions used by many
 * operations on SparseVec structs.
 */

#define IS_BG_DOUBLE(val, na_background) \
	((na_background) ? R_IsNA(val) : ((val) == double0))

#define IS_BG_CHARSXP(val, na_background) \
	((na_background) ? ((val) == NA_STRING) \
			 : ((val) != NA_STRING && LENGTH(val) == 0))

#define APPEND_TO_NZVALS_NZOFFS(val, off, out_nzvals, out_nzoffs, out_nzcount) \
{									\
	(out_nzvals)[out_nzcount] = (val);				\
	(out_nzoffs)[out_nzcount] = (off);				\
	(out_nzcount)++;						\
}

/* PROPAGATE_NZOFFS is a special value returned by _Arith_sv1_scalar(),
   _Compare_sv1_scalar(), and other functions that take a single input
   SparseVec to indicate that the result of the operation is a sparse vector
   with the same nzoffs as the input ones and with a single nzval shared by
   all the nzoffs.
   IMPORTANT: If this is the case then the function doesn't write anything
   to output buffer 'out_nzoffs' and writes the single shared nzval to
   'out_nzvals[0]'. */
#define PROPAGATE_NZOFFS   -1  /* must be a **negative** int */

/* Maybe move this to Rvector_utils.h */
static inline int int_equal(int x, int y)
{
	return x == y;
}
static inline int double_equal(double x, double y)
{
	return x == y;
}
static inline int Rcomplex_equal(Rcomplex x, Rcomplex y)
{
	return x.r == y.r && x.i == y.i;
}


/****************************************************************************
 * Inline function toSparseVec()
 */

/* 'Rtype' **must** be set to 'TYPEOF(nzvals)' if 'nzvals' is not R_NilValue.
   The only reason we have the 'Rtype' argument is so that we can still store
   the 'Rtype' in the SparseVec even when the supplied 'nzvals' is R_NilValue
   (lacunar case). */
static inline SparseVec toSparseVec(SEXP nzvals, SEXP nzoffs,
		SEXPTYPE Rtype, int len, int na_background)
{
	/* Sanity checks (should never fail). */
	if (!IS_INTEGER(nzoffs))
		goto on_error;
	R_xlen_t nzcount = XLENGTH(nzoffs);
	if (nzcount == 0 || nzcount > INT_MAX)
		goto on_error;

	if (na_background && Rtype == RAWSXP)
		error("SparseArray internal error in toSparseVec():\n"
		      "    NaArray objects of type \"raw\" are not supported");

	SparseVec sv;
	sv.Rtype = Rtype;
	if (nzvals == R_NilValue) {
		sv.nzvals = NULL;
	} else {
		if (Rtype != INTSXP && Rtype != LGLSXP && Rtype != REALSXP &&
		    Rtype != CPLXSXP && Rtype != RAWSXP &&
		    Rtype != STRSXP && Rtype != VECSXP)
			error("SparseArray internal error in toSparseVec():\n"
			      "    type \"%s\" is not supported",
			      type2char(Rtype));
		if (TYPEOF(nzvals) != Rtype)
			error("SparseArray internal error in toSparseVec():\n"
			      "    TYPEOF(nzvals) != Rtype");
		if (XLENGTH(nzvals) != nzcount)
			goto on_error;
		if (IS_STRSXP_OR_VECSXP(Rtype)) {
			sv.nzvals = nzvals;
		} else {
			sv.nzvals = DATAPTR(nzvals);
		}
	}
	sv.nzoffs = INTEGER(nzoffs);
	sv.nzcount = LENGTH(nzoffs);
	sv.len = len;
	sv.na_background = na_background;
	return sv;

    on_error:
	error("SparseArray internal error in toSparseVec():\n"
	      "    supplied 'nzvals' and/or 'nzoffs' "
	      "are invalid or incompatible");
}


/****************************************************************************
 * SparseVec getters as inline functions
 */

static inline SEXPTYPE get_SV_Rtype(const SparseVec *sv)
{
	return sv->Rtype;
}

static inline int get_SV_nzcount(const SparseVec *sv)
{
	return sv->nzcount;
}

static inline const int *get_intSV_nzvals_p(const SparseVec *sv)
{
	return sv->nzvals;
}

static inline const double *get_doubleSV_nzvals_p(const SparseVec *sv)
{
	return sv->nzvals;
}

static inline const Rcomplex *get_RcomplexSV_nzvals_p(const SparseVec *sv)
{
	return sv->nzvals;
}

static inline const Rbyte *get_RbyteSV_nzvals_p(const SparseVec *sv)
{
	return sv->nzvals;
}

static inline SEXP get_characterSV_nzvals_p(const SparseVec *sv)
{
	return sv->nzvals;
}

static inline SEXP get_listSV_nzvals_p(const SparseVec *sv)
{
	return sv->nzvals;
}

static inline Rbyte get_RbyteSV_nzval(const SparseVec *sv, int k)
{
	const Rbyte *nzvals_p = get_RbyteSV_nzvals_p(sv);
	return nzvals_p == NULL ? Rbyte1 : nzvals_p[k];
}

static inline int get_intSV_nzval(const SparseVec *sv, int k)
{
	const int *nzvals_p = get_intSV_nzvals_p(sv);
	return nzvals_p == NULL ? int1 : nzvals_p[k];
}

static inline double get_doubleSV_nzval(const SparseVec *sv, int k)
{
	const double *nzvals_p = get_doubleSV_nzvals_p(sv);
	return nzvals_p == NULL ? double1 : nzvals_p[k];
}

static inline Rcomplex get_RcomplexSV_nzval(const SparseVec *sv, int k)
{
	const Rcomplex *nzvals_p = get_RcomplexSV_nzvals_p(sv);
	return nzvals_p == NULL ? Rcomplex1 : nzvals_p[k];
}

static inline SEXP get_characterSV_nzval(const SparseVec *sv, int k)
{
	SEXP nzvals_p = get_characterSV_nzvals_p(sv);
	return nzvals_p == NULL ? character1 : STRING_ELT(nzvals_p, k);
}

static inline SEXP get_listSV_nzval(const SparseVec *sv, int k)
{
	SEXP nzvals_p = get_listSV_nzvals_p(sv);
	if (nzvals_p == NULL)
		error("SparseArray internal error in get_listSV_nzval():\n"
		      "    lacunar SparseVec of type \"list\" not supported");
	return VECTOR_ELT(nzvals_p, k);
}


/****************************************************************************
 * Inline function next_offset()
 */

static inline int next_offset(
		const int *offs1, int n1,
		const int *offs2, int n2,
		int k1, int k2, int *off)
{
	if (k1 < n1 && k2 < n2) {
		int off1 = offs1[k1];
		int off2 = offs2[k2];
		if (off1 < off2) {
			*off = off1;
			return 1;
		}
		if (off1 > off2) {
			*off = off2;
			return 2;
		}
		*off = off1;  /* same as 'off2' */
		return 3;
	}
	if (k1 < n1) {
		*off = offs1[k1];
		return 1;
	}
	if (k2 < n2) {
		*off = offs2[k2];
		return 2;
	}
	return 0;
}


/****************************************************************************
 * The next_<Ltype>SV_<Rtype>SV_vals() inline functions (11 in total)
 */

static inline int next_RbyteSV_RbyteSV_vals(
	const SparseVec *sv1, const SparseVec *sv2,
	int *k1, int *k2, int *off, Rbyte *val1, Rbyte *val2)
{
	if (sv1->na_background || sv2->na_background)
		error("SparseArray internal error in "
		      "next_RbyteSV_RbyteSV_vals():\n"
		      "    NaArray objects of type \"raw\" are not supported");
	int ret = next_offset(sv1->nzoffs, get_SV_nzcount(sv1),
			      sv2->nzoffs, get_SV_nzcount(sv2),
			      *k1, *k2, off);
	switch (ret) {
	    case 1: {
		*val1 = get_RbyteSV_nzval(sv1, *k1);
		*val2 = Rbyte0;
		(*k1)++;
		break;
	    }
	    case 2: {
		*val1 = Rbyte0;
		*val2 = get_RbyteSV_nzval(sv2, *k2);
		(*k2)++;
		break;
	    }
	    case 3: {
		*val1 = get_RbyteSV_nzval(sv1, *k1);
		*val2 = get_RbyteSV_nzval(sv2, *k2);
		(*k1)++;
		(*k2)++;
		break;
	    }
	}
	return ret;
}

#define DEFINE_next_RbyteSV_RtypeSV_vals_FUN(Rtype)			\
static inline int next_RbyteSV_ ## Rtype ## SV_vals(			\
	const SparseVec *sv1, const SparseVec *sv2,			\
	int *k1, int *k2, int *off, Rbyte *val1, Rtype *val2)		\
{									\
	if (sv1->na_background)						\
		error("SparseArray internal error in "			\
		      "next_RbyteSV_<Rtype>SV_vals():\n"		\
		      "    NaArray objects of type \"raw\" "		\
		      "are not supported");				\
	int ret = next_offset(sv1->nzoffs, get_SV_nzcount(sv1),		\
			      sv2->nzoffs, get_SV_nzcount(sv2),		\
			      *k1, *k2, off);				\
	switch (ret) {							\
	    case 1: {							\
		*val1 = get_RbyteSV_nzval(sv1, *k1);			\
		*val2 = sv2->na_background ? Rtype ## NA : Rtype ## 0;	\
		(*k1)++;						\
		break;							\
	    }								\
	    case 2: {							\
		*val1 = Rbyte0;						\
		*val2 = get_ ## Rtype ## SV_nzval(sv2, *k2);		\
		(*k2)++;						\
		break;							\
	    }								\
	    case 3: {							\
		*val1 = get_RbyteSV_nzval(sv1, *k1);			\
		*val2 = get_ ## Rtype ## SV_nzval(sv2, *k2);		\
		(*k1)++;						\
		(*k2)++;						\
		break;							\
	    }								\
	}								\
	return ret;							\
}

#define DEFINE_next_LtypeSV_RtypeSV_vals_FUN(Ltype, Rtype)		\
static inline int next_ ## Ltype ## SV_ ## Rtype ## SV_vals(		\
	const SparseVec *sv1, const SparseVec *sv2,			\
	int *k1, int *k2, int *off, Ltype *val1, Rtype *val2)		\
{									\
	int ret = next_offset(sv1->nzoffs, get_SV_nzcount(sv1),		\
			      sv2->nzoffs, get_SV_nzcount(sv2),		\
			      *k1, *k2, off);				\
	switch (ret) {							\
	    case 1: {							\
		*val1 = get_ ## Ltype ## SV_nzval(sv1, *k1);		\
		*val2 = sv2->na_background ? Rtype ## NA : Rtype ## 0;	\
		(*k1)++;						\
		break;							\
	    }								\
	    case 2: {							\
		*val1 = sv1->na_background ? Ltype ## NA : Ltype ## 0;	\
		*val2 = get_ ## Rtype ## SV_nzval(sv2, *k2);		\
		(*k2)++;						\
		break;							\
	    }								\
	    case 3: {							\
		*val1 = get_ ## Ltype ## SV_nzval(sv1, *k1);		\
		*val2 = get_ ## Rtype ## SV_nzval(sv2, *k2);		\
		(*k1)++;						\
		(*k2)++;						\
		break;							\
	    }								\
	}								\
	return ret;							\
}

DEFINE_next_RbyteSV_RtypeSV_vals_FUN(int)
DEFINE_next_RbyteSV_RtypeSV_vals_FUN(double)
DEFINE_next_RbyteSV_RtypeSV_vals_FUN(Rcomplex)
DEFINE_next_LtypeSV_RtypeSV_vals_FUN(int, int)
DEFINE_next_LtypeSV_RtypeSV_vals_FUN(int, double)
DEFINE_next_LtypeSV_RtypeSV_vals_FUN(int, Rcomplex)
DEFINE_next_LtypeSV_RtypeSV_vals_FUN(double, int)
DEFINE_next_LtypeSV_RtypeSV_vals_FUN(double, double)
DEFINE_next_LtypeSV_RtypeSV_vals_FUN(double, Rcomplex)
DEFINE_next_LtypeSV_RtypeSV_vals_FUN(Rcomplex, Rcomplex)


/****************************************************************************
 * Function prototypes
 */

SparseVec _alloc_buf_SparseVec(
	SEXPTYPE Rtype,
	int len,
	int na_background
);

void _write_Rsubvec_to_SV(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	const int *offs,
	int n,
	SparseVec *out_sv
);

void _expand_intSV(
	const SparseVec *sv,
	int *out,
	int set_background
);

void _expand_doubleSV(
	const SparseVec *sv,
	double *out,
	int set_background
);

#endif  /* _SPARSEVEC_H_ */

