#ifndef _SPARSEVEC_H_
#define _SPARSEVEC_H_

#include <Rdefines.h>

#include "Rvector_utils.h"

#include <limits.h>  /* for INT_MAX */

/* Set 'nzvals' to R_NilValue to represent a lacunar SparseVec. */
typedef struct sparse_vec_t {
	SEXP nzvals;         /* vector type (atomic or list) or R_NilValue */
	const int *nzoffs;   /* offsets of nonzero values */
	int nzcount;         /* nb of nonzero values */
	SEXPTYPE Rtype;      /* = TYPEOF(nzvals) if nzvals != R_NilValue */
	int len;             /* vector length (= nzcount + nb of zeros) */
} SparseVec;

/* 'Rtype' **must** be set to 'TYPEOF(nzvals)' if 'nzvals' is not R_NilValue.
   The only reason we have the 'Rtype' argument (and SparseVec has the 'Rtype'
   member) is so that we can still store the 'Rtype' in the SparseVec even
   when the supplied 'nzvals' is R_NilValue (lacunar case). */
static inline SparseVec toSparseVec(SEXP nzvals, SEXP nzoffs,
				    SEXPTYPE Rtype, int len)
{
	/* Sanity checks (should never fail). */
	if (!IS_INTEGER(nzoffs))
		goto on_error;
	R_xlen_t nzcount = XLENGTH(nzoffs);
	if (nzcount == 0 || nzcount > INT_MAX)
		goto on_error;
	if (nzvals != R_NilValue) {
		if (TYPEOF(nzvals) != Rtype)
			error("SparseArray internal error in toSparseVec():\n"
			      "    TYPEOF(nzvals) != Rtype");
		if (XLENGTH(nzvals) != nzcount)
			goto on_error;
	}

	SparseVec sv;
	sv.nzvals = nzvals;
	sv.nzoffs = INTEGER(nzoffs);
	sv.nzcount = LENGTH(nzoffs);
	sv.Rtype = Rtype;
	sv.len = len;
	return sv;

    on_error:
	error("SparseArray internal error in toSparseVec():\n"
	      "    supplied 'nzvals' and/or 'nzoffs' "
	      "are invalid or incompatible");

}

static inline SEXPTYPE get_SV_Rtype(const SparseVec *sv)
{
	return sv->Rtype;
}

static inline int get_SV_nzcount(const SparseVec *sv)
{
	return sv->nzcount;
}

static inline const Rbyte *get_RbyteSV_nzvals_p(const SparseVec *sv)
{
	return RAW(sv->nzvals);
}

static inline const int *get_intSV_nzvals_p(const SparseVec *sv)
{
	return INTEGER(sv->nzvals);
}

static inline const double *get_doubleSV_nzvals_p(const SparseVec *sv)
{
	return REAL(sv->nzvals);
}

static inline const Rcomplex *get_RcomplexSV_nzvals_p(const SparseVec *sv)
{
	return COMPLEX(sv->nzvals);
}

static inline Rbyte get_RbyteSV_nzval(const SparseVec *sv, int k)
{
	return sv->nzvals == R_NilValue ? Rbyte1 : RAW(sv->nzvals)[k];
}

static inline int get_intSV_nzval(const SparseVec *sv, int k)
{
	return sv->nzvals == R_NilValue ? int1 : INTEGER(sv->nzvals)[k];
}

static inline double get_doubleSV_nzval(const SparseVec *sv, int k)
{
	return sv->nzvals == R_NilValue ? double1 : REAL(sv->nzvals)[k];
}

static inline Rcomplex get_RcomplexSV_nzval(const SparseVec *sv, int k)
{
	return sv->nzvals == R_NilValue ? Rcomplex1 : COMPLEX(sv->nzvals)[k];
}

static inline int smallest_offset(
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

#define FUNDEF_next_2SV_vals(Ltype, Rtype)				\
	(const SparseVec *sv1,						\
	 const SparseVec *sv2,						\
	 int *k1, int *k2, int *off, Ltype *val1, Rtype *val2)		\
{									\
	int ret = smallest_offset(sv1->nzoffs, get_SV_nzcount(sv1),	\
				  sv2->nzoffs, get_SV_nzcount(sv2),	\
				  *k1, *k2, off);			\
	switch (ret) {							\
	    case 1: {							\
		*val1 = get_ ## Ltype  ## SV_nzval(sv1, *k1);		\
		*val2 = Rtype ## 0;					\
		(*k1)++;						\
		break;							\
	    }								\
	    case 2: {							\
		*val1 = Ltype ## 0;					\
		*val2 = get_ ## Rtype  ## SV_nzval(sv2, *k2);		\
		(*k2)++;						\
		break;							\
	    }								\
	    case 3: {							\
		*val1 = get_ ## Ltype  ## SV_nzval(sv1, *k1);		\
		*val2 = get_ ## Rtype  ## SV_nzval(sv2, *k2);		\
		(*k1)++;						\
		(*k2)++;						\
		break;							\
	    }								\
	}								\
	return ret;							\
}

static inline int next_2SV_vals_Rbyte_Rbyte
	FUNDEF_next_2SV_vals(Rbyte, Rbyte)

static inline int next_2SV_vals_Rbyte_int
	FUNDEF_next_2SV_vals(Rbyte, int)

static inline int next_2SV_vals_Rbyte_double
	FUNDEF_next_2SV_vals(Rbyte, double)

static inline int next_2SV_vals_Rbyte_Rcomplex
	FUNDEF_next_2SV_vals(Rbyte, Rcomplex)

static inline int next_2SV_vals_int_int
	FUNDEF_next_2SV_vals(int, int)

static inline int next_2SV_vals_int_double
	FUNDEF_next_2SV_vals(int, double)

static inline int next_2SV_vals_int_Rcomplex
	FUNDEF_next_2SV_vals(int, Rcomplex)

static inline int next_2SV_vals_double_int
	FUNDEF_next_2SV_vals(double, int)

static inline int next_2SV_vals_double_double
	FUNDEF_next_2SV_vals(double, double)

static inline int next_2SV_vals_double_Rcomplex
	FUNDEF_next_2SV_vals(double, Rcomplex)

static inline int next_2SV_vals_Rcomplex_Rcomplex
	FUNDEF_next_2SV_vals(Rcomplex, Rcomplex)

/* PROPAGATE_NZOFFS is a special value returned by _Arith_sv1_scalar(),
   _Compare_sv1_scalar(), and other functions that take a single input
   SparseVec to indicate that the result of the operation is a sparse vector
   with the same nzoffs as the input ones and with a single nzval shared by
   all the nzoffs.
   IMPORTANT: If this is the case then the function doesn't write anything
   to output buffer 'out_nzoffs' and writes the single shared nzval to
   'out_nzvals[0]'. */
#define	PROPAGATE_NZOFFS   -1  /* must be a **negative** int */

/* Another special value used in SparseVec_Arith.c. */
#define	NZCOUNT_IS_NOT_SET -2  /* must be < 0 and != PROPAGATE_NZOFFS */

#endif  /* _SPARSEVEC_H_ */

