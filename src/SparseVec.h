#ifndef _SPARSEVEC_H_
#define _SPARSEVEC_H_

#include <Rdefines.h>

#include "Rvector_utils.h"

#include <limits.h>  /* for INT_MAX */

typedef struct sparse_vec_t {
	SEXP nzvals;         /* vector type (atomic or list) */
	const int *nzoffs;   /* offsets of nonzero values */
	int len;             /* vector length (-1 if unknown) */
} SparseVec;

static inline SparseVec make_SparseVec(
		const SEXP nzvals, const int *nzoffs, int len)
{
	SparseVec sv;

	R_xlen_t nzcount = XLENGTH(nzvals);
	if (nzcount > INT_MAX)
		error("SparseArray internal error in make_SparseVec():\n"
		      "    number of nonzero values must be <= INT_MAX");
	sv.nzvals = nzvals;
	sv.nzoffs = nzoffs;
	sv.len = len;
	return sv;
}

static inline SEXPTYPE get_SV_Rtype(const SparseVec *sv)
{
	return TYPEOF(sv->nzvals);
}

static inline int get_SV_nzcount(const SparseVec *sv)
{
	return LENGTH(sv->nzvals);
}

static inline const Rbyte *get_RbyteSV_nzvals(const SparseVec *sv)
{
	return RAW(sv->nzvals);
}

static inline const int *get_intSV_nzvals(const SparseVec *sv)
{
	return INTEGER(sv->nzvals);
}

static inline const double *get_doubleSV_nzvals(const SparseVec *sv)
{
	return REAL(sv->nzvals);
}

static inline const Rcomplex *get_RcomplexSV_nzvals(const SparseVec *sv)
{
	return COMPLEX(sv->nzvals);
}

#define FUNDEF_next_nzvals(Ltype, Rtype)				\
	(const SparseVec *sv1,						\
	 const SparseVec *sv2,						\
	 int *k1, int *k2, int *off, Ltype *val1, Rtype *val2)		\
{									\
	const Ltype *nzvals1 = get_ ## Ltype  ## SV_nzvals(sv1);	\
	const Rtype *nzvals2 = get_ ## Rtype  ## SV_nzvals(sv2);	\
	int nzcount1 = get_SV_nzcount(sv1);				\
	int nzcount2 = get_SV_nzcount(sv2);				\
	const int *nzoffs1 = sv1->nzoffs;				\
	const int *nzoffs2 = sv2->nzoffs;				\
	if (*k1 < nzcount1 && *k2 < nzcount2) {				\
		int nzoff1 = nzoffs1[*k1];				\
		int nzoff2 = nzoffs2[*k2];				\
		if (nzoff1 < nzoff2) {					\
			*off = nzoff1;					\
			*val1 = nzvals1[(*k1)++];			\
			*val2 = Rtype ## 0;				\
			return 1;					\
		}							\
		if (nzoff1 > nzoff2) {					\
			*off = nzoff2;					\
			*val1 = Ltype ## 0;				\
			*val2 = nzvals2[(*k2)++];			\
			return 2;					\
		}							\
		*off = nzoff1;						\
		*val1 = nzvals1[(*k1)++];				\
		*val2 = nzvals2[(*k2)++];				\
		return 3;						\
	}								\
	if (*k1 < nzcount1) {						\
		*off = nzoffs1[*k1];					\
		*val1 = nzvals1[(*k1)++];				\
		*val2 = Rtype ## 0;					\
		return 1;						\
	}								\
	if (*k2 < nzcount2) {						\
		*off = nzoffs2[*k2];					\
		*val1 = Ltype ## 0;					\
		*val2 = nzvals2[(*k2)++];				\
		return 2;						\
	}								\
	return 0;							\
}

static inline int next_nzvals_Rbyte_Rbyte
	FUNDEF_next_nzvals(Rbyte, Rbyte)

static inline int next_nzvals_Rbyte_int
	FUNDEF_next_nzvals(Rbyte, int)

static inline int next_nzvals_Rbyte_double
	FUNDEF_next_nzvals(Rbyte, double)

static inline int next_nzvals_Rbyte_Rcomplex
	FUNDEF_next_nzvals(Rbyte, Rcomplex)

static inline int next_nzvals_int_int
	FUNDEF_next_nzvals(int, int)

static inline int next_nzvals_int_double
	FUNDEF_next_nzvals(int, double)

static inline int next_nzvals_int_Rcomplex
	FUNDEF_next_nzvals(int, Rcomplex)

static inline int next_nzvals_double_int
	FUNDEF_next_nzvals(double, int)

static inline int next_nzvals_double_double
	FUNDEF_next_nzvals(double, double)

static inline int next_nzvals_double_Rcomplex
	FUNDEF_next_nzvals(double, Rcomplex)

static inline int next_nzvals_Rcomplex_Rcomplex
	FUNDEF_next_nzvals(Rcomplex, Rcomplex)

#endif  /* _SPARSEVEC_H_ */

