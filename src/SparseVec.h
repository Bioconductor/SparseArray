#ifndef _SPARSEVEC_H_
#define _SPARSEVEC_H_

#include <Rdefines.h>

typedef struct sparse_vec_t {
	SEXPTYPE Rtype;      /* type of values */
	const void *nzvals;  /* nonzero values */
	const int *nzoffs;   /* offsets of nonzero values */
	int nzcount;         /* number of nonzero values */
	int len;             /* vector length (-1 if unknown) */
} SparseVec;

static inline SparseVec make_SparseVec(
		SEXPTYPE Rtype, const void *nzvals,
		const int *nzoffs, int nzcount, int len)
{
	SparseVec sv;

	sv.Rtype = Rtype;
	sv.nzvals = nzvals;
	sv.nzoffs = nzoffs;
	sv.nzcount = nzcount;
	sv.len = len;
	return sv;
}

static inline SparseVec make_RbyteSV(const Rbyte *nzvals,
		const int *nzoffs, int nzcount, int len)
{
	return make_SparseVec(RAWSXP, nzvals, nzoffs, nzcount, len);
}

static inline SparseVec make_intSV(const int *nzvals,
		const int *nzoffs, int nzcount, int len)
{
	return make_SparseVec(INTSXP, nzvals, nzoffs, nzcount, len);
}

static inline SparseVec make_doubleSV(const double *nzvals,
		const int *nzoffs, int nzcount, int len)
{
	return make_SparseVec(REALSXP, nzvals, nzoffs, nzcount, len);
}

static inline SparseVec make_RcomplexSV(const Rcomplex *nzvals,
		const int *nzoffs, int nzcount, int len)
{
	return make_SparseVec(CPLXSXP, nzvals, nzoffs, nzcount, len);
}

static inline const Rbyte *get_RbyteSV_nzvals(const SparseVec *sv)
{
        if (sv->Rtype != RAWSXP)
		error("SparseArray internal error in "
		      "get_RbyteSV_nzvals():\n"
		      "    sv->Rtype != RAWSXP");
	return sv->nzvals;
}

static inline const int *get_intSV_nzvals(const SparseVec *sv)
{
        if (sv->Rtype != INTSXP && sv->Rtype != LGLSXP)
		error("SparseArray internal error in "
		      "get_intSV_nzvals():\n"
		      "    sv->Rtype != INTSXP && sv->Rtype != LGLSXP");
	return sv->nzvals;
}

static inline const double *get_doubleSV_nzvals(const SparseVec *sv)
{
        if (sv->Rtype != REALSXP)
		error("SparseArray internal error in "
		      "get_doubleSV_nzvals():\n"
		      "    sv->Rtype != REALSXP");
	return sv->nzvals;
}

static inline const Rcomplex *get_RcomplexSV_nzvals(const SparseVec *sv)
{
        if (sv->Rtype != CPLXSXP)
		error("SparseArray internal error in "
		      "get_RcomplexSV_nzvals():\n"
		      "    sv->Rtype != CPLXSXP");
	return sv->nzvals;
}

static const Rbyte Rbyte0 = 0;
static const int int0 = 0;
static const double double0 = 0.0;
static const Rcomplex Rcomplex0 = {{0.0, 0.0}};

#define FUNDEF_next_nzvals(Ltype, Rtype)				\
	(const SparseVec *sv1,						\
	 const SparseVec *sv2,						\
	 int *k1, int *k2, int *off, Ltype *val1, Rtype *val2)		\
{									\
	int nzcount1 = sv1->nzcount;					\
	int nzcount2 = sv2->nzcount;					\
	const int *nzoffs1 = sv1->nzoffs;				\
	const int *nzoffs2 = sv2->nzoffs;				\
	const Ltype *nzvals1 = get_ ## Ltype  ## SV_nzvals(sv1);	\
	const Rtype *nzvals2 = get_ ## Rtype  ## SV_nzvals(sv2);	\
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

static inline SparseVec leaf2SV(SEXP leaf, int len)
{
	SEXP nzoffs, nzvals;
	R_xlen_t nzcount;

	/* Sanity checks (should never fail). */
	if (!isVectorList(leaf))  // IS_LIST() is broken
		goto on_error;
	if (LENGTH(leaf) != 2)
		goto on_error;
	nzoffs = VECTOR_ELT(leaf, 0);
	if (!IS_INTEGER(nzoffs))
		goto on_error;
	nzcount = XLENGTH(nzoffs);
	if (nzcount == 0 || nzcount > INT_MAX)
		goto on_error;
	nzvals = VECTOR_ELT(leaf, 1);
	if (XLENGTH(nzvals) != nzcount)
		goto on_error;
	return make_SparseVec(TYPEOF(nzvals), DATAPTR(nzvals),
			      INTEGER(nzoffs), (int) nzcount, len);

    on_error:
	error("SparseArray internal error in leaf2SV():\n"
	      "    invalid tree leaf");
}

#endif  /* _SPARSEVEC_H_ */

