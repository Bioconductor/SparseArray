#ifndef _SPARSE_VEC_H_
#define _SPARSE_VEC_H_

#include <Rdefines.h>

struct sparse_vec {
	SEXPTYPE Rtype;      /* type of values */
	const void *nzvals;  /* nonzero values */
	const int *nzoffs;   /* offsets of nonzero values */
	int nzcount;         /* number of nonzero values */
	int len;             /* vector length (-1 if unknown) */
};

static inline struct sparse_vec new_sparse_vec(
		SEXPTYPE Rtype, const void *nzvals,
		const int *nzoffs, int nzcount, int len)
{
	struct sparse_vec sv;

	sv.Rtype = Rtype;
	sv.nzvals = nzvals;
	sv.nzoffs = nzoffs;
	sv.nzcount = nzcount;
	sv.len = len;
	return sv;
}

static inline struct sparse_vec new_Rbyte_sparse_vec(const Rbyte *nzvals,
		const int *nzoffs, int nzcount, int len)
{
	return new_sparse_vec(RAWSXP, nzvals, nzoffs, nzcount, len);
}

static inline struct sparse_vec new_int_sparse_vec(const int *nzvals,
		const int *nzoffs, int nzcount, int len)
{
	return new_sparse_vec(INTSXP, nzvals, nzoffs, nzcount, len);
}

static inline struct sparse_vec new_double_sparse_vec(const double *nzvals,
		const int *nzoffs, int nzcount, int len)
{
	return new_sparse_vec(REALSXP, nzvals, nzoffs, nzcount, len);
}

static inline struct sparse_vec new_Rcomplex_sparse_vec(const Rcomplex *nzvals,
		const int *nzoffs, int nzcount, int len)
{
	return new_sparse_vec(CPLXSXP, nzvals, nzoffs, nzcount, len);
}

static inline const Rbyte *get_Rbyte_nzvals(const struct sparse_vec *sv)
{
        if (sv->Rtype != RAWSXP)
		error("SparseArray internal error in "
		      "get_Rbyte_nzvals():\n"
		      "    sv->Rtype != RAWSXP");
	return sv->nzvals;
}

static inline const int *get_int_nzvals(const struct sparse_vec *sv)
{
        if (sv->Rtype != INTSXP && sv->Rtype != LGLSXP)
		error("SparseArray internal error in "
		      "get_int_nzvals():\n"
		      "    sv->Rtype != INTSXP && sv->Rtype != LGLSXP");
	return sv->nzvals;
}

static inline const double *get_double_nzvals(const struct sparse_vec *sv)
{
        if (sv->Rtype != REALSXP)
		error("SparseArray internal error in "
		      "get_double_nzvals():\n"
		      "    sv->Rtype != REALSXP");
	return sv->nzvals;
}

static inline const Rcomplex *get_Rcomplex_nzvals(const struct sparse_vec *sv)
{
        if (sv->Rtype != CPLXSXP)
		error("SparseArray internal error in "
		      "get_Rcomplex_nzvals():\n"
		      "    sv->Rtype != CPLXSXP");
	return sv->nzvals;
}

static const Rbyte Rbyte0 = 0;
static const int int0 = 0;
static const double double0 = 0.0;
static const Rcomplex Rcomplex0 = {{0.0, 0.0}};

#define FUNDEF_next_nzvals(Ltype, Rtype)			\
	(const struct sparse_vec *sv1,				\
	 const struct sparse_vec *sv2,				\
	 int *k1, int *k2, int *off, Ltype *val1, Rtype *val2)	\
{								\
	int nzcount1 = sv1->nzcount;				\
	int nzcount2 = sv2->nzcount;				\
	const int *nzoffs1 = sv1->nzoffs;			\
	const int *nzoffs2 = sv2->nzoffs;			\
	const Ltype *nzvals1 = get_ ## Ltype  ## _nzvals(sv1);	\
	const Rtype *nzvals2 = get_ ## Rtype  ## _nzvals(sv2);	\
	if (*k1 < nzcount1 && *k2 < nzcount2) {			\
		int nzoff1 = nzoffs1[*k1];			\
		int nzoff2 = nzoffs2[*k2];			\
		if (nzoff1 < nzoff2) {				\
			*off = nzoff1;				\
			*val1 = nzvals1[(*k1)++];		\
			*val2 = Rtype ## 0;			\
			return 1;				\
		}						\
		if (nzoff1 > nzoff2) {				\
			*off = nzoff2;				\
			*val1 = Ltype ## 0;			\
			*val2 = nzvals2[(*k2)++];		\
			return 2;				\
		}						\
		*off = nzoff1;					\
		*val1 = nzvals1[(*k1)++];			\
		*val2 = nzvals2[(*k2)++];			\
		return 3;					\
	}							\
	if (*k1 < nzcount1) {					\
		*off = nzoffs1[*k1];				\
		*val1 = nzvals1[(*k1)++];			\
		*val2 = Rtype ## 0;				\
		return 1;					\
	}							\
	if (*k2 < nzcount2) {					\
		*off = nzoffs2[*k2];				\
		*val1 = Ltype ## 0;				\
		*val2 = nzvals2[(*k2)++];			\
		return 2;					\
	}							\
	return 0;						\
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

static inline struct sparse_vec leaf2sv(SEXP leaf, int len)
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
	return new_sparse_vec(TYPEOF(nzvals), DATAPTR(nzvals),
			      INTEGER(nzoffs), (int) nzcount, len);

    on_error:
	error("SparseArray internal error in "
	      "leaf2sv():\n"
	      "    invalid tree leaf");
}

#endif  /* _SPARSE_VEC_H_ */

