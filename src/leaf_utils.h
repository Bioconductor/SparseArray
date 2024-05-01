#ifndef _LEAF_UTILS_H_
#define _LEAF_UTILS_H_

#include <Rdefines.h>

#include "SparseVec.h"

#include <limits.h>  /* for INT_MAX */

/* A "leaf vector" is a vector of offset/value pairs sorted by strictly
   ascending offset. It is represented by a list of 2 parallel vectors:
   an integer vector of offsets (i.e. 0-based positions) and a vector
   (atomic or list) of values that are typically but not necessarily nonzeros.
   The length of a "leaf vector" is the number of offset/value pairs in it.

   IMPORTANT: Empty leaf vectors are ILLEGAL! A "leaf vector" should **never**
   be empty i.e. it must **always** contain at least one offset/value pair.
   That's because you should always use a R_NilValue if you need to represent
   an empty "leaf vector". Furthermore, the length of a "leaf vector" is
   **always** >= 1 and <= INT_MAX. */

static inline SEXP zip_leaf(SEXP nzoffs, SEXP nzvals)
{
	/* Sanity checks (should never fail). */
	if (!IS_INTEGER(nzoffs))
		goto on_error;
	R_xlen_t nzcount = XLENGTH(nzoffs);
	if (nzcount == 0 || nzcount > INT_MAX || XLENGTH(nzvals) != nzcount)
		goto on_error;
	SEXP leaf = PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(leaf, 0, nzoffs);
	SET_VECTOR_ELT(leaf, 1, nzvals);
	UNPROTECT(1);
	return leaf;

    on_error:
	error("SparseArray internal error in zip_leaf():\n"
	      "    invalid 'nzoffs' or 'nzvals'");
}

static inline int unzip_leaf(SEXP leaf, SEXP *nzoffs, SEXP *nzvals)
{
	/* Sanity checks (should never fail). */
	if (!isVectorList(leaf))  // IS_LIST() is broken
		goto on_error;
	if (LENGTH(leaf) != 2)
		goto on_error;
	*nzoffs = VECTOR_ELT(leaf, 0);
	if (!IS_INTEGER(*nzoffs))
		goto on_error;
	R_xlen_t nzcount = XLENGTH(*nzoffs);
	if (nzcount == 0 || nzcount > INT_MAX)
		goto on_error;
	*nzvals = VECTOR_ELT(leaf, 1);
	if (XLENGTH(*nzvals) != nzcount)
		goto on_error;
	return (int) nzcount;

    on_error:
	error("SparseArray internal error in unzip_leaf():\n"
	      "    invalid SVT leaf");
}

static inline SparseVec leaf2SV(SEXP leaf, int len)
{
	SEXP nzoffs, nzvals;

	unzip_leaf(leaf, &nzoffs, &nzvals);
	return make_SparseVec(nzvals, INTEGER(nzoffs), len);
}

typedef struct apply_2double_funs_t {
	double (*Rbyte2double_FUN)(Rbyte);
	double (*int2double_FUN)(int);
	double (*double2double_FUN)(double);
	double (*Rcomplex2double_FUN)(Rcomplex);
} apply_2double_FUNS;

SEXP _alloc_leaf_vector(
	int lv_len,
	SEXPTYPE Rtype
);

SEXP _alloc_and_unzip_leaf(
	int lv_len,
	SEXPTYPE Rtype,
	SEXP *lv_offs,
	SEXP *lv_vals
);

SEXP _make_leaf_from_bufs(
	SEXPTYPE Rtype,
	const int *nzoffs_buf,
	const void *nzvals_buf,
	int buf_len
);

SEXP _make_leaf_from_Rsubvec(
	SEXP Rvector,
	R_xlen_t vec_offset,
	int subvec_len,
	int *offs_buf,
	int avoid_copy_if_all_nonzeros
);

int _expand_leaf_vector(
	SEXP lv,
	SEXP out_Rvector,
	R_xlen_t out_offset
);

SEXP _remove_zeros_from_leaf_vector(
	SEXP lv,
	int *offs_buf
);

SEXP _coerce_leaf_vector(
	SEXP lv,
	SEXPTYPE new_Rtype,
	int *warn,
	int *offs_buf
);

SEXP _subassign_leaf_vector_with_Rvector(
	SEXP lv,
	SEXP index,
	SEXP Rvector
);

SEXP _leaf_apply_to_REALSXP(
	SEXP leaf,
	apply_2double_FUNS *funs,
	int *nzoffs_buf,
	double *nzvals_buf
);

#endif  /* _LEAF_UTILS_H_ */

