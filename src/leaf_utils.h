#ifndef _LEAF_UTILS_H_
#define _LEAF_UTILS_H_

#include <Rdefines.h>

#include "SparseVec.h"

#include <limits.h>  /* for INT_MAX */

/* SVT leaves
   ----------
   The leaves of a Sparse Vector Tree (SVT) represent sparse vectors along
   the first dimension (a.k.a. innermost or fastest moving dimension) of the
   sparse array. They contain a collection of offset/value pairs sorted by
   strictly ascending offset.
   A leaf is represented by an R_NilValue if it's empty, or by a list of 2
   parallel dense vectors:
     - nzoffs: an integer vector of offsets (i.e. 0-based positions);
     - nzvals: a vector (atomic or list) of nonzero values (zeroes are
               not allowed).
   The common length of 'nzoffs' and 'nzvals' is called the "nonzero count"
   (nzcount) and is guaranteed to be >= 1. Also we don't support "long leaves"
   so 'nzcount' must always be <= INT_MAX. */

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

SEXP _alloc_leaf(
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

int _expand_leaf(
	SEXP lv,
	SEXP out_Rvector,
	R_xlen_t out_offset
);

SEXP _remove_zeros_from_leaf(
	SEXP lv,
	int *offs_buf
);

SEXP _coerce_leaf(
	SEXP lv,
	SEXPTYPE new_Rtype,
	int *warn,
	int *offs_buf
);

SEXP _subassign_leaf_with_Rvector(
	SEXP lv,
	SEXP index,
	SEXP Rvector
);

#endif  /* _LEAF_UTILS_H_ */

