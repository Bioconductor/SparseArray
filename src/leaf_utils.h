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
     - nzvals: a vector (atomic or list) of nonzero values (zeros are
               not allowed).
     - nzoffs: an integer vector of offsets (i.e. 0-based positions);
   The common length of 'nzvals' and 'nzoffs' is called the "nonzero count"
   (nzcount) and is guaranteed to be >= 1. Also we don't support "long leaves"
   so 'nzcount' must always be <= INT_MAX.
   Note that a leaf represents a 1D SVT. */

/* Is it ok to produce lacunar leaves? Proper handling of lacunar leaves
   is still work-in-progress so we can turn this off any time if things
   go wrong. */
#define LACUNAR_MODE_IS_ON 0

/* In-place replacement. Supplied 'nzvals' is trusted! */
static inline void replace_leaf_nzvals(SEXP leaf, SEXP nzvals)
{
	SET_VECTOR_ELT(leaf, 1, nzvals);
}

/* In-place replacement. Supplied 'nzoffs' is trusted! */
static inline void replace_leaf_nzoffs(SEXP leaf, SEXP nzoffs)
{
	SET_VECTOR_ELT(leaf, 0, nzoffs);
}

static inline SEXP zip_leaf(SEXP nzvals, SEXP nzoffs)
{
	/* Sanity checks (should never fail). */
	if (!IS_INTEGER(nzoffs))
		goto on_error;
	R_xlen_t nzcount = XLENGTH(nzoffs);
	if (nzcount == 0 || nzcount > INT_MAX)
		goto on_error;
	if (nzvals != R_NilValue && XLENGTH(nzvals) != nzcount)
		goto on_error;

	SEXP leaf = PROTECT(NEW_LIST(2));
	replace_leaf_nzvals(leaf, nzvals);
	replace_leaf_nzoffs(leaf, nzoffs);
	UNPROTECT(1);
	return leaf;

    on_error:
	error("SparseArray internal error in zip_leaf():\n"
	      "    supplied 'nzvals' and/or 'nzoffs' "
	      "are invalid or incompatible");
}

static inline SEXP get_leaf_nzvals(SEXP leaf)
{
	if (!isVectorList(leaf))  // IS_LIST() is broken
		goto on_error;
	/* A regular leaf is a list of length 2 but we don't test for
	   LENGTH(leaf) == 2 because we want this to work on an "extended
	   leaf" which is represented by a list of length 3. See
	   SparseArray_subassignment.c where "extended leaves" are explained
	   and used. */
	if (LENGTH(leaf) < 2)
		goto on_error;
	return VECTOR_ELT(leaf, 1);

    on_error:
	error("SparseArray internal error in get_leaf_nzvals():\n"
	      "    invalid SVT leaf");
}

static inline SEXP get_leaf_nzoffs(SEXP leaf)
{
	if (!isVectorList(leaf))  // IS_LIST() is broken
		goto on_error;
	if (LENGTH(leaf) < 2)  /* see why we don't do LENGTH(leaf) == 2 above */
		goto on_error;
	SEXP nzoffs = VECTOR_ELT(leaf, 0);
	if (!IS_INTEGER(nzoffs))
		goto on_error;
	R_xlen_t nzcount = XLENGTH(nzoffs);
	if (nzcount == 0 || nzcount > INT_MAX)
		goto on_error;
	return nzoffs;

    on_error:
	error("SparseArray internal error in get_leaf_nzoffs():\n"
	      "    invalid SVT leaf");
}

static inline int get_leaf_nzcount(SEXP leaf)
{
	return LENGTH(get_leaf_nzoffs(leaf));
}

static inline int unzip_leaf(SEXP leaf, SEXP *nzvals, SEXP *nzoffs)
{
	*nzvals = get_leaf_nzvals(leaf);
	*nzoffs = get_leaf_nzoffs(leaf);
	R_xlen_t nzcount = XLENGTH(*nzoffs);
	if (*nzvals != R_NilValue && XLENGTH(*nzvals) != nzcount)
		error("SparseArray internal error in unzip_leaf():\n"
		      "    invalid SVT leaf ('nzvals' and 'nzoffs' "
		      "are not parallel)");
	return (int) nzcount;
}

static inline SparseVec leaf2SV(SEXP leaf, SEXPTYPE Rtype, int len)
{
	SEXP nzvals, nzoffs;

	unzip_leaf(leaf, &nzvals, &nzoffs);
	return toSparseVec(nzvals, nzoffs, Rtype, len);
}

SEXP C_lacunar_mode_is_on(void);

SEXP _alloc_leaf(
	SEXPTYPE Rtype,
	int nzcount
);

SEXP _alloc_and_unzip_leaf(
	SEXPTYPE Rtype,
	int nzcount,
	SEXP *nzvals,
	SEXP *nzoffs
);

void _expand_leaf(
	SEXP leaf,
	SEXP out_Rvector,
	R_xlen_t out_offset
);

SEXP _make_leaf_from_two_arrays(
	SEXPTYPE Rtype,
	const void *nzvals_p,
	const int *nzoffs_p,
	int nzcount
);

SEXP _make_leaf_from_Rsubvec(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	int subvec_len,
	int *selection_buf,
	int avoid_copy_if_all_nonzeros
);

SEXP _INPLACE_remove_zeros_from_leaf(
	SEXP leaf,
	int *selection_buf
);

void _INPLACE_turn_into_lacunar_leaf_if_all_ones(SEXP leaf);

SEXP _coerce_leaf(
	SEXP leaf,
	SEXPTYPE new_Rtype,
	int *warn,
	int *selection_buf
);

SEXP _subassign_leaf_with_Rvector(
	SEXP leaf,
	SEXP index,
	SEXP Rvector
);

#endif  /* _LEAF_UTILS_H_ */

