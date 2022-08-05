#ifndef _LEAF_VECTOR_UTILS_H_
#define _LEAF_VECTOR_UTILS_H_

#include <Rdefines.h>
#include "Rvector_utils.h"
#include "Rvector_summarize.h"

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

SEXP _new_leaf_vector(
	SEXP lv_offs,
	SEXP lv_vals
);

static inline int _split_leaf_vector(SEXP lv, SEXP *lv_offs, SEXP *lv_vals)
{
	R_xlen_t lv_offs_len;

	/* Sanity checks (should never fail). */
	if (!isVectorList(lv))  // IS_LIST() is broken
		return -1;
	if (LENGTH(lv) != 2)
		return -1;
	*lv_offs = VECTOR_ELT(lv, 0);
	*lv_vals = VECTOR_ELT(lv, 1);
	if (!IS_INTEGER(*lv_offs))
		return -1;
	lv_offs_len = XLENGTH(*lv_offs);
	if (lv_offs_len > INT_MAX)
		return -1;
	if (XLENGTH(*lv_vals) != lv_offs_len)
		return -1;
	return (int) lv_offs_len;
}

SEXP _alloc_leaf_vector(
	int lv_len,
	SEXPTYPE Rtype
);

SEXP _alloc_and_split_leaf_vector(
	int lv_len,
	SEXPTYPE Rtype,
	SEXP *lv_offs,
	SEXP *lv_vals
);

SEXP _make_leaf_vector_from_Rsubvec(
	SEXP Rvector,
	R_xlen_t subvec_offset,
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

int _summarize_leaf_vector(
	SEXP lv,
	int d,
	const SummarizeOp *summarize_op,
	void *init,
	R_xlen_t *na_rm_count,
	int status
);

#endif  /* _LEAF_VECTOR_UTILS_H_ */

