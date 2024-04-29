#ifndef _LEAF_VECTOR_UTILS_H_
#define _LEAF_VECTOR_UTILS_H_

#include <Rdefines.h>

#include "sparse_vec.h"

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

#define FUNDEF_next_nzval(Ltype, Rtype)				\
	(const int *offs1, const Ltype *vals1, int n1,		\
	 const int *offs2, const Rtype *vals2, int n2,		\
	 int *k1, int *k2, int *off, Ltype *val1, Rtype *val2)	\
{								\
	struct sparse_vec sv1 =					\
		new_ ## Ltype ## _sparse_vec(vals1, offs1, n1);	\
	struct sparse_vec sv2 =					\
		new_ ## Rtype ## _sparse_vec(vals2, offs2, n2);	\
	return next_nzvals_ ## Ltype ## _ ## Rtype		\
			(&sv1, &sv2, k1, k2, off, val1, val2);	\
}

static inline int next_nzval_Rbyte_Rbyte
	FUNDEF_next_nzval(Rbyte, Rbyte)

static inline int next_nzval_Rbyte_int
	FUNDEF_next_nzval(Rbyte, int)

static inline int next_nzval_Rbyte_double
	FUNDEF_next_nzval(Rbyte, double)

static inline int next_nzval_Rbyte_Rcomplex
	FUNDEF_next_nzval(Rbyte, Rcomplex)

static inline int next_nzval_int_int
	FUNDEF_next_nzval(int, int)

static inline int next_nzval_int_double
	FUNDEF_next_nzval(int, double)

static inline int next_nzval_int_Rcomplex
	FUNDEF_next_nzval(int, Rcomplex)

static inline int next_nzval_double_int
	FUNDEF_next_nzval(double, int)

static inline int next_nzval_double_double
	FUNDEF_next_nzval(double, double)

static inline int next_nzval_double_Rcomplex
	FUNDEF_next_nzval(double, Rcomplex)

static inline int next_nzval_Rcomplex_Rcomplex
	FUNDEF_next_nzval(Rcomplex, Rcomplex)

SEXP _new_leaf_vector(
	SEXP lv_offs,
	SEXP lv_vals
);

static inline int _split_leaf_vector(SEXP lv, SEXP *lv_offs, SEXP *lv_vals)
{
	R_xlen_t lv_offs_len;

	/* Sanity checks (should never fail). */
	if (!isVectorList(lv))  // IS_LIST() is broken
		goto on_error;
	if (LENGTH(lv) != 2)
		goto on_error;
	*lv_offs = VECTOR_ELT(lv, 0);
	if (!IS_INTEGER(*lv_offs))
		goto on_error;
	lv_offs_len = XLENGTH(*lv_offs);
	if (lv_offs_len == 0 || lv_offs_len > INT_MAX)
		goto on_error;
	*lv_vals = VECTOR_ELT(lv, 1);
	if (XLENGTH(*lv_vals) != lv_offs_len)
		goto on_error;
	return (int) lv_offs_len;

    on_error:
	error("SparseArray internal error in "
              "_split_leaf_vector():\n"
	      "    invalid leaf vector");
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

SEXP _alloc_and_split_leaf_vector(
	int lv_len,
	SEXPTYPE Rtype,
	SEXP *lv_offs,
	SEXP *lv_vals
);

SEXP _make_leaf_vector_from_bufs(
	SEXPTYPE Rtype,
	const int *offs_buf,
	const void * vals_buf,
	int buf_len
);

SEXP _make_leaf_vector_from_Rsubvec(
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

SEXP _lv_apply_to_REALSXP(
	SEXP lv,
	apply_2double_FUNS *funs,
	int *offs_buf,
	double *vals_buf
);

#endif  /* _LEAF_VECTOR_UTILS_H_ */

