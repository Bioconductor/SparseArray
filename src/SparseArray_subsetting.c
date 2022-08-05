/****************************************************************************
 *                       Extract a SparseArray subset                       *
 ****************************************************************************/
#include "SparseArray_subsetting.h"

#include "Rvector_utils.h"
#include "leaf_vector_utils.h"  /* for _split_leaf_vector() */

#include <limits.h>  /* for INT_MAX */
#include <string.h>  /* for memcpy() */


/****************************************************************************
 * C_drop_SVT_SparseArray_ineffective_dims()
 */

static int count_effective_dims(SEXP x_dim)
{
	int count, x_ndim, along, d;

	count = 0;
	x_ndim = LENGTH(x_dim);
	for (along = 0; along < x_ndim; along++) {
		d = INTEGER(x_dim)[along];
		if (d != 1)
			count++;
	}
	return count;
}

static SEXP drop_ineffective_dims(SEXP x_dim)
{
	int ans_ndim, x_ndim, along1, along2, d;
	SEXP ans_dim;

	ans_ndim = count_effective_dims(x_dim);
	if (ans_ndim == 0)
		return ScalarInteger(1);

	ans_dim = PROTECT(NEW_INTEGER(ans_ndim));
	along1 = 0;
	x_ndim = LENGTH(x_dim);
	for (along2 = 0; along2 < x_ndim; along2++) {
		d = INTEGER(x_dim)[along2];
		if (d != 1) {
			INTEGER(ans_dim)[along1] = d;
			along1++;
		}
	}
	UNPROTECT(1);
	return ans_dim;
}

static SEXP drop_dimnames_along_ineffective_dims(SEXP x_dimnames, SEXP x_dim)
{
	int ans_ndim, x_ndim, along1, along2, d;
	SEXP ans_dimnames;

	ans_ndim = count_effective_dims(x_dim);
	if (ans_ndim == 0)
		return NEW_LIST(1);

	ans_dimnames = PROTECT(NEW_LIST(ans_ndim));
	along1 = 0;
	x_ndim = LENGTH(x_dim);
	for (along2 = 0; along2 < x_ndim; along2++) {
		d = INTEGER(x_dim)[along2];
		if (d != 1) {
			SET_VECTOR_ELT(ans_dimnames, along1,
				       VECTOR_ELT(x_dimnames, along2));
			along1++;
		}
	}
	UNPROTECT(1);
	return ans_dimnames;
}

/* 'SVT' must be a depth-1 tree representing a 1xN matrix.
   Turn it into a "leaf vector". */
static SEXP make_leaf_vector_from_depth_one_SVT(SEXP SVT,
		SEXPTYPE Rtype, CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	int SVT_len, ans_len, i, lv_len;
	SEXP subSVT, ans_offs, ans_vals, lv_vals, ans;

	SVT_len = LENGTH(SVT);
	ans_len = 0;
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue)
			continue;
		/* 'subSVT' is a "leaf vector" of length 1. */
		lv_len = LENGTH(VECTOR_ELT(subSVT, 0));
		if (lv_len != 1)
			error("SparseArray internal error in "
			      "make_leaf_vector_from_depth_one_SVT():\n"
			      "    lv_len != 1");
		ans_len++;
	}
	if (ans_len == 0)
		error("SparseArray internal error in "
		      "make_leaf_vector_from_depth_one_SVT():\n"
		      "    ans_len == 0");
	ans_offs = PROTECT(NEW_INTEGER(ans_len));
	ans_vals = PROTECT(allocVector(Rtype, ans_len));
	ans_len = 0;
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue)
			continue;
		/* 'subSVT' is a "leaf vector" of length 1. */
		INTEGER(ans_offs)[ans_len] = i;
		lv_vals = VECTOR_ELT(subSVT, 1);
		copy_Rvector_elt_FUN(lv_vals, 0, ans_vals, ans_len);
		ans_len++;
	}
	ans = _new_leaf_vector(ans_offs, ans_vals);
	UNPROTECT(2);
	return ans;
}

/* Recursive. */
static SEXP REC_remove_SVT_unary_nodes(SEXP SVT, const int *dim, int ndim,
		const int *all_effective_dim,
		const int *any_effective_dim,
		SEXPTYPE Rtype, CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	int SVT_len, i;
	SEXP subSVT, ans, ans_elt;

	if (SVT == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		return SVT;
	}

	/* 'SVT' is a regular node (list). */

	/* Sanity check (should never fail). */
	if (!isVectorList(SVT))  // IS_LIST() is broken
		 error("SparseArray internal error in "
		       "REC_remove_SVT_unary_nodes():\n"
		       "    SVT node is not a list");

	SVT_len = LENGTH(SVT);

	/* Sanity check (should never fail). */
	if (SVT_len != dim[ndim - 1])
		 error("SparseArray internal error in "
		       "REC_remove_SVT_unary_nodes():\n"
		       "    SVT_len != dim[ndim - 1]");

	if (SVT_len == 1)
		return REC_remove_SVT_unary_nodes(
				VECTOR_ELT(SVT, 0), dim, ndim - 1,
				all_effective_dim, any_effective_dim,
				Rtype, copy_Rvector_elt_FUN);

	if (all_effective_dim[ndim - 2])
		return SVT;  /* nothing else to do */

	ans = PROTECT(NEW_LIST(SVT_len));
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue)
			continue;
		ans_elt = PROTECT(
			REC_remove_SVT_unary_nodes(
				subSVT, dim, ndim - 1,
				all_effective_dim, any_effective_dim,
				Rtype, copy_Rvector_elt_FUN)
		);
		SET_VECTOR_ELT(ans, i, ans_elt);
		UNPROTECT(1);
	}

	if (any_effective_dim[ndim - 2]) {
		UNPROTECT(1);
		return ans;  /* nothing else to do */
	}

	/* 'ans' represents a 1xN matrix. Turn it into a "leaf vector". */
	ans = make_leaf_vector_from_depth_one_SVT(ans, Rtype,
						  copy_Rvector_elt_FUN);
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_drop_SVT_SparseArray_ineffective_dims(SEXP x_dim, SEXP x_dimnames,
		SEXP x_type, SEXP x_SVT)
{
	SEXPTYPE Rtype;
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;
	int x_ndim, ok1, ok2, along, d;
	int *all_effective_dim, *any_effective_dim;
	SEXP ans_dim, ans_dimnames, ans_SVT, ans;

	Rtype = _get_Rtype_from_Rstring(x_type);
	copy_Rvector_elt_FUN = _select_copy_Rvector_elt_FUN(Rtype);
	if (copy_Rvector_elt_FUN == NULL)
		error("SparseArray internal error in "
		      "C_drop_SVT_SparseArray_ineffective_dims():\n"
		      "    SVT_SparseArray object has invalid type");

	x_ndim = LENGTH(x_dim);
	all_effective_dim = (int *) R_alloc(x_ndim, sizeof(int));
	any_effective_dim = (int *) R_alloc(x_ndim, sizeof(int));
	ok1 = 1;
	ok2 = 0;
	for (along = 0; along < x_ndim; along++) {
		d = INTEGER(x_dim)[along];
		if (d == 1)
			ok1 = 0;
		else
			ok2 = 1;
		all_effective_dim[along] = ok1;
		any_effective_dim[along] = ok2;
	}

	ans_dim = PROTECT(drop_ineffective_dims(x_dim));
	ans_dimnames = PROTECT(
		drop_dimnames_along_ineffective_dims(x_dimnames, x_dim)
	);

	if (x_SVT != R_NilValue)
		ans_SVT = PROTECT(
			REC_remove_SVT_unary_nodes(x_SVT,
					INTEGER(x_dim), x_ndim,
					all_effective_dim,
					any_effective_dim,
					Rtype, copy_Rvector_elt_FUN)
		);

	ans = PROTECT(NEW_LIST(3));
	SET_VECTOR_ELT(ans, 0, ans_dim);
	SET_VECTOR_ELT(ans, 1, ans_dimnames);
	if (x_SVT != R_NilValue) {
		SET_VECTOR_ELT(ans, 2, ans_SVT);
		UNPROTECT(1);
	}
	UNPROTECT(3);
	return ans;
}


/****************************************************************************
 * C_subset_SVT_SparseArray()
 */

static SEXP compute_subset_dim(SEXP index, SEXP x_dim)
{
	int ndim, along;
	SEXP ans_dim, index_elt;
	R_xlen_t d;

	ndim = LENGTH(x_dim);
	if (!isVectorList(index) || LENGTH(index) != ndim)
		error("'index' must be a list with one "
		      "list element per dimension in 'x'");

	ans_dim = PROTECT(duplicate(x_dim));

	for (along = 0; along < ndim; along++) {
		index_elt = VECTOR_ELT(index, along);
		if (index_elt == R_NilValue)
			continue;
		if (!IS_INTEGER(index_elt)) {
			UNPROTECT(1);
			error("each list element in 'index' must "
			      "be either NULL or an integer vector");
		}
		d = XLENGTH(index_elt);
		if (d > INT_MAX) {
			UNPROTECT(1);
			error("cannot select more than INT_MAX array "
			      "slice along any of the dimension");
		}
		INTEGER(ans_dim)[along] = (int) d;
	}

	UNPROTECT(1);
	return ans_dim;
}

static inline int get_i2(const int *idx, int i1, int i2max)
{
	int i2;

	i2 = idx[i1];
	if (i2 == NA_INTEGER) {
		UNPROTECT(1);
		error("'index' cannot contain NAs");
	}
	if (i2 < 1 || i2 > i2max) {
		UNPROTECT(1);
		error("'index' contains out-of-bound "
		      "indices");
	}
	return --i2;
}

static void build_lookup_table(int *lookup_table,
		const int *lv_offs, int lv_len)
{
	for (int k = 0; k < lv_len; k++)
		lookup_table[*(lv_offs++)] = k;
	return;
}

static void reset_lookup_table(int *lookup_table,
		const int *lv_offs, int lv_len)
{
	for (int k = 0; k < lv_len; k++)
		lookup_table[*(lv_offs)++] = -1;
	return;
}

static inline int map_i2_to_k2_with_lookup_table(int i2,
		const int *lookup_table)
{
	return lookup_table[i2];
}

/* Returns a value >= 0 and < 'lv_len' if success, or -1 if failure. */
static inline int map_i2_to_k2_with_bsearch(int i2,
		const int *lv_offs, int lv_len)
{
	int k1, k2, k, off;

	/* Compare with first offset. */
	k1 = 0;
	off = lv_offs[k1];
	if (i2 < off)
		return -1;
	if (i2 == off)
		return k1;

	/* Compare with last offset. */
	k2 = lv_len - 1;
	off = lv_offs[k2];
	if (i2 > off)
		return -1;
	if (i2 == off)
		return k2;

	/* Binary search.
	   Seems that using >> 1 instead of / 2 is faster, even when compiling
	   with 'gcc -O2' (one would hope that the optimizer is able to do that
	   kind of optimization). */
	while ((k = (k1 + k2) >> 1) != k1) {
		off = lv_offs[k];
		if (i2 == off)
			return k;
		if (i2 > off)
			k1 = k;
		else
			k2 = k;
	}
	return -1;
}

static SEXP subset_leaf_vector(SEXP lv, SEXP idx, int i2max,
		int *i1_buf, int *k2_buf, int *lookup_table)
{
	int idx_len, lv_len, ans_len, i1, i2, k2;
	SEXP lv_offs, lv_vals, ans_offs, ans_vals, ans;

	if (idx == R_NilValue)
		return lv;

	idx_len = LENGTH(idx);
	if (idx_len == 0)
		return R_NilValue;

	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	build_lookup_table(lookup_table, INTEGER(lv_offs), lv_len);
	ans_len = 0;
	for (i1 = 0; i1 < idx_len; i1++) {
		i2 = get_i2(INTEGER(idx), i1, i2max);
		//k2 = map_i2_to_k2_with_bsearch(i2, INTEGER(lv_offs), lv_len);
		//k2 = map_i2_to_k2_with_lookup_table(i2, lookup_table);
		k2 = lookup_table[i2];
		if (k2 >= 0) {
			i1_buf[ans_len] = i1;
			k2_buf[ans_len] = k2;
			ans_len++;
		}
	}
	reset_lookup_table(lookup_table, INTEGER(lv_offs), lv_len);
	if (ans_len == 0)
		return R_NilValue;

	ans_offs = PROTECT(NEW_INTEGER(ans_len));
	ans_vals = PROTECT(allocVector(TYPEOF(lv_vals), ans_len));
	memcpy(INTEGER(ans_offs), i1_buf, sizeof(int) * ans_len);
	_copy_selected_Rsubvec_elts(lv_vals, 0, k2_buf, ans_vals);
	ans = _new_leaf_vector(ans_offs, ans_vals);
	UNPROTECT(2);
	return ans;
}

/* Recursive.
   Returns R_NilValue or a list of length 'ans_dim[ndim - 1]'. */
static SEXP REC_subset_SVT(SEXP SVT, SEXP index,
		const int *x_dim, const int *ans_dim, int ndim,
		int *i1_buf, int *k2_buf, int *lookup_table)
{
	SEXP idx, ans, SVT_elt, ans_elt;
	int SVT_len, ans_len, is_empty, i1, i2;

	if (SVT == R_NilValue)
		return R_NilValue;

	idx = VECTOR_ELT(index, ndim - 1);

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		return subset_leaf_vector(SVT, idx, x_dim[ndim - 1],
					  i1_buf, k2_buf, lookup_table);
	}

	/* 'SVT' is a regular node (list). */
	SVT_len = LENGTH(SVT);
	ans_len = ans_dim[ndim - 1];
	ans = PROTECT(NEW_LIST(ans_len));
	is_empty = 1;
	for (i1 = 0; i1 < ans_len; i1++) {
		if (idx == R_NilValue) {
			i2 = i1;
		} else {
			i2 = get_i2(INTEGER(idx), i1, SVT_len);
		}
		SVT_elt = VECTOR_ELT(SVT, i2);
		ans_elt = REC_subset_SVT(SVT_elt, index,
					 x_dim, ans_dim, ndim - 1,
					 i1_buf, k2_buf, lookup_table);
		if (ans_elt != R_NilValue) {
			PROTECT(ans_elt);
			SET_VECTOR_ELT(ans, i1, ans_elt);
			UNPROTECT(1);
			is_empty = 0;
		}
	}
	UNPROTECT(1);
	return is_empty ? R_NilValue : ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_subset_SVT_SparseArray(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		SEXP index)
{
	SEXPTYPE Rtype;
	SEXP ans_dim, ans_SVT, ans;
	int *i1_buf, *k2_buf, *lookup_table, i;

	Rtype = _get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("SparseArray internal error in "
		      "C_subset_SVT_SparseArray():\n"
		      "    SVT_SparseArray object has invalid type");

	ans_dim = PROTECT(compute_subset_dim(index, x_dim));
	i1_buf = (int *) R_alloc(INTEGER(ans_dim)[0], sizeof(int));
	k2_buf = (int *) R_alloc(INTEGER(ans_dim)[0], sizeof(int));
	lookup_table = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	for (i = 0; i < INTEGER(x_dim)[0]; i++)
		lookup_table[i] = -1;
	ans_SVT = REC_subset_SVT(x_SVT, index,
				 INTEGER(x_dim),
				 INTEGER(ans_dim), LENGTH(ans_dim),
				 i1_buf, k2_buf, lookup_table);
	if (ans_SVT != R_NilValue)
		PROTECT(ans_SVT);

	ans = PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(ans, 0, ans_dim);
	if (ans_SVT != R_NilValue) {
		SET_VECTOR_ELT(ans, 1, ans_SVT);
		UNPROTECT(1);
	}
	UNPROTECT(2);
	return ans;
}

