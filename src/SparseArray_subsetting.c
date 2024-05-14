/****************************************************************************
 *                       Extract a SparseArray subset                       *
 ****************************************************************************/
#include "SparseArray_subsetting.h"

#include "Rvector_utils.h"
#include "leaf_utils.h"

#include <limits.h>  /* for INT_MAX */
#include <string.h>  /* for memcpy() */


static SEXP compute_subset_dim(SEXP index, SEXP x_dim)
{
	int ndim, along;
	SEXP ans_dim, index_elt;
	R_xlen_t d;

	ndim = LENGTH(x_dim);
	if (!isVectorList(index) || LENGTH(index) != ndim)
		error("'index' must be a list with one list "
		      "element along each dimension in 'x'");

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
		const int *nzoffs, int nzcount)
{
	for (int k = 0; k < nzcount; k++)
		lookup_table[*(nzoffs++)] = k;
	return;
}

static void reset_lookup_table(int *lookup_table,
		const int *nzoffs, int nzcount)
{
	for (int k = 0; k < nzcount; k++)
		lookup_table[*(nzoffs)++] = -1;
	return;
}

static inline int map_i2_to_k2_with_lookup_table(int i2,
		const int *lookup_table)
{
	return lookup_table[i2];
}

/* Returns a value >= 0 and < 'nzcount' if success, or -1 if failure. */
static inline int map_i2_to_k2_with_bsearch(int i2,
		const int *nzoffs, int nzcount)
{
	int k1, k2, k, off;

	/* Compare with first offset. */
	k1 = 0;
	off = nzoffs[k1];
	if (i2 < off)
		return -1;
	if (i2 == off)
		return k1;

	/* Compare with last offset. */
	k2 = nzcount - 1;
	off = nzoffs[k2];
	if (i2 > off)
		return -1;
	if (i2 == off)
		return k2;

	/* Binary search.
	   Seems that using >> 1 instead of / 2 is faster, even when compiling
	   with 'gcc -O2' (one would hope that the optimizer is able to do that
	   kind of optimization). */
	while ((k = (k1 + k2) >> 1) != k1) {
		off = nzoffs[k];
		if (i2 == off)
			return k;
		if (i2 > off)
			k1 = k;
		else
			k2 = k;
	}
	return -1;
}

/* Takes a non-NULL leaf (standard or lacunar), and returns a leaf that can
   be NULL, standard, or lacunar. */
static SEXP subset_leaf(SEXP leaf, SEXP idx, int i2max,
		int *i1_buf, int *k2_buf, int *lookup_table)
{
	if (idx == R_NilValue)
		return leaf;

	int idx_len = LENGTH(idx);
	if (idx_len == 0)
		return R_NilValue;

	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
	build_lookup_table(lookup_table, INTEGER(nzoffs), nzcount);
	int ans_nzcount = 0;
	for (int i1 = 0; i1 < idx_len; i1++) {
		int i2 = get_i2(INTEGER(idx), i1, i2max);
		//k2 = map_i2_to_k2_with_bsearch(i2, INTEGER(nzoffs), nzcount);
		//k2 = map_i2_to_k2_with_lookup_table(i2, lookup_table);
		int k2 = lookup_table[i2];
		if (k2 >= 0) {
			i1_buf[ans_nzcount] = i1;
			k2_buf[ans_nzcount] = k2;
			ans_nzcount++;
		}
	}
	reset_lookup_table(lookup_table, INTEGER(nzoffs), nzcount);
	if (ans_nzcount == 0)
		return R_NilValue;

	SEXP ans_nzoffs = PROTECT(NEW_INTEGER(ans_nzcount));
	memcpy(INTEGER(ans_nzoffs), i1_buf, sizeof(int) * ans_nzcount);
	if (nzvals == R_NilValue) {
		/* Leaf to subset is lacunar --> subsetting preserves that. */
		SEXP ans = _make_lacunar_leaf(ans_nzoffs);
		UNPROTECT(1);
		return ans;
	}
	if (LACUNAR_MODE_IS_ON) {
		int all_ones = _all_selected_Rsubvec_elts_equal_one(nzvals, 0,
						     k2_buf, ans_nzcount);
		if (all_ones) {
			SEXP ans = _make_lacunar_leaf(ans_nzoffs);
			UNPROTECT(1);
			return ans;
		}
	}
	SEXP ans_nzvals =
		PROTECT(_subset_Rsubvec(nzvals, 0, k2_buf, ans_nzcount));
	SEXP ans = zip_leaf(ans_nzvals, ans_nzoffs);
	UNPROTECT(2);
	return ans;
}

/* Recursive.
   Returns R_NilValue or a list of length 'ans_dim[ndim - 1]'. */
static SEXP REC_subset_SVT(SEXP SVT, SEXP index,
		const int *x_dim, const int *ans_dim, int ndim,
		int *i1_buf, int *k2_buf, int *lookup_table)
{
	if (SVT == R_NilValue)
		return R_NilValue;

	SEXP idx = VECTOR_ELT(index, ndim - 1);

	if (ndim == 1) {
		/* 'SVT' is a leaf. */
		return subset_leaf(SVT, idx, x_dim[ndim - 1],
				   i1_buf, k2_buf, lookup_table);
	}

	/* 'SVT' is a regular node (list). */
	int SVT_len = LENGTH(SVT);
	int ans_len = ans_dim[ndim - 1];
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	int is_empty = 1;
	for (int i1 = 0; i1 < ans_len; i1++) {
		int i2;
		if (idx == R_NilValue) {
			i2 = i1;
		} else {
			i2 = get_i2(INTEGER(idx), i1, SVT_len);
		}
		SEXP SVT_elt = VECTOR_ELT(SVT, i2);
		SEXP ans_elt = REC_subset_SVT(SVT_elt, index,
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

/* --- .Call ENTRY POINT ---
   'index' must be an N-index, that is, a list of integer vectors (or NULLs),
   one along each dimension in the array. */
SEXP C_subset_SVT_SparseArray(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		SEXP index)
{
	SEXPTYPE Rtype = _get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("SparseArray internal error in "
		      "C_subset_SVT_SparseArray():\n"
		      "    SVT_SparseArray object has invalid type");

	SEXP ans_dim = PROTECT(compute_subset_dim(index, x_dim));
	int x_dim0 = INTEGER(x_dim)[0];
	int *lookup_table = (int *) R_alloc(x_dim0, sizeof(int));
	int ans_dim0 = INTEGER(ans_dim)[0];
	int *i1_buf = (int *) R_alloc(ans_dim0, sizeof(int));
	int *k2_buf = (int *) R_alloc(ans_dim0, sizeof(int));
	for (int i = 0; i < x_dim0; i++)
		lookup_table[i] = -1;
	SEXP ans_SVT = REC_subset_SVT(x_SVT, index,
				 INTEGER(x_dim),
				 INTEGER(ans_dim), LENGTH(ans_dim),
				 i1_buf, k2_buf, lookup_table);
	if (ans_SVT != R_NilValue)
		PROTECT(ans_SVT);

	SEXP ans = PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(ans, 0, ans_dim);
	if (ans_SVT != R_NilValue) {
		SET_VECTOR_ELT(ans, 1, ans_SVT);
		UNPROTECT(1);
	}
	UNPROTECT(2);
	return ans;
}

