/****************************************************************************
 *                     Basic manipulation of SVT leaves                     *
 ****************************************************************************/
#include "leaf_utils.h"

#include "S4Vectors_interface.h"

#include "Rvector_utils.h"
#include "coerceVector2.h"
#include "SparseArray_class.h"  /* for _coercion_can_introduce_zeros() */

#include <string.h>  /* for memcpy() */


/****************************************************************************
 * _alloc_and_unzip_leaf()
 * _make_leaf_from_bufs()
 */

SEXP _alloc_leaf(int lv_len, SEXPTYPE Rtype)
{
	SEXP lv_offs, lv_vals, ans;

	lv_offs = PROTECT(NEW_INTEGER(lv_len));
	lv_vals = PROTECT(allocVector(Rtype, lv_len));
	ans = zip_leaf(lv_offs, lv_vals);
	UNPROTECT(2);
	return ans;
}

/* Do NOT use when 'lv_len' is 0. Leaf vectors of length 0 are ILLEGAL! Always
   use a R_NilValue instead. See leaf_utils.h */
SEXP _alloc_and_unzip_leaf(int lv_len, SEXPTYPE Rtype,
		SEXP *lv_offs, SEXP *lv_vals)
{
	SEXP ans;

	ans = PROTECT(_alloc_leaf(lv_len, Rtype));
	unzip_leaf(ans, lv_offs, lv_vals);
	UNPROTECT(1);
	return ans;
}

/* Returns R_NilValue (if 'buf_len' is 0) or a "leaf vector".
   Does NOT work for 'Rtype' == STRSXP or VECSXP. */
SEXP _make_leaf_from_bufs(SEXPTYPE Rtype,
		const int *nzoffs_buf, const void *nzvals_buf, int buf_len)
{
	size_t Rtype_size;
	SEXP ans, ans_offs, ans_vals;

	if (buf_len == 0)
		return R_NilValue;
	Rtype_size = _get_Rtype_size(Rtype);
	if (Rtype_size == 0)
		error("SparseArray internal error in "
		      "_make_leaf_from_bufs():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));
	ans = PROTECT(_alloc_and_unzip_leaf(buf_len, Rtype,
					    &ans_offs, &ans_vals));
	memcpy(INTEGER(ans_offs), nzoffs_buf, sizeof(int) * buf_len);
	memcpy(DATAPTR(ans_vals), nzvals_buf, Rtype_size * buf_len);
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * _make_leaf_from_Rsubvec()
 */

/* 'offs_buf' must be of length 'subvec_len' (or more).
   Returns R_NilValue or a "leaf vector". */
SEXP _make_leaf_from_Rsubvec(
		SEXP Rvector, R_xlen_t vec_offset, int subvec_len,
		int *offs_buf, int avoid_copy_if_all_nonzeros)
{
	int ans_len;
	SEXP ans_offs, ans_vals, ans;

	ans_len = _collect_offsets_of_nonzero_Rsubvec_elts(
				Rvector, vec_offset, subvec_len,
				offs_buf);
	if (ans_len == 0)
		return R_NilValue;

	ans_offs = PROTECT(NEW_INTEGER(ans_len));
	memcpy(INTEGER(ans_offs), offs_buf, sizeof(int) * ans_len);

	if (avoid_copy_if_all_nonzeros && ans_len == subvec_len &&
	    vec_offset == 0 && XLENGTH(Rvector) == subvec_len)
	{
		/* The full 'Rvector' contains no zeros and can be reused
		   as-is without the need to copy its nonzero elements to
		   a new SEXP. */
		ans = zip_leaf(ans_offs, Rvector);
		UNPROTECT(1);
		return ans;
	}

	ans_vals = PROTECT(allocVector(TYPEOF(Rvector), ans_len));
	_copy_selected_Rsubvec_elts(Rvector, vec_offset, offs_buf, ans_vals);
	ans = zip_leaf(ans_offs, ans_vals);
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * _expand_leaf()
 */

/* Assumes that 'out_Rvector' is long enough, has the right type, and
   is already filled with zeros e.g. it was created with
   _new_Rvector0(TYPEOF(lv), d). */
int _expand_leaf(SEXP lv, SEXP out_Rvector, R_xlen_t out_offset)
{
	SEXP lv_offs, lv_vals;
	int ret;

	ret = unzip_leaf(lv, &lv_offs, &lv_vals);
	if (ret < 0)
		return -1;
	_copy_Rvector_elts_to_offsets(lv_vals, INTEGER(lv_offs),
				      out_Rvector, out_offset);
	return 0;
}


/****************************************************************************
 * _remove_zeros_from_leaf()
 *
 * Returns R_NilValue or a "leaf vector" of the same length or shorter than
 * the input "leaf vector".
 */

SEXP _remove_zeros_from_leaf(SEXP lv, int *offs_buf)
{
	SEXP lv_offs, lv_vals, ans_offs, ans_vals, ans;
	int lv_len, ans_len;

	lv_len = unzip_leaf(lv, &lv_offs, &lv_vals);
	ans_len = _collect_offsets_of_nonzero_Rsubvec_elts(
				lv_vals, 0, lv_len,
				offs_buf);
	if (ans_len == 0)       /* all values in 'lv' are zeros */
		return R_NilValue;
	if (ans_len == lv_len)  /* all values in 'lv' are nonzeros */
		return lv;  /* no-op */

	ans_offs = PROTECT(NEW_INTEGER(ans_len));
	_copy_selected_ints(INTEGER(lv_offs), offs_buf, ans_len,
			    INTEGER(ans_offs));
	ans_vals = PROTECT(allocVector(TYPEOF(lv_vals), ans_len));
	_copy_selected_Rsubvec_elts(lv_vals, 0, offs_buf, ans_vals);
	ans = zip_leaf(ans_offs, ans_vals);
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * _coerce_leaf()
 *
 * Returns R_NilValue or a "leaf vector" of the same length or shorter than
 * the input "leaf vector".
 */

/* Note that, in practice, _coerce_leaf() is always called to actually
   change the type of 'lv', so the code below does not bother to check
   for the (trivial) no-op case. */
SEXP _coerce_leaf(SEXP lv, SEXPTYPE new_Rtype, int *warn, int *offs_buf)
{
	SEXP lv_offs, lv_vals, ans_vals, ans;

	unzip_leaf(lv, &lv_offs, &lv_vals);
	ans_vals = PROTECT(_coerceVector2(lv_vals, new_Rtype, warn));
	ans = PROTECT(zip_leaf(lv_offs, ans_vals));
	/* The above coercion can introduce zeros in 'ans_vals' e.g. when
	   going from double/complex to int/raw. We need to remove them. */
	if (_coercion_can_introduce_zeros(TYPEOF(lv_vals), new_Rtype))
		ans = _remove_zeros_from_leaf(ans, offs_buf);
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * _subassign_leaf_with_Rvector()
 *
 * 'lv' must be a "leaf vector" (cannot be R_NilValue).
 * 'index' must be an integer vector containing valid zero-based indices
 * (a.k.a. offsets) into 'lv'. They are expected to be already sorted in
 * strictly ascending order.
 * 'Rvector' must be a vector (atomic or list) of the same length as 'index'.
 * Returns a "leaf vector" whose length is guaranteed to not exceed
 * min(length(lv) + length(index), INT_MAX). Note that the zeros in 'Rvector'
 * will be injected in the returned "leaf vector". This function does NOT
 * remove them!
 */

SEXP _subassign_leaf_with_Rvector(SEXP lv, SEXP index, SEXP Rvector)
{
	int lv_len, index_len, ans_len, k1, k2, k, n;
	SEXP lv_offs, lv_vals, ans, ans_offs, ans_vals;
	SEXPTYPE Rtype;
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;
	CopyRVectorElts_FUNType copy_Rvector_elts_FUN;
	const int *offs1_p, *offs2_p;
	int *ans_offs_p;

	lv_len = unzip_leaf(lv, &lv_offs, &lv_vals);
	Rtype = TYPEOF(lv_vals);

	copy_Rvector_elt_FUN = _select_copy_Rvector_elt_FUN(Rtype);
	copy_Rvector_elts_FUN = _select_copy_Rvector_elts_FUN(Rtype);
	if (copy_Rvector_elt_FUN == NULL || copy_Rvector_elts_FUN == NULL)
		error("SparseArray internal error in "
		      "_subassign_leaf_with_Rvector():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));

	if (TYPEOF(Rvector) != Rtype)
		error("SparseArray internal error in "
		      "_subassign_leaf_with_Rvector():\n"
		      "    'lv' and 'Rvector' have different types");

	index_len = LENGTH(index);
	if (LENGTH(Rvector) != index_len)
		error("SparseArray internal error in "
		      "_subassign_leaf_with_Rvector():\n"
		      "    'index' and 'Rvector' have different lengths");
	if (index_len == 0)
		return zip_leaf(lv_offs, lv_vals);

	offs1_p = INTEGER(lv_offs);
	offs2_p = INTEGER(index);
	ans_len = k1 = k2 = 0;
	while (k1 < lv_len && k2 < index_len) {
		if (*offs1_p < *offs2_p) {
			offs1_p++;
			k1++;
		} else if (*offs1_p > *offs2_p) {
			offs2_p++;
			k2++;
		} else {
			/* *offs1_p == *offs2_p */
			offs1_p++;
			k1++;
			offs2_p++;
			k2++;
		}
		ans_len++;
	}
	if (k1 < lv_len) {
		ans_len += lv_len - k1;
	} else if (k2 < index_len) {
		ans_len += index_len - k2;
	}
	ans = PROTECT(_alloc_and_unzip_leaf(ans_len, Rtype,
					    &ans_offs, &ans_vals));
	offs1_p = INTEGER(lv_offs);
	offs2_p = INTEGER(index);
	ans_offs_p = INTEGER(ans_offs);
	k = k1 = k2 = 0;
	while (k1 < lv_len && k2 < index_len) {
		if (*offs1_p < *offs2_p) {
			*ans_offs_p = *offs1_p;
			copy_Rvector_elt_FUN(lv_vals, (R_xlen_t) k1,
					     ans_vals, (R_xlen_t) k);
			offs1_p++;
			k1++;
		} else if (*offs1_p > *offs2_p) {
			*ans_offs_p = *offs2_p;
			copy_Rvector_elt_FUN(Rvector, (R_xlen_t) k2,
					     ans_vals, (R_xlen_t) k);
			offs2_p++;
			k2++;
		} else {
			/* *offs1_p == *offs2_p */
			*ans_offs_p = *offs2_p;
			copy_Rvector_elt_FUN(Rvector, (R_xlen_t) k2,
					     ans_vals, (R_xlen_t) k);
			offs1_p++;
			k1++;
			offs2_p++;
			k2++;
		}
		ans_offs_p++;
		k++;
	}
	if (k1 < lv_len) {
		n = lv_len - k1;
		memcpy(ans_offs_p, offs1_p, sizeof(int) * n);
		copy_Rvector_elts_FUN(lv_vals, (R_xlen_t) k1,
				      ans_vals, (R_xlen_t) k,
				      (R_xlen_t) n);
	} else if (k2 < index_len) {
		n = index_len - k2;
		memcpy(ans_offs_p, offs2_p, sizeof(int) * n);
		copy_Rvector_elts_FUN(Rvector, (R_xlen_t) k2,
				      ans_vals, (R_xlen_t) k,
				      (R_xlen_t) n);
	}
	UNPROTECT(1);
	return ans;
}

