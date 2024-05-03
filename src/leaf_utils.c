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
 * _alloc_leaf()
 * _alloc_and_unzip_leaf()
 * _make_leaf_from_bufs()
 */

/* Do NOT use when 'nzcount' is 0!
   Always use an R_NilValue to represent an empty leaf. See leaf_utils.h */
SEXP _alloc_leaf(SEXPTYPE Rtype, int nzcount)
{
	if (nzcount == 0)
		error("SparseArray internal error in _alloc_leaf():\n"
		      "    nzcount == 0");
	SEXP nzoffs = PROTECT(NEW_INTEGER(nzcount));
	SEXP nzvals = PROTECT(allocVector(Rtype, nzcount));
	SEXP ans = zip_leaf(nzoffs, nzvals);
	UNPROTECT(2);
	return ans;
}

SEXP _alloc_and_unzip_leaf(SEXPTYPE Rtype, int nzcount,
		SEXP *nzoffs, SEXP *nzvals)
{
	SEXP ans = PROTECT(_alloc_leaf(Rtype, nzcount));
	unzip_leaf(ans, nzoffs, nzvals);
	UNPROTECT(1);
	return ans;
}

/* Does NOT work for 'Rtype' == STRSXP or VECSXP. */
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
	ans = PROTECT(_alloc_and_unzip_leaf(Rtype, buf_len,
					    &ans_offs, &ans_vals));
	memcpy(INTEGER(ans_offs), nzoffs_buf, sizeof(int) * buf_len);
	memcpy(DATAPTR(ans_vals), nzvals_buf, Rtype_size * buf_len);
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * _make_leaf_from_Rsubvec()
 */

/* 'nzoffs_buf' must be of length 'subvec_len' (or more). */
SEXP _make_leaf_from_Rsubvec(
		SEXP Rvector, R_xlen_t vec_offset, int subvec_len,
		int *nzoffs_buf, int avoid_copy_if_all_nonzeros)
{
	int ans_nzcount;
	SEXP ans_nzoffs, ans_nzvals, ans;

	ans_nzcount = _collect_offsets_of_nonzero_Rsubvec_elts(
				Rvector, vec_offset, subvec_len,
				nzoffs_buf);
	if (ans_nzcount == 0)
		return R_NilValue;

	ans_nzoffs = PROTECT(NEW_INTEGER(ans_nzcount));
	memcpy(INTEGER(ans_nzoffs), nzoffs_buf, sizeof(int) * ans_nzcount);

	if (avoid_copy_if_all_nonzeros && ans_nzcount == subvec_len &&
	    vec_offset == 0 && XLENGTH(Rvector) == subvec_len)
	{
		/* The full 'Rvector' contains no zeros and can be reused
		   as-is without the need to copy its nonzero elements to
		   a new SEXP. */
		ans = zip_leaf(ans_nzoffs, Rvector);
		UNPROTECT(1);
		return ans;
	}

	ans_nzvals = PROTECT(allocVector(TYPEOF(Rvector), ans_nzcount));
	_copy_selected_Rsubvec_elts(Rvector, vec_offset,
				    nzoffs_buf, ans_nzvals);
	ans = zip_leaf(ans_nzoffs, ans_nzvals);
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * _expand_leaf()
 */

/* Assumes that 'out_Rvector' is long enough, has the right type, and
   is already filled with zeros e.g. it was created with
   _new_Rvector0(TYPEOF(leaf), d). */
int _expand_leaf(SEXP leaf, SEXP out_Rvector, R_xlen_t out_offset)
{
	SEXP nzoffs, nzvals;

	unzip_leaf(leaf, &nzoffs, &nzvals);  /* ignore returned nzcount */
	_copy_Rvector_elts_to_offsets(nzvals, INTEGER(nzoffs),
				      out_Rvector, out_offset);
	return 0;
}


/****************************************************************************
 * _remove_zeros_from_leaf()
 */

SEXP _remove_zeros_from_leaf(SEXP leaf, int *nzoffs_buf)
{
	SEXP nzoffs, nzvals, ans_nzoffs, ans_nzvals, ans;
	int nzcount, ans_nzcount;

	nzcount = unzip_leaf(leaf, &nzoffs, &nzvals);
	ans_nzcount = _collect_offsets_of_nonzero_Rsubvec_elts(
				nzvals, 0, nzcount,
				nzoffs_buf);
	if (ans_nzcount == 0)        /* all values in 'leaf' are zeros */
		return R_NilValue;
	if (ans_nzcount == nzcount)  /* all values in 'leaf' are nonzeros */
		return leaf;  /* no-op */

	ans_nzoffs = PROTECT(NEW_INTEGER(ans_nzcount));
	_copy_selected_ints(INTEGER(nzoffs), nzoffs_buf, ans_nzcount,
			    INTEGER(ans_nzoffs));
	ans_nzvals = PROTECT(allocVector(TYPEOF(nzvals), ans_nzcount));
	_copy_selected_Rsubvec_elts(nzvals, 0, nzoffs_buf, ans_nzvals);
	ans = zip_leaf(ans_nzoffs, ans_nzvals);
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * _coerce_leaf()
 *
 * Returns a leaf whose nzcount is guaranteed to not exceed the nzcount of
 * the input leaf.
 */

/* Note that, in practice, _coerce_leaf() is always called to actually
   change the type of 'leaf', so the code below does not bother to check
   for the (trivial) no-op case. */
SEXP _coerce_leaf(SEXP leaf, SEXPTYPE new_Rtype, int *warn, int *nzoffs_buf)
{
	SEXP nzoffs, nzvals, ans_nzvals, ans;

	unzip_leaf(leaf, &nzoffs, &nzvals);
	ans_nzvals = PROTECT(_coerceVector2(nzvals, new_Rtype, warn));
	ans = PROTECT(zip_leaf(nzoffs, ans_nzvals));
	/* The above coercion can introduce zeros in 'ans_nzvals' e.g. when
	   going from double/complex to int/raw. We need to remove them. */
	if (_coercion_can_introduce_zeros(TYPEOF(nzvals), new_Rtype))
		ans = _remove_zeros_from_leaf(ans, nzoffs_buf);
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * _subassign_leaf_with_Rvector()
 *
 * 'index' must be an integer vector containing valid zero-based indices
 * (a.k.a. offsets) into 'leaf' seen as dense. They are expected to be
 * already sorted in strictly ascending order.
 * 'Rvector' must be a vector (atomic or list) of the same length as 'index'.
 * Returns a leaf whose nzcount is guaranteed to not exceed
 * min(length(leaf) + length(index), INT_MAX). Note that the zeros in 'Rvector'
 * will be injected in the returned leaf. This function does NOT remove them!
 */

SEXP _subassign_leaf_with_Rvector(SEXP leaf, SEXP index, SEXP Rvector)
{
	SEXP nzoffs, nzvals;
	int nzcount = unzip_leaf(leaf, &nzoffs, &nzvals);
	SEXPTYPE Rtype = TYPEOF(nzvals);

	CopyRVectorElt_FUNType copy_Rvector_elt_FUN =
		_select_copy_Rvector_elt_FUN(Rtype);
	CopyRVectorElts_FUNType copy_Rvector_elts_FUN =
		_select_copy_Rvector_elts_FUN(Rtype);
	if (copy_Rvector_elt_FUN == NULL || copy_Rvector_elts_FUN == NULL)
		error("SparseArray internal error in "
		      "_subassign_leaf_with_Rvector():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));

	if (TYPEOF(Rvector) != Rtype)
		error("SparseArray internal error in "
		      "_subassign_leaf_with_Rvector():\n"
		      "    'leaf' and 'Rvector' have different types");

	int index_len = LENGTH(index);
	if (LENGTH(Rvector) != index_len)
		error("SparseArray internal error in "
		      "_subassign_leaf_with_Rvector():\n"
		      "    'index' and 'Rvector' have different lengths");
	if (index_len == 0)
		return leaf;  /* no-op */

	const int *nzoffs1_p = INTEGER(nzoffs);
	const int *offs2_p = INTEGER(index);
	int ans_nzcount = 0, k1 = 0, k2 = 0;
	while (k1 < nzcount && k2 < index_len) {
		if (*nzoffs1_p < *offs2_p) {
			nzoffs1_p++;
			k1++;
		} else if (*nzoffs1_p > *offs2_p) {
			offs2_p++;
			k2++;
		} else {
			/* *nzoffs1_p == *offs2_p */
			nzoffs1_p++;
			k1++;
			offs2_p++;
			k2++;
		}
		ans_nzcount++;
	}
	if (k1 < nzcount) {
		ans_nzcount += nzcount - k1;
	} else if (k2 < index_len) {
		ans_nzcount += index_len - k2;
	}

	SEXP ans_nzoffs, ans_nzvals;
	SEXP ans = PROTECT(_alloc_and_unzip_leaf(Rtype, ans_nzcount,
						 &ans_nzoffs, &ans_nzvals));
	nzoffs1_p = INTEGER(nzoffs);
	offs2_p = INTEGER(index);
	int *ans_nzoffs_p = INTEGER(ans_nzoffs);
	int k = 0;
	k1 = k2 = 0;
	while (k1 < nzcount && k2 < index_len) {
		if (*nzoffs1_p < *offs2_p) {
			*ans_nzoffs_p = *nzoffs1_p;
			copy_Rvector_elt_FUN(nzvals, (R_xlen_t) k1,
					     ans_nzvals, (R_xlen_t) k);
			nzoffs1_p++;
			k1++;
		} else if (*nzoffs1_p > *offs2_p) {
			*ans_nzoffs_p = *offs2_p;
			copy_Rvector_elt_FUN(Rvector, (R_xlen_t) k2,
					     ans_nzvals, (R_xlen_t) k);
			offs2_p++;
			k2++;
		} else {
			/* *nzoffs1_p == *offs2_p */
			*ans_nzoffs_p = *offs2_p;
			copy_Rvector_elt_FUN(Rvector, (R_xlen_t) k2,
					     ans_nzvals, (R_xlen_t) k);
			nzoffs1_p++;
			k1++;
			offs2_p++;
			k2++;
		}
		ans_nzoffs_p++;
		k++;
	}
	int n;
	if (k1 < nzcount) {
		n = nzcount - k1;
		memcpy(ans_nzoffs_p, nzoffs1_p, sizeof(int) * n);
		copy_Rvector_elts_FUN(nzvals, (R_xlen_t) k1,
				      ans_nzvals, (R_xlen_t) k,
				      (R_xlen_t) n);
	} else if (k2 < index_len) {
		n = index_len - k2;
		memcpy(ans_nzoffs_p, offs2_p, sizeof(int) * n);
		copy_Rvector_elts_FUN(Rvector, (R_xlen_t) k2,
				      ans_nzvals, (R_xlen_t) k,
				      (R_xlen_t) n);
	}
	UNPROTECT(1);
	return ans;
}

