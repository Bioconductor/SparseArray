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
 * _expand_leaf()
 */

/* Do NOT use when 'nzcount' is 0!
   Always use an R_NilValue to represent an empty leaf. See leaf_utils.h */
SEXP _alloc_leaf(SEXPTYPE Rtype, int nzcount)
{
	if (nzcount == 0)
		error("SparseArray internal error in _alloc_leaf():\n"
		      "    nzcount == 0");
	SEXP nzvals = PROTECT(allocVector(Rtype, nzcount));
	SEXP nzoffs = PROTECT(NEW_INTEGER(nzcount));
	SEXP ans = zip_leaf(nzvals, nzoffs);
	UNPROTECT(2);
	return ans;
}

SEXP _alloc_and_unzip_leaf(SEXPTYPE Rtype, int nzcount,
		SEXP *nzvals, SEXP *nzoffs)
{
	SEXP ans = PROTECT(_alloc_leaf(Rtype, nzcount));
	unzip_leaf(ans, nzvals, nzoffs);
	UNPROTECT(1);
	return ans;
}

/* Does NOT work for 'Rtype' == STRSXP or VECSXP. */
SEXP _make_leaf_from_bufs(SEXPTYPE Rtype,
		const void *nzvals_buf, const int *nzoffs_buf, int buf_len)
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
					    &ans_vals, &ans_offs));
	memcpy(DATAPTR(ans_vals), nzvals_buf, Rtype_size * buf_len);
	memcpy(INTEGER(ans_offs), nzoffs_buf, sizeof(int) * buf_len);
	UNPROTECT(1);
	return ans;
}

/* Assumes that 'out_Rvector' is long enough, has the right type,
   and is already filled with zeros e.g. it was created with
   _new_Rvector0(TYPEOF(leaf), dim0). */
void _expand_leaf(SEXP leaf, SEXP out_Rvector, R_xlen_t out_offset)
{
	SEXP nzvals, nzoffs;
	unzip_leaf(leaf, &nzvals, &nzoffs);  /* ignore returned nzcount */
	if (nzvals == R_NilValue) {
		/* lacunar leaf */
		_set_selected_Rsubvec_elts_to_one(out_Rvector, out_offset,
					INTEGER(nzoffs), LENGTH(nzoffs));

	} else {
		/* standard leaf */
		_copy_Rvector_elts_to_offsets(nzvals, INTEGER(nzoffs),
					out_Rvector, out_offset);
	}
	return;
}


/****************************************************************************
 * _make_leaf_from_Rsubvec()
 */

/* 'nzoffs_buf' must be of length 'subvec_len' (or more).
   The returned leaf can be lacunar. */
SEXP _make_leaf_from_Rsubvec(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *selection_buf, int avoid_copy_if_all_nonzeros)
{
	/* 'n' will always be >= 0 and <= subvec_len. */
	int n = _collect_offsets_of_nonzero_Rsubvec_elts(
				Rvector, subvec_offset, subvec_len,
				selection_buf);
	if (n == 0)
		return R_NilValue;

	SEXP ans_nzoffs = PROTECT(NEW_INTEGER(n));
	memcpy(INTEGER(ans_nzoffs), selection_buf, sizeof(int) * n);

	if (OK_TO_MAKE_LACUNAR_LEAVES) {
		int all_ones = _all_selected_Rsubvec_elts_equal_one(Rvector,
					subvec_offset, selection_buf, n);
		if (all_ones) {
			SEXP ans = zip_leaf(R_NilValue, ans_nzoffs);
			UNPROTECT(1);
			return ans;
		}
	}

	if (avoid_copy_if_all_nonzeros &&
	    subvec_offset == 0 && n == XLENGTH(Rvector))
	{
		/* The full 'Rvector' contains no zeros and can be reused
		   as-is without the need to copy its nonzero elements to
		   a new SEXP. */
		SEXP ans = zip_leaf(Rvector, ans_nzoffs);
		UNPROTECT(1);
		return ans;
	}

	SEXP ans_nzvals = PROTECT(
		_subset_Rsubvec(Rvector, subvec_offset, selection_buf, n)
	);
	SEXP ans = zip_leaf(ans_nzvals, ans_nzoffs);
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * _INPLACE_remove_zeros_from_leaf()
 * _INPLACE_turn_into_lacunar_leaf_if_all_ones()
 */

/* Do NOT use on a NULL or lacunar leaf. */
SEXP _INPLACE_remove_zeros_from_leaf(SEXP leaf, int *selection_buf)
{
	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
	/* 'n' will always be >= 0 and <= nzcount. */
	int n = _collect_offsets_of_nonzero_Rsubvec_elts(
				nzvals, 0, nzcount, selection_buf);

	if (n == 0)	      /* all values in 'nzvals' are zeros */
		return R_NilValue;

	if (n == nzcount)     /* all values in 'nzvals' are nonzeros */
		return leaf;  /* no-op */

	/* Shrink 'nzvals'. */
	SEXP new_nzvals = PROTECT(_subset_Rsubvec(nzvals, 0, selection_buf, n));
	replace_leaf_nzvals(leaf, new_nzvals);
	UNPROTECT(1);

	/* Shrink 'nzoffs'. */
	SEXP new_nzoffs = PROTECT(_subset_Rsubvec(nzoffs, 0, selection_buf, n));
	replace_leaf_nzoffs(leaf, new_nzoffs);
	UNPROTECT(1);
	return leaf;
}

/* Do NOT use on a NULL or lacunar leaf. */
void _INPLACE_turn_into_lacunar_leaf_if_all_ones(SEXP leaf)
{
	SEXP nzvals = get_leaf_nzvals(leaf);
	int nzcount = LENGTH(nzvals);
	int all_ones = _all_Rsubvec_elts_equal_one(nzvals, 0, nzcount);
	if (all_ones)
		replace_leaf_nzvals(leaf, R_NilValue);
	return;
}


/****************************************************************************
 * _coerce_leaf()
 *
 * Returns a leaf whose nzcount is guaranteed to not exceed the nzcount of
 * the input leaf.
 */

static SEXP coerce_lacunar_leaf(SEXP leaf, SEXPTYPE new_Rtype)
{
	if (new_Rtype != STRSXP && new_Rtype != VECSXP)
		return leaf;  /* no-op */
	error("SparseArray internal error in coerce_lacunar_leaf():"
	      "    coercing a lacunar leaf to \"character\" or \"double\" "
	      "    is not supported yet");
}

/* Note that, in practice, _coerce_leaf() is always called to actually
   change the type of 'leaf', so the code below does not bother to check
   for the (trivial) no-op case. */
SEXP _coerce_leaf(SEXP leaf, SEXPTYPE new_Rtype, int *warn, int *selection_buf)
{
	SEXP nzvals, nzoffs;
	unzip_leaf(leaf, &nzvals, &nzoffs);  /* ignore returned nzcount */
	if (nzvals == R_NilValue)  /* lacunar leaf */
		return coerce_lacunar_leaf(leaf, new_Rtype);
	SEXP ans_nzvals = PROTECT(_coerceVector2(nzvals, new_Rtype, warn));
	SEXP ans = PROTECT(zip_leaf(ans_nzvals, nzoffs));
	/* The above coercion can introduce zeros in 'ans_nzvals' e.g. when
	   going from double/complex to int/raw. We need to remove them. */
	if (_coercion_can_introduce_zeros(TYPEOF(nzvals), new_Rtype))
		ans = _INPLACE_remove_zeros_from_leaf(ans, selection_buf);
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * _subassign_leaf_with_Rvector()
 */

/* Do NOT use on a NULL leaf. Can be used on a lacunar leaf.
   'index' must be an integer vector containing valid zero-based indices
   (a.k.a. offsets) into 'leaf' seen as dense. They are expected to be
   already sorted in strictly ascending order.
   'Rvector' must be a vector (atomic or list) of the same length as 'index'.
   Returns a leaf POSSIBLY CONTAMINATED WITH ZEROS in its 'nzvals'.
   More precisely: the zeros in 'Rvector' will be injected in the 'nzvals'
   of the returned leaf. This function does NOT remove them!
   'nzcount' of the returned leaf is guaranteed to not exceed
   min(nzcount(leaf) + length(index), INT_MAX).
   Will NEVER return a NULL. Can ONLY return a lacunar leaf if input leaf
   is already lacunar **and** 'index' has length 0 (no-op). */
SEXP _subassign_leaf_with_Rvector(SEXP leaf, SEXP index, SEXP Rvector)
{
	int index_len = LENGTH(index);
	if (LENGTH(Rvector) != index_len)
		error("SparseArray internal error in "
		      "_subassign_leaf_with_Rvector():\n"
		      "    'index' and 'Rvector' have different lengths");
	if (index_len == 0)
		return leaf;  /* no-op */

	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
	SEXPTYPE Rtype = TYPEOF(Rvector);
	if (nzvals != R_NilValue && TYPEOF(nzvals) != Rtype)
		error("SparseArray internal error in "
		      "_subassign_leaf_with_Rvector():\n"
		      "    'Rvector' and 'leaf' have different types");

	/* Compute 'ans_nzcount'. */
	const int *nzoffs_p = INTEGER(nzoffs);
	const int *index_p = INTEGER(index);
	int ans_nzcount = 0, k1 = 0, k2 = 0;
	while (k1 < nzcount && k2 < index_len) {
		if (*nzoffs_p < *index_p) {
			nzoffs_p++;
			k1++;
		} else if (*nzoffs_p > *index_p) {
			index_p++;
			k2++;
		} else {
			/* *nzoffs_p == *index_p */
			nzoffs_p++;
			k1++;
			index_p++;
			k2++;
		}
		ans_nzcount++;
	}
	if (k1 < nzcount) {
		ans_nzcount += nzcount - k1;
	} else if (k2 < index_len) {
		ans_nzcount += index_len - k2;
	}

	CopyRVectorElt_FUNType copy_Rvector_elt_FUN =
		_select_copy_Rvector_elt_FUN(Rtype);
	CopyRVectorElts_FUNType copy_Rvector_elts_FUN =
		_select_copy_Rvector_elts_FUN(Rtype);
	if (copy_Rvector_elt_FUN == NULL || copy_Rvector_elts_FUN == NULL)
		error("SparseArray internal error in "
		      "_subassign_leaf_with_Rvector():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));

	SEXP ans_nzvals, ans_nzoffs;
	SEXP ans = PROTECT(_alloc_and_unzip_leaf(Rtype, ans_nzcount,
						 &ans_nzvals, &ans_nzoffs));

	/* Fill 'ans_nzvals' and 'ans_nzoffs'. */
	nzoffs_p = INTEGER(nzoffs);
	index_p = INTEGER(index);
	int *ans_nzoffs_p = INTEGER(ans_nzoffs);
	int k = 0;
	k1 = k2 = 0;
	while (k1 < nzcount && k2 < index_len) {
		if (*nzoffs_p < *index_p) {
			*ans_nzoffs_p = *nzoffs_p;
			if (nzvals == R_NilValue) {
				/* lacunar leaf */
				_set_selected_Rsubvec_elts_to_one(
						     ans_nzvals, 0, &k, 1);
			} else {
				copy_Rvector_elt_FUN(nzvals, (R_xlen_t) k1,
						     ans_nzvals, (R_xlen_t) k);
			}
			nzoffs_p++;
			k1++;
		} else if (*nzoffs_p > *index_p) {
			*ans_nzoffs_p = *index_p;
			copy_Rvector_elt_FUN(Rvector, (R_xlen_t) k2,
					     ans_nzvals, (R_xlen_t) k);
			index_p++;
			k2++;
		} else {
			/* *nzoffs_p == *index_p */
			*ans_nzoffs_p = *index_p;
			copy_Rvector_elt_FUN(Rvector, (R_xlen_t) k2,
					     ans_nzvals, (R_xlen_t) k);
			nzoffs_p++;
			k1++;
			index_p++;
			k2++;
		}
		ans_nzoffs_p++;
		k++;
	}
	int n;
	if (k1 < nzcount) {
		n = nzcount - k1;
		memcpy(ans_nzoffs_p, nzoffs_p, sizeof(int) * n);
		if (nzvals == R_NilValue) {
			/* lacunar leaf */
			_set_Rsubvec_to_one(ans_nzvals, (R_xlen_t) k, n);
		} else {
			copy_Rvector_elts_FUN(nzvals, (R_xlen_t) k1,
					      ans_nzvals, (R_xlen_t) k,
					      (R_xlen_t) n);
		}
	} else if (k2 < index_len) {
		n = index_len - k2;
		memcpy(ans_nzoffs_p, index_p, sizeof(int) * n);
		copy_Rvector_elts_FUN(Rvector, (R_xlen_t) k2,
				      ans_nzvals, (R_xlen_t) k,
				      (R_xlen_t) n);
	}
	UNPROTECT(1);
	return ans;
}

