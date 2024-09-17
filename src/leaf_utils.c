/****************************************************************************
 *                     Basic manipulation of SVT leaves                     *
 ****************************************************************************/
#include "leaf_utils.h"

#include "S4Vectors_interface.h"

#include "Rvector_utils.h"
#include "coerceVector2.h"

#include <string.h>  /* for memcpy() */


/* --- .Call ENTRY POINT --- */
SEXP C_lacunar_mode_is_on(void)
{
	return ScalarLogical(LACUNAR_MODE_IS_ON);
}


/****************************************************************************
 * _alloc_leaf()
 * _alloc_and_unzip_leaf()
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
	SEXP ans = zip_leaf(nzvals, nzoffs, 0);
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

/* Assumes that 'out_Rvector' is long enough, has the right type,
   and is already filled with zeros e.g. it was created with
   _new_Rvector0(TYPEOF(leaf), dim0). */
void _expand_leaf(SEXP leaf, SEXP out_Rvector, R_xlen_t out_offset)
{
	SEXP nzvals, nzoffs;
	unzip_leaf(leaf, &nzvals, &nzoffs);  /* ignore returned nzcount */
	if (nzvals == R_NilValue) {  /* lacunar leaf */
		_set_selected_Rsubvec_elts_to_one(out_Rvector, out_offset,
					INTEGER(nzoffs), LENGTH(nzoffs));

	} else {  /* standard leaf */
		_copy_Rvector_elts_to_offsets(nzvals, INTEGER(nzoffs),
					out_Rvector, out_offset);
	}
	return;
}


/****************************************************************************
 * _make_lacunar_leaf()
 * _make_leaf_with_single_shared_nzval()
 * _make_leaf_from_two_arrays()
 * _make_leaf_from_Rsubvec()
 * _make_naleaf_from_Rsubvec()
 */

SEXP _make_lacunar_leaf(SEXP nzoffs)
{
	return zip_leaf(R_NilValue, nzoffs, 0);
}

/* When 'Rtype' is STRSXP or VECSXP, 'shared_nzval' must be an SEXP.
   Otherwise, it must be a pointer to an int, double, Rcomplex, or Rbyte. */
SEXP _make_leaf_with_single_shared_nzval(SEXPTYPE Rtype,
		const void *shared_nzval, SEXP nzoffs)
{
	if (LACUNAR_MODE_IS_ON && _all_elts_equal_one(Rtype, shared_nzval, 1))
		return _make_lacunar_leaf(nzoffs);
	SEXP nzvals = PROTECT(allocVector(Rtype, LENGTH(nzoffs)));
	_set_Rvector_elts_to_val(nzvals, shared_nzval);
	SEXP ans = zip_leaf(nzvals, nzoffs, 0);
	UNPROTECT(1);
	return ans;
}

/* Does NOT work if 'Rtype' is STRSXP or VECSXP.
   Each of 'nzvals_p' and 'nzoffs_p' must be a pointer to an array of length
   'nzcount'. 'nzvals_p' is **trusted** to not contain any zeros. This is NOT
   checked! The returned leaf can be lacunar. */
SEXP _make_leaf_from_two_arrays(SEXPTYPE Rtype,
		const void *nzvals_p, const int *nzoffs_p, int nzcount)
{
	if (nzcount == 0)
		return R_NilValue;

	size_t Rtype_size = _get_Rtype_size(Rtype);
	if (Rtype_size == 0)
		error("SparseArray internal error in "
		      "_make_leaf_from_two_arrays():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));

	SEXP ans_nzoffs = PROTECT(NEW_INTEGER(nzcount));
	memcpy(INTEGER(ans_nzoffs), nzoffs_p, sizeof(int) * nzcount);

	if (LACUNAR_MODE_IS_ON) {
		int all_ones = _all_elts_equal_one(Rtype, nzvals_p, nzcount);
		if (all_ones) {
			SEXP ans = _make_lacunar_leaf(ans_nzoffs);
			UNPROTECT(1);
			return ans;
		}
	}
	SEXP ans_nzvals = PROTECT(allocVector(Rtype, nzcount));
	memcpy(DATAPTR(ans_nzvals), nzvals_p, Rtype_size * nzcount);
	SEXP ans = zip_leaf(ans_nzvals, ans_nzoffs, 0);
	UNPROTECT(2);
	return ans;
}

static SEXP make_leaf_from_selected_Rsubvec_elts(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		const int *selection, int n, int avoid_copy_if_all_selected)
{
	if (n == 0)
		return R_NilValue;

	SEXP ans_nzoffs = PROTECT(NEW_INTEGER(n));
	memcpy(INTEGER(ans_nzoffs), selection, sizeof(int) * n);

	if (LACUNAR_MODE_IS_ON) {
		int all_ones = _all_selected_Rsubvec_elts_equal_one(Rvector,
					subvec_offset, selection, n);
		if (all_ones) {
			SEXP ans = _make_lacunar_leaf(ans_nzoffs);
			UNPROTECT(1);
			return ans;
		}
	}

	if (avoid_copy_if_all_selected &&
	    subvec_offset == 0 && n == XLENGTH(Rvector))
	{
		/* The full 'Rvector' is selected so can be reused as-is
		   with no need to copy the selected elements to a new SEXP. */
		SEXP ans = zip_leaf(Rvector, ans_nzoffs, 0);
		UNPROTECT(1);
		return ans;
	}

	SEXP ans_nzvals = PROTECT(
		_subset_Rsubvec(Rvector, subvec_offset, selection, n)
	);
	SEXP ans = zip_leaf(ans_nzvals, ans_nzoffs, 0);
	UNPROTECT(2);
	return ans;
}


/* 'selection_buf' must be of length 'subvec_len' (at least).
   The returned leaf can be lacunar. */
SEXP _make_leaf_from_Rsubvec(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *selection_buf, int avoid_copy_if_all_nonzeros)
{
	/* 'n' will always be >= 0 and <= subvec_len. */
	int n = _collect_offsets_of_nonzero_Rsubvec_elts(
				Rvector, subvec_offset, subvec_len,
				selection_buf);
	return make_leaf_from_selected_Rsubvec_elts(
				Rvector, subvec_offset, subvec_len,
				selection_buf, n, avoid_copy_if_all_nonzeros);
}

/* 'selection_buf' must be of length 'subvec_len' (at least).
   The returned leaf can be lacunar. */
SEXP _make_naleaf_from_Rsubvec(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *selection_buf, int avoid_copy_if_all_nonNAs)
{
	/* 'n' will always be >= 0 and <= subvec_len. */
	int n = _collect_offsets_of_nonNA_Rsubvec_elts(
				Rvector, subvec_offset, subvec_len,
				selection_buf);
	return make_leaf_from_selected_Rsubvec_elts(
				Rvector, subvec_offset, subvec_len,
				selection_buf, n, avoid_copy_if_all_nonNAs);
}


/****************************************************************************
 * _INPLACE_turn_into_lacunar_leaf_if_all_ones()
 * _INPLACE_remove_zeros_from_leaf()
 * _INPLACE_remove_NAs_from_leaf()
 * _INPLACE_order_leaf_by_nzoff()
 */

/* Do NOT use on a NULL or lacunar leaf.
   Returns 1 if leaf was effectively turned into lacunar leaf. Otherwise
   returns 0 (no-op). */
int _INPLACE_turn_into_lacunar_leaf_if_all_ones(SEXP leaf)
{
	SEXP nzvals = get_leaf_nzvals(leaf);
	int nzcount = LENGTH(nzvals);
	int all_ones = _all_Rsubvec_elts_equal_one(nzvals, 0, nzcount);
	if (all_ones)
		replace_leaf_nzvals(leaf, R_NilValue);
	return all_ones;
}

/* Do NOT use on a NULL or lacunar leaf.
   Returns a code between 0 and 3, with the following meaning:
     0: Nothing is selected (i.e. n = 0). So everything **would** need to be
        removed from the input leaf. However, in this case, please note that
        the function does NOT touch the input leaf and it's the responsibility
        of the caller to replace the original leaf with a NULL leaf.
     1: No-op i.e. everything is selected (so nothing needs to be removed).
        Note that unlike with code 0 above, the caller does not need to do
        anything. However, not that the caller could still try to turn the
        original leaf into a lacunar leaf if that wasn't attempted earlier,
        by calling _INPLACE_turn_into_lacunar_leaf_if_all_ones() on it.
     2: The selection is neither nothing or everything, and the modified leaf
        is lacunar.
     3: The selection is neither nothing or everything, and the modified leaf
        is not lacunar (and it cannot be turned into one). */
static int INPLACE_extract_selection_from_leaf(SEXP leaf,
		const int *selection, int n)
{
	if (n == 0)
		return 0;

	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
	if (n == nzcount)
		return 1;  /* no-op */

	/* Shrink 'nzoffs'. */
	SEXP new_nzoffs = PROTECT(_subset_Rsubvec(nzoffs, 0, selection, n));
	replace_leaf_nzoffs(leaf, new_nzoffs);
	UNPROTECT(1);

	/* Shrink 'nzvals'. */
	if (LACUNAR_MODE_IS_ON) {
		int all_ones =
			_all_selected_Rsubvec_elts_equal_one(nzvals, 0,
							     selection, n);
		if (all_ones) {
			replace_leaf_nzvals(leaf, R_NilValue);
			return 2;
		}
	}
	SEXP new_nzvals = PROTECT(_subset_Rsubvec(nzvals, 0, selection, n));
	replace_leaf_nzvals(leaf, new_nzvals);
	UNPROTECT(1);
	return 3;
}

/* Do NOT use on a NULL or lacunar leaf.
   See INPLACE_extract_selection_from_leaf() above for the returned code. */
int _INPLACE_remove_zeros_from_leaf(SEXP leaf, int *selection_buf)
{
	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
	/* 'new_nzcount' will always be >= 0 and <= nzcount. */
	int new_nzcount = _collect_offsets_of_nonzero_Rsubvec_elts(
					nzvals, 0, nzcount, selection_buf);
	return INPLACE_extract_selection_from_leaf(leaf,
					selection_buf, new_nzcount);
}

/* Do NOT use on a NULL or lacunar leaf.
   See INPLACE_extract_selection_from_leaf() above for the returned code. */
int _INPLACE_remove_NAs_from_leaf(SEXP leaf, int *selection_buf)
{
	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
	/* 'new_nzcount' will always be >= 0 and <= nzcount. */
	int new_nzcount = _collect_offsets_of_nonNA_Rsubvec_elts(
					nzvals, 0, nzcount, selection_buf);
	return INPLACE_extract_selection_from_leaf(leaf,
					selection_buf, new_nzcount);
}

/* Do NOT use on a NULL leaf. Can be used on a lacunar leaf. */
void _INPLACE_order_leaf_by_nzoff(SEXP leaf, int *order_buf,
		unsigned short int *rxbuf1, int *rxbuf2)
{
	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
	for (int k = 0; k < nzcount; k++)
		order_buf[k] = k;
	int ret = sort_ints(order_buf, nzcount, INTEGER(nzoffs), 0, 1,
			    rxbuf1, rxbuf2);
	/* Note that ckecking the value returned by sort_ints() is not really
	   necessary here because sort_ints() should never fail when 'rxbuf1'
	   and 'rxbuf2' are supplied (see implementation of _sort_ints() in
	   S4Vectors/src/sort_utils.c for the details). We perform this check
	   nonetheless just to be on the safe side in case the implementation
	   of sort_ints() changes in the future. */
	if (ret < 0)
		error("SparseArray internal error in "
		      "_INPLACE_order_leaf_by_nzoff():\n"
		      "    sort_ints() returned an error");
	if (ret == 0)
		return;  /* no-op */

	SEXP new_nzoffs = PROTECT(NEW_INTEGER(nzcount));
	_copy_selected_int_elts(INTEGER(nzoffs),
				order_buf, nzcount,
				INTEGER(new_nzoffs));
	replace_leaf_nzoffs(leaf, new_nzoffs);
	UNPROTECT(1);

	if (nzvals == R_NilValue)  /* lacunar leaf */
		return;

	/* regular leaf */
	SEXP new_nzvals = PROTECT(allocVector(TYPEOF(nzvals), nzcount));
	_copy_selected_Rsubvec_elts(nzvals, 0, order_buf, new_nzvals);
	replace_leaf_nzvals(leaf, new_nzvals);
	UNPROTECT(1);
	return;
}


/****************************************************************************
 * _coerce_leaf()
 * _coerce_naleaf()
 *
 * Returns a leaf whose nzcount is guaranteed to not exceed the nzcount of
 * the input leaf.
 */

static SEXP coerce_lacunar_leaf(SEXP leaf, SEXPTYPE new_Rtype)
{
	if (new_Rtype != STRSXP && new_Rtype != VECSXP)
		return leaf;  /* no-op */
	error("SparseArray internal error in coerce_lacunar_leaf():\n"
	      "    coercing a lacunar leaf to \"character\" or \"list\" "
	      "is not supported yet");
}

/* Note that, in practice, _coerce_leaf() is always called to actually
   change the type of 'leaf', so the code below does not bother to check
   for the (trivial) no-op case. */
SEXP _coerce_leaf(SEXP leaf, SEXPTYPE new_Rtype, int *warn,
		  int *selection_buf)
{
	SEXP nzvals, nzoffs;
	unzip_leaf(leaf, &nzvals, &nzoffs);  /* ignore returned nzcount */
	if (nzvals == R_NilValue)  /* lacunar leaf */
		return coerce_lacunar_leaf(leaf, new_Rtype);
	/* standard leaf */
	SEXP ans_nzvals = PROTECT(_coerceVector2(nzvals, new_Rtype, warn));
	SEXP ans = PROTECT(zip_leaf(ans_nzvals, nzoffs, 0));
	/* The above coercion can introduce zeros in 'ans_nzvals' e.g. when
	   going from double/complex to int/raw. We need to remove them. */
	if (_coercion_can_introduce_zeros(TYPEOF(nzvals), new_Rtype)) {
		int ret = _INPLACE_remove_zeros_from_leaf(ans, selection_buf);
		if (ret == 0) {
			ans = R_NilValue;
		} if (ret == 1 && LACUNAR_MODE_IS_ON) {
			_INPLACE_turn_into_lacunar_leaf_if_all_ones(ans);
		}
	}
	UNPROTECT(2);
	return ans;
}

SEXP _coerce_naleaf(SEXP leaf, SEXPTYPE new_Rtype, int *warn,
		    int *selection_buf)
{
	SEXP nzvals, nzoffs;
	unzip_leaf(leaf, &nzvals, &nzoffs);  /* ignore returned nzcount */
	if (nzvals == R_NilValue)  /* lacunar leaf */
		return coerce_lacunar_leaf(leaf, new_Rtype);
	/* standard leaf */
	int w = 0;
	SEXP ans_nzvals = PROTECT(_coerceVector2(nzvals, new_Rtype, &w));
	SEXP ans = PROTECT(zip_leaf(ans_nzvals, nzoffs, 0));
	if (w) {
		*warn = 1;
		int ret = _INPLACE_remove_NAs_from_leaf(ans, selection_buf);
		if (ret == 0) {
			ans = R_NilValue;
		} if (ret == 1 && LACUNAR_MODE_IS_ON) {
			_INPLACE_turn_into_lacunar_leaf_if_all_ones(ans);
		}
	}
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
	if (copy_Rvector_elt_FUN == NULL)
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
			copy_Rvector_elt_FUN(nzvals, (R_xlen_t) k1,
					     ans_nzvals, (R_xlen_t) k);
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
		if (nzvals == R_NilValue) {  /* lacunar leaf */
			_set_Rsubvec_elts_to_one(ans_nzvals, (R_xlen_t) k,
						 (R_xlen_t) n);
		} else {  /* standard leaf */
			_copy_Rvector_elts(nzvals, (R_xlen_t) k1,
					   ans_nzvals, (R_xlen_t) k,
					   (R_xlen_t) n);
		}
	} else if (k2 < index_len) {
		n = index_len - k2;
		memcpy(ans_nzoffs_p, index_p, sizeof(int) * n);
		_copy_Rvector_elts(Rvector, (R_xlen_t) k2,
				   ans_nzvals, (R_xlen_t) k,
				   (R_xlen_t) n);
	}
	UNPROTECT(1);
	return ans;
}

