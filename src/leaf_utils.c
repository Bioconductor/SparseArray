/****************************************************************************
 *                     Basic manipulation of SVT leaves                     *
 ****************************************************************************/
#include "leaf_utils.h"

#include "S4Vectors_interface.h"

#include "Rvector_utils.h"
#include "coerceVector2.h"
#include "SparseVec_subassignment.h"

#include <string.h>  /* for memcpy() */


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
		_fill_Rvector_subset_with_ones(out_Rvector,
					INTEGER(nzoffs), LENGTH(nzoffs),
					out_offset);

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
 * _make_leaf_from_Rvector_block()
 * _make_naleaf_from_Rvector_block()
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
	if (_all_elts_equal_one(Rtype, shared_nzval, 1))
		return _make_lacunar_leaf(nzoffs);
	SEXP nzvals = PROTECT(allocVector(Rtype, LENGTH(nzoffs)));
	_fill_Rvector_with_val(nzvals, shared_nzval);
	SEXP ans = zip_leaf(nzvals, nzoffs, 0);
	UNPROTECT(1);
	return ans;
}

/* Each of 'nzvals_p' and 'nzoffs_p' must be a pointer to an array of
   length 'nzcount'. 'nzvals_p' is **trusted** to not contain any zeros (this
   is NOT checked).
   The returned leaf can be lacunar. */
SEXP _make_leaf_from_two_arrays(SEXPTYPE Rtype,
		const void *nzvals_p, const int *nzoffs_p, int nzcount)
{
	if (nzcount == 0)
		return R_NilValue;
	SEXP ans_nzoffs = PROTECT(NEW_INTEGER(nzcount));
	memcpy(INTEGER(ans_nzoffs), nzoffs_p, sizeof(int) * nzcount);
	SEXP ans_nzvals;
	if (Rtype == STRSXP) {
		/* Lacunar leaves of Rtype STRSXP are not supported yet
		   (they will be in the future). */
		ans_nzvals = PROTECT(NEW_CHARACTER(nzcount));
		for (int k = 0; k < nzcount; k++) {
			SEXP nzval = STRING_ELT((SEXP) nzvals_p, k);
			SET_STRING_ELT(ans_nzvals, k, nzval);
		}
	} else if (Rtype == VECSXP) {
		/* Lacunar leaves of Rtype VECSXP are not supported. */
		ans_nzvals = PROTECT(NEW_LIST(nzcount));
		for (int k = 0; k < nzcount; k++) {
			SEXP nzval = VECTOR_ELT((SEXP) nzvals_p, k);
			SET_VECTOR_ELT(ans_nzvals, k, nzval);
		}
	} else {
		size_t Rtype_size = _get_Rtype_size(Rtype);
		if (Rtype_size == 0)
			error("SparseArray internal error in "
			      "_make_leaf_from_two_arrays():\n    type "
			      "\"%s\" is not supported", type2char(Rtype));
		int all_ones = _all_elts_equal_one(Rtype, nzvals_p, nzcount);
		if (all_ones) {
			SEXP ans = _make_lacunar_leaf(ans_nzoffs);
			UNPROTECT(1);
			return ans;
		}
		ans_nzvals = PROTECT(allocVector(Rtype, nzcount));
		memcpy(DATAPTR(ans_nzvals), nzvals_p, Rtype_size * nzcount);
	}
	SEXP ans = zip_leaf(ans_nzvals, ans_nzoffs, 0);
	UNPROTECT(2);
	return ans;
}

static SEXP make_leaf_from_Rvector_subset(SEXP Rvector,
		const int *selection, int selection_len,
		R_xlen_t selection_offset, int avoid_copy_if_all_selected)
{
	if (selection_len == 0)
		return R_NilValue;

	SEXP ans_nzoffs = PROTECT(NEW_INTEGER(selection_len));
	memcpy(INTEGER(ans_nzoffs), selection, sizeof(int) * selection_len);

	int all_ones = _Rvector_subset_is_filled_with_ones(Rvector,
				selection, selection_len, selection_offset);
	if (all_ones) {
		SEXP ans = _make_lacunar_leaf(ans_nzoffs);
		UNPROTECT(1);
		return ans;
	}

	if (avoid_copy_if_all_selected &&
	    selection_offset == 0 && selection_len == XLENGTH(Rvector) &&
	    ATTRIB(Rvector) == R_NilValue)
	{
		/* The full 'Rvector' is selected so can be reused as-is
		   with no need to copy the selected elements to a new SEXP. */
		SEXP ans = zip_leaf(Rvector, ans_nzoffs, 0);
		UNPROTECT(1);
		return ans;
	}

	SEXP ans_nzvals = PROTECT(
	  _subset_Rvector(Rvector, selection, selection_len, selection_offset)
	);
	SEXP ans = zip_leaf(ans_nzvals, ans_nzoffs, 0);
	UNPROTECT(2);
	return ans;
}


/* 'selection_buf' must be of length 'block_len' (at least).
   The returned leaf can be lacunar. */
SEXP _make_leaf_from_Rvector_block(
		SEXP Rvector, R_xlen_t block_offset, int block_len,
		int *selection_buf, int avoid_copy_if_all_nonzeros)
{
	/* 'selection_len' will always be >= 0 and <= block_len. */
	int selection_len = _collect_offsets_of_nonzero_elts_in_Rvector_block(
				Rvector, block_offset, block_len,
				selection_buf);
	return make_leaf_from_Rvector_subset(Rvector,
				selection_buf, selection_len, block_offset,
				avoid_copy_if_all_nonzeros);
}

/* 'selection_buf' must be of length 'block_len' (at least).
   The returned leaf can be lacunar. */
SEXP _make_naleaf_from_Rvector_block(
		SEXP Rvector, R_xlen_t block_offset, int block_len,
		int *selection_buf, int avoid_copy_if_all_nonNAs)
{
	/* 'selection_len' will always be >= 0 and <= block_len. */
	int selection_len = _collect_offsets_of_nonNA_elts_in_Rvector_block(
				Rvector, block_offset, block_len,
				selection_buf);
	return make_leaf_from_Rvector_subset(Rvector,
				selection_buf, selection_len, block_offset,
				avoid_copy_if_all_nonNAs);
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
	int all_ones = _Rvector_is_filled_with_ones(nzvals);
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
	SEXP new_nzoffs = PROTECT(_subset_Rvector(nzoffs, selection, n, 0));
	replace_leaf_nzoffs(leaf, new_nzoffs);
	UNPROTECT(1);

	/* Shrink 'nzvals'. */
	int all_ones = _Rvector_subset_is_filled_with_ones(nzvals,
							   selection, n, 0);
	if (all_ones) {
		replace_leaf_nzvals(leaf, R_NilValue);
		return 2;
	}
	SEXP new_nzvals = PROTECT(_subset_Rvector(nzvals, selection, n, 0));
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
	int new_nzcount = _collect_offsets_of_nonzero_elts_in_Rvector_block(
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
	int new_nzcount = _collect_offsets_of_nonNA_elts_in_Rvector_block(
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
	_copy_Rvector_subset(nzvals, order_buf, 0, new_nzvals);
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
	if (!IS_STRSXP_OR_VECSXP(new_Rtype))
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
		} else if (ret == 1) {
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
		} else if (ret == 1) {
			_INPLACE_turn_into_lacunar_leaf_if_all_ones(ans);
		}
	}
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * _subassign_leaf_with_Rvector_block()
 * _subassign_leaf_with_Rvector_subset()
 * _subassign_leaf_with_Rvector_xsubset()
 * _subassign_leaf_with_leaf()
 */

static const int *get_offs0(SEXP offs, int n, int dim0)
{
	if (offs == R_NilValue) {
		if (n != dim0)
			error("SparseArray internal error in get_offs0():\n"
			      "    'offs' is NULL but 'n != dim0'");
		return NULL;
	}
	if (n != LENGTH(offs))
		error("SparseArray internal error in get_offs0():\n"
		      "    n != LENGTH(offs)");
	return INTEGER(offs);
}

/* TODO: Maybe move this to SparseVec.c next to _alloc_buf_SparseVec(). */
static int alloc_SV_nzvals_if_needed(SparseVec *sv)
{
	if (sv->nzvals != NULL)
		return 0;
	size_t Rtype_size = _get_Rtype_size(sv->Rtype);
	if (Rtype_size == 0)
		error("SparseArray internal error in "
		      "alloc_SV_nzvals_if_needed():\n"
		      "    type \"%s\" is not supported", type2char(sv->Rtype));
	sv->nzvals = R_alloc(sv->len, Rtype_size);
	return 1;
}

/* Can be used on a NULL or lacunar leaf.
   'offs' must be NULL or an array of 'n' offsets (non-negative integers)
   that are sorted in strictly ascending order. The last offset in the
   array must be < 'buf_sv->len'.
   About the length of 'Rvector':
   - If 'block_offset != 0' then its length must be >= 'block_offset + n'.
   - If 'block_offset == 0' then 'Rvector' can be of any length (except 0
     unless 'n' is also 0) and will be recycled if its length is < 'n'. */
SEXP _subassign_leaf_with_Rvector_block(SEXP leaf, SEXP offs, int n,
		SEXP Rvector, R_xlen_t block_offset, SparseVec *buf_sv)
{
	const int *offs0 = get_offs0(offs, n, buf_sv->len);
	if (leaf == R_NilValue) {
		/* Taking the shortcut is not really needed because using
		   _write_Rvector_block_to_SV() can handle this case. It's
		   just that _make_leaf_from_Rvector_block() or
		   _make_naleaf_from_Rvector_block() don't
		   need 'buf_sv->nzvals'. This gives us the opportunity to
		   use a 'buf_sv' where 'buf_sv->nzvals' is not allocated
		   yet (i.e. is set to NULL) and to delay allocation until
		   it's actually needed.
		   Note that there's not measurable speed difference
		   between taking the shortcut or not, despite the shortcut
		   copying the nzvals data only once (directly from 'Rvector'
		   to the returned leaf) instead of twice (first by
		   _write_Rvector_block_to_SV() from 'Rvector'
		   to 'buf_sv->nzvals', then by SV2leaf()
		   from 'buf_sv->nzvals' to the returned leaf). */
		int use_shortcut = offs0 == NULL &&
				   (block_offset != 0 || LENGTH(Rvector) >= n);
		if (use_shortcut) {
			SEXP ans;
			if (buf_sv->na_background) {
				ans = _make_naleaf_from_Rvector_block(
						Rvector, block_offset, n,
						buf_sv->nzoffs, 1);
			} else {
				ans = _make_leaf_from_Rvector_block(
						Rvector, block_offset, n,
						buf_sv->nzoffs, 1);
			}
			return ans;
		}
		alloc_SV_nzvals_if_needed(buf_sv);
		_write_Rvector_block_to_SV(Rvector, block_offset,
					   offs0, n, buf_sv);
	} else {
		alloc_SV_nzvals_if_needed(buf_sv);
		const SparseVec sv = leaf2SV(leaf, buf_sv->Rtype,
					     buf_sv->len,
					     buf_sv->na_background);
		int neffrep;
		if (n == buf_sv->len) {
			/* Full replacement.
			   We don't need 'offs0' even if it's not NULL. */
			neffrep = _subassign_full_SV_with_Rvector_block(&sv,
					     Rvector, block_offset, buf_sv);
		} else {
			/* 'offs0' guaranteed to be != NULL. */
			neffrep = _subassign_SV_with_Rvector_block(&sv,
					     offs0, n,
					     Rvector, block_offset, buf_sv);
		}
		//printf("n = %d / neffrep = %d\n", n, neffrep);
		if (neffrep == 0)
			return leaf;  /* no-op */
	}
	return SV2leaf(buf_sv);
}

/* Can be used on a NULL or lacunar leaf.
   Both 'offs' and 'selection' must be arrays of 'n' offsets (non-negative
   integers).
   The offsets in 'offs' must be sorted in strictly ascending order. The
   last offset in the array must be < 'buf_sv->len'.
   The offsets in 'selection' must be < 'LENGTH(Rvector)'. They are typically
   unsorted. */
SEXP _subassign_leaf_with_Rvector_subset(SEXP leaf, const int *offs, int n,
		SEXP Rvector, const int *selection, SparseVec *buf_sv)
{
	if (leaf == R_NilValue) {
		_write_Rvector_subset_to_SV(Rvector, selection, offs, n,
					    buf_sv);
	} else {
		const SparseVec sv = leaf2SV(leaf, buf_sv->Rtype,
					     buf_sv->len,
					     buf_sv->na_background);
		int neffrep = _subassign_SV_with_Rvector_subset(&sv,
					     offs, n,
					     Rvector, selection, buf_sv);
		//printf("n = %d / neffrep = %d\n", n, neffrep);
		if (neffrep == 0)
			return leaf;  /* no-op */
	}
	return SV2leaf(buf_sv);
}

/* Can be used on a NULL or lacunar leaf.
   'offs' and 'xselection' must have length 'n' (they cannot be NULL). */
SEXP _subassign_leaf_with_Rvector_xsubset(SEXP leaf, const int *offs, int n,
		SEXP Rvector, const R_xlen_t *xselection, SparseVec *buf_sv)
{
	error("_subassign_leaf_with_Rvector_xsubset() is not ready yet");
	return R_NilValue;
}

/* Can be used on NULL or lacunar leaves.
   'offs' must be NULL or an array of 'n' offsets (non-negative integers)
   that are sorted in strictly ascending order. The last offset in the
   array must be < 'buf_sv->len'. */
SEXP _subassign_leaf_with_leaf(SEXP leaf1, SEXP offs, int n,
		SEXP leaf2, SparseVec *buf_sv)
{
	const int *offs0 = get_offs0(offs, n, buf_sv->len);
	/* 'n == buf_sv->len' indicates full replacement. */
	if ((leaf1 == R_NilValue && leaf2 == R_NilValue) || n == buf_sv->len)
		return leaf2;
	/* 'offs0' guaranteed to be != NULL. */
	const SparseVec sv1 = leaf2SV(leaf1, buf_sv->Rtype,
				      buf_sv->len,
				      buf_sv->na_background);
	const SparseVec sv2 = leaf2SV(leaf2, buf_sv->Rtype,
				      n,
				      buf_sv->na_background);
	int neffrep = _subassign_SV_with_SV(&sv1, offs0, &sv2, buf_sv);
	//printf("n = %d / neffrep = %d\n", n, neffrep);
	if (neffrep == 0)
		return leaf1;  /* no-op */
	return SV2leaf(buf_sv);
}

