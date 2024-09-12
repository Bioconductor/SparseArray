/****************************************************************************
 ****************************************************************************
 **									   **
 **             Subassignment (`[<-`) to a SparseArray object              **
 **									   **
 ****************************************************************************
 ****************************************************************************/
#include "SparseArray_subassignment.h"

#include "argcheck_utils.h"
#include "OPBufTree.h"
#include "thread_control.h"  /* for which_max() */
#include "leaf_utils.h"

#include <limits.h>  /* for INT_MAX */
//#include <time.h>


/* Copied from S4Arrays/src/array_selection.h */
#define INVALID_COORD(coord, maxcoord) \
	((coord) == NA_INTEGER || (coord) < 1 || (coord) > (maxcoord))


static inline R_xlen_t get_Lidx(SEXP Lindex, long long atid_lloff)
{
	R_xlen_t Lidx;

	if (IS_INTEGER(Lindex)) {
		int i = INTEGER(Lindex)[atid_lloff];
		if (i == NA_INTEGER || i < 1)
			error("'Lindex' contains invalid linear indices");
		Lidx = (R_xlen_t) i;
	} else {
		double x = REAL(Lindex)[atid_lloff];
		/* ISNAN(): True for *both* NA and NaN. See <R_ext/Arith.h> */
		if (ISNAN(x) || x < 1 || x >= 1.00 + R_XLEN_T_MAX)
			error("'Lindex' contains invalid linear indices");
		Lidx = (R_xlen_t) x;
	}
	return Lidx;
}


/****************************************************************************
 * Basic manipulation of "extended leaves"
 *
 * An "extended leaf" is used to temporarily attach a subset of the incoming
 * data (represented by 'Mindex' and 'vals', or by 'Lindex' and 'vals') to
 * an SVT leaf.
 *
 * There are 3 types of extended leaves:
 *
 * - type 1: A standalone Incoming Data Subset (IDS). An IDS is simply a
 *           set of offsets w.r.t. 'Mindex' (or 'Lindex') and 'vals'.
 *           These offsets get stored in an IntAE or LLongAE buffer placed
 *           behind an external pointer, and are referred to as "atid" offsets
 *           (offsets along the incoming data).
 *           Note that using an IntAE buffer would be ok for now because we're
 *           not dealing with _long_ incoming data yet. However, this will
 *           change when we start supporting _long_ incoming data e.g. when
 *           C_subassign_SVT_by_Lindex() will get passed a _long_ linear index.
 *
 * - type 2: Just a regular leaf (possibly lacunar) so not really "extended"
 *           in that case.
 *
 * - type 3: A regular leaf with an IDS on it. This is represented by a
 *           list of length 3: the 2 list elements of a regular leaf (nzvals
 *           and nzoffs) + the IDS.
 *
 * IMPORTANT NOTE: We don't allow the length of an IDS to be more than INT_MAX
 * at the moment. This is because we use sort_ints() in compute_offs_order()
 * below to sort a vector of 'IDS_len' integers and sort_ints() only handles
 * a vector of length <= INT_MAX!
 * Note however that 'IDS_len' > INT_MAX can't happen at the moment anyway
 * because 'IDS_len' is necessarily <= 'nrow(Mindex)' which is guaranteed
 * to be <= INT_MAX. However, this will change when we start supporting
 * _long_ incoming data e.g. when C_subassign_SVT_by_Lindex() is called
 * with a _long_ linear index (Lindex). Then it will be possible that more
 * than INT_MAX incoming values land on the same SVT leaf but only in some
 * crazy and rather unlikely situations. More precisely this will be possible
 * only if the supplied Lindex is _long_ and contains duplicates. Like here:
 *
 *     svt[sample(nrow(svt), 3e9, replace=TRUE)] <- 2.5
 *
 * where 3e9 incoming values are landing on the SVT leaf associated with
 * the first column of the sparse matrix! A very atypical situation.
 */
#include "S4Vectors_interface.h"

typedef SEXP (*NewIDS_FUNType)(void);


/****************************************************************************
 * REC_postprocess_SVT_using_[M|L]index()
 *
 * This implements the 2nd pass of C_subassign_SVT_by_[M|L]index().
 */

typedef struct sort_bufs_t {
	int *order;
	unsigned short int *rxbuf1;
	int *rxbuf2;
	int *offs;
} SortBufs;

/* All buffers are made of length 'max_IDS_len' except 'sort_bufs.offs'
   which we must make of length 'max(max_IDS_len, max_postsubassign_nzcount)'
   so that we can use it in the call to _INPLACE_remove_zeros_from_leaf()
   in the subassign_xleaf3_with_offval_pairs() function below. */
static SortBufs alloc_sort_bufs(int max_IDS_len, int max_postsubassign_nzcount)
{
	SortBufs sort_bufs;
	int offs_len;

	sort_bufs.order = (int *) R_alloc(max_IDS_len, sizeof(int));
	sort_bufs.rxbuf1 = (unsigned short int *)
			R_alloc(max_IDS_len, sizeof(unsigned short int));
	sort_bufs.rxbuf2 = (int *) R_alloc(max_IDS_len, sizeof(int));
	offs_len = max_postsubassign_nzcount > max_IDS_len ?
			max_postsubassign_nzcount : max_IDS_len;
	sort_bufs.offs = (int *) R_alloc(offs_len, sizeof(int));
	return sort_bufs;
}

static void compute_offs_order(SortBufs *sort_bufs, int n)
{
	int k, ret;

	for (k = 0; k < n; k++)
		sort_bufs->order[k] = k;
	ret = sort_ints(sort_bufs->order, n, sort_bufs->offs, 0, 1,
			sort_bufs->rxbuf1, sort_bufs->rxbuf2);
	/* Note that ckecking the value returned by sort_ints() is not really
	   necessary here because sort_ints() should never fail when 'rxbuf1'
	   and 'rxbuf2' are supplied (see implementation of _sort_ints() in
	   S4Vectors/src/sort_utils.c for the details). We perform this check
	   nonetheless just to be on the safe side in case the implementation
	   of sort_ints() changes in the future. */
	if (ret < 0)
		error("SparseArray internal error in compute_offs_order():\n"
		      "    sort_ints() returned an error");
	return;
}

/* Returns number of offsets after removal of the duplicates. */
static int remove_offs_dups(int *order_buf, int n, const int *offs)
{
	int *p1, k2;
	const int *p2;

	if (n <= 1)
		return n;
	p1 = order_buf;
	for (k2 = 1, p2 = p1 + 1; k2 < n; k2++, p2++) {
		if (offs[*p1] != offs[*p2])
			p1++;
		*p1 = *p2;
	}
	return p1 - order_buf + 1;
}


/****************************************************************************
 * subassign_leaf_by_Lindex()
 *
 * This is the 1D case and it needs special treatment.
 */

/* 'Lindex' and 'vals' are assumed to have the same length. This length
   is assumed to be >= 1 and <= INT_MAX.
   Returns a set of offset/value pairs sorted by strictly ascending offset.
   It is returned as a list of 2 parallel vectors: an integer vector of
   strictly sorted offsets and a subset of 'vals'. Their common length
   is >= 1 and <= length(vals).
   Note that this is the "leaf representation", that is, the representation
   that we use for a 1D SVT. With an important gotcha: in the case of these
   off/val pairs the values are allowed to be zero! Also let's keep in mind
   that they are conceptually really different: in this case the 2 parallel
   vectors in the returned list are the 'index' and 'value' vectors of a
   subassignment operation that we will perform later on. They do NOT
   represent a 1D SVT!
   Anyways, we still use the "leaf representation" because it's convenient
   e.g. this will allow us to use things like _INPLACE_remove_zeros_from_leaf()
   later on it etc.. */
static SEXP make_offval_pairs_from_Lindex_vals(SEXP Lindex, SEXP vals,
		int dim0, SortBufs *sort_bufs)
{
	int nvals = LENGTH(vals);  /* we know 'length(vals)' is <= INT_MAX */
	/* Walk along the incoming data. */
	for (int atid_off = 0; atid_off < nvals; atid_off++) {
		R_xlen_t Lidx = get_Lidx(Lindex, atid_off);
		if (Lidx > dim0)
			error("subassignment subscript contains "
			      "invalid indices");
		sort_bufs->offs[atid_off] = Lidx - 1;
	}
	compute_offs_order(sort_bufs, nvals);
	int num_pairs = remove_offs_dups(sort_bufs->order, nvals,
					 sort_bufs->offs);
	SEXP ans_offs = PROTECT(NEW_INTEGER(num_pairs));
	_copy_selected_int_elts(sort_bufs->offs, sort_bufs->order, num_pairs,
				INTEGER(ans_offs));
	SEXP ans_vals = PROTECT(allocVector(TYPEOF(vals), num_pairs));
	_copy_selected_Rsubvec_elts(vals, 0, sort_bufs->order, ans_vals);
	/* Use the "leaf representation" even though this is NOT a 1D SVT!
	   See above. */
	SEXP ans = PROTECT(zip_leaf(ans_vals, ans_offs, 0));
	UNPROTECT(3);
	return ans;
}

/* 'Lindex' and 'vals' are assumed to have the same nonzero length.
   The returned leaf can be NULL or lacunar. */
static SEXP subassign_leaf_by_Lindex(SEXP leaf, int dim0, int na_background,
		SEXP Lindex, SEXP vals)
{
	if (na_background)
		error("subassignment of 1D NaArray objects "
		      "is not supported yet");
	R_xlen_t nvals = XLENGTH(vals);
	if (nvals > INT_MAX)
		error("assigning more than INT_MAX values to "
		      "a monodimensional SVT_SparseArray object "
		      "is not supported");
	size_t worst_nzcount;
	if (leaf == R_NilValue) {
		worst_nzcount = nvals;
	} else {
		int nzcount = get_leaf_nzcount(leaf);
		worst_nzcount = nzcount + nvals;
		if (worst_nzcount > dim0)
			worst_nzcount = dim0;
	}
	SortBufs sort_bufs = alloc_sort_bufs((int) nvals, (int) worst_nzcount);
	SEXP offval_pairs = PROTECT(
		make_offval_pairs_from_Lindex_vals(Lindex, vals,
						   dim0, &sort_bufs)
	);
	if (leaf != R_NilValue) {
		offval_pairs = PROTECT(
			_subassign_leaf_with_Rvector(leaf,
					get_leaf_nzoffs(offval_pairs),
					get_leaf_nzvals(offval_pairs))
		);
	}
	/* We use the "leaf representation" for 'offval_pairs' so it
	   should be safe to use _INPLACE_remove_zeros_from_leaf() on it.
	   Also we've made sure that 'sort_bufs.offs' is big enough for this
	   (its length is at least 'worst_nzcount'). */
	int ret = _INPLACE_remove_zeros_from_leaf(offval_pairs,
						  sort_bufs.offs);
	if (ret == 0) {
		offval_pairs = R_NilValue;
	} if (ret == 1 && LACUNAR_MODE_IS_ON) {
		_INPLACE_turn_into_lacunar_leaf_if_all_ones(offval_pairs);
	}
	UNPROTECT(leaf != R_NilValue ? 2 : 1);
	return offval_pairs;
}


/****************************************************************************
 * subassign_leaf_by_OPBuf()
 */

static void init_idx0_to_k_map(int *idx0_to_k_map, const int *idx0s, int nelt)
{
	for (int k = 0; k < nelt; k++)
		idx0_to_k_map[idx0s[k]] = k;
	return;
}

static void reset_idx0_to_k_map(int *idx0_to_k_map, const int *idx0s, int nelt)
{
	for (int k = 0; k < nelt; k++)
		idx0_to_k_map[idx0s[k]] = -1;
	return;
}

/*
static void print_idx0_to_k_map(const int *idx0_to_k_map, int dim0)
{
	printf("idx0_to_k_map:");
	for (int i = 0; i < dim0; i++)
		printf(" %4d", idx0_to_k_map[i]);
	printf("\n");
	return;
}
*/

/* TODO: Maybe add this to OPBufTree.h as inline functions. */
#define	GET_LOFF(Loffs, xLoffs, k) \
	((Loffs) != NULL ? (R_xlen_t) ((Loffs)[(k)]) : (xLoffs)[(k)])
#define	GET_OPBUF_LOFF(opbuf, k) GET_LOFF(opbuf->Loffs, opbuf->xLoffs, k)

/* TODO: Move all this to Rvector_utils.h. */
static inline int Rvector_elt_is_int0(SEXP Rvector, R_xlen_t i)
{
	return INTEGER(Rvector)[i] == int0;
}
static inline int Rvector_elt_is_intNA(SEXP Rvector, R_xlen_t i)
{
	return INTEGER(Rvector)[i] == NA_INTEGER;
}

static inline int Rvector_elt_is_double0(SEXP Rvector, R_xlen_t i)
{
	return REAL(Rvector)[i] == double0;
}
static inline int Rvector_elt_is_doubleNA(SEXP Rvector, R_xlen_t i)
{
	return R_IsNA(REAL(Rvector)[i]);
}

static inline int Rvector_elt_is_Rcomplex0(SEXP Rvector, R_xlen_t i)
{
	const Rcomplex *z = COMPLEX(Rvector) + i;
	return z->r == Rcomplex0.r && z->i == Rcomplex0.i;
}
static inline int Rvector_elt_is_RcomplexNA(SEXP Rvector, R_xlen_t i)
{
	const Rcomplex *z = COMPLEX(Rvector) + i;
	return R_IsNA(z->r) || R_IsNA(z->i);
}

static inline int Rvector_elt_is_Rbyte0(SEXP Rvector, R_xlen_t i)
{
	return RAW(Rvector)[i] == Rbyte0;
}

static inline int Rvector_elt_is_Rstring0(SEXP Rvector, R_xlen_t i)
{
	return IS_EMPTY_CHARSXP(STRING_ELT(Rvector, i));
}
static inline int Rvector_elt_is_RstringNA(SEXP Rvector, R_xlen_t i)
{
	return STRING_ELT(Rvector, i) == NA_STRING;
}

static inline int Rvector_elt_is_R_NilValue(SEXP Rvector, R_xlen_t i)
{
	return VECTOR_ELT(Rvector, i) == R_NilValue;
}

typedef int (*RVectorEltIsZero_FUNType)(SEXP Rvector, R_xlen_t i);

static RVectorEltIsZero_FUNType select_Rvector_elt_is_zero_FUN(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case INTSXP: case LGLSXP: return Rvector_elt_is_int0;
	    case REALSXP:             return Rvector_elt_is_double0;
	    case CPLXSXP:             return Rvector_elt_is_Rcomplex0;
	    case RAWSXP:              return Rvector_elt_is_Rbyte0;
	    case STRSXP:              return Rvector_elt_is_Rstring0;
	    case VECSXP:              return Rvector_elt_is_R_NilValue;
	}
	error("SparseArray internal error in "
	      "select_Rvector_elt_is_zero_FUN():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
}

static RVectorEltIsZero_FUNType select_Rvector_elt_is_NA_FUN(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case INTSXP: case LGLSXP: return Rvector_elt_is_intNA;
	    case REALSXP:             return Rvector_elt_is_doubleNA;
	    case CPLXSXP:             return Rvector_elt_is_RcomplexNA;
	    case STRSXP:              return Rvector_elt_is_RstringNA;
	}
	error("SparseArray internal error in "
	      "select_Rvector_elt_is_NA_FUN():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
}

static inline int same_INTEGER_vals(
		SEXP Rvector1, R_xlen_t i1,
		SEXP Rvector2, R_xlen_t i2)
{
	int val1 = Rvector1 == R_NilValue ? int1 : INTEGER(Rvector1)[i1];
	return val1 == INTEGER(Rvector2)[i2];
}

static inline int same_NUMERIC_vals(
		SEXP Rvector1, R_xlen_t i1,
		SEXP Rvector2, R_xlen_t i2)
{
	double val1 = Rvector1 == R_NilValue ? double1 : REAL(Rvector1)[i1];
	return val1 == REAL(Rvector2)[i2];
}

static inline int same_COMPLEX_vals(
		SEXP Rvector1, R_xlen_t i1,
		SEXP Rvector2, R_xlen_t i2)
{
	const Rcomplex *z1 = Rvector1 == R_NilValue ? &Rcomplex1
						    : COMPLEX(Rvector1) + i1;
	const Rcomplex *z2 = COMPLEX(Rvector2) + i2;
	return z1->r == z2->r && z1->i == z2->i;
}

static inline int same_RAW_vals(
		SEXP Rvector1, R_xlen_t i1,
		SEXP Rvector2, R_xlen_t i2)
{
	Rbyte val1 = Rvector1 == R_NilValue ? Rbyte1 : RAW(Rvector1)[i1];
	return val1 == RAW(Rvector2)[i2];
}

static inline int same_CHARACTER_vals(
		SEXP Rvector1, R_xlen_t i1,
		SEXP Rvector2, R_xlen_t i2)
{
	if (Rvector1 == R_NilValue)
		error("SparseArray internal error in same_CHARACTER_vals():\n"
		      "    lacunar leaf found in an SVT_SparseArray object "
		      "of type \"character\"");
	/* Compares the addresses, not the actual values. Doesn't matter as
	   long as our primary use case is covered. Primary use case is that
	       svt[Lindex] <- svt[Lindex]
	   is a no-op that triggers no copy. */
	return STRING_ELT(Rvector1, i1) == STRING_ELT(Rvector2, i2);
}

static inline int same_LIST_vals(
		SEXP Rvector1, R_xlen_t i1,
		SEXP Rvector2, R_xlen_t i2)
{
	if (Rvector1 == R_NilValue)
		error("SparseArray internal error in same_LIST_vals():\n"
		      "    lacunar leaf found in an SVT_SparseArray object "
		      "of type \"list\"");
	/* Compares the addresses, not the actual values. Doesn't matter as
	   long as our primary use case is covered. Primary use case is that
	       svt[Lindex] <- svt[Lindex]
	   is a no-op that triggers no copy. */
	return VECTOR_ELT(Rvector1, i1) == VECTOR_ELT(Rvector2, i2);
}

typedef int (*SameRVectorVals_FUNType)(SEXP Rvector1, R_xlen_t i1,
				       SEXP Rvector2, R_xlen_t i2);

static SameRVectorVals_FUNType select_same_Rvector_vals_FUN(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case INTSXP: case LGLSXP: return same_INTEGER_vals;
	    case REALSXP:             return same_NUMERIC_vals;
	    case CPLXSXP:             return same_COMPLEX_vals;
	    case RAWSXP:              return same_RAW_vals;
	    case STRSXP:              return same_CHARACTER_vals;
	    case VECSXP:              return same_LIST_vals;
	}
	return NULL;
}

static SEXP subassign_NULL_by_OPBuf(int dim0,
		const OPBuf *opbuf, SEXP vals,
		RVectorEltIsZero_FUNType Rvector_elt_is_zero_FUN,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN,
		int *idx0_order_buf, unsigned short int *rxbuf1, int *rxbuf2,
		int *idx0_to_k_map)
{
	int ans_nzcount = 0;
	for (int k = 0; k < opbuf->nelt; k++) {
		int idx0 = opbuf->idx0s[k];
		R_xlen_t Loff = GET_OPBUF_LOFF(opbuf, k);
		int val_is_zero = Rvector_elt_is_zero_FUN(vals, Loff);
		if (val_is_zero) {
			if (idx0_to_k_map[idx0] == -1)
				continue;
			idx0_to_k_map[idx0] = -1;
			ans_nzcount--;
		} else {
			if (idx0_to_k_map[idx0] == -1)
				ans_nzcount++;
			idx0_to_k_map[idx0] = k;
		}
	}
	if (ans_nzcount == 0)
		return R_NilValue;

	for (int k = 0; k < opbuf->nelt; k++)
		idx0_order_buf[k] = k;
	int ret = sort_ints(idx0_order_buf, opbuf->nelt, opbuf->idx0s, 0, 1,
			    rxbuf1, rxbuf2);
	/* Note that ckecking the value returned by sort_ints() is not really
	   necessary here because sort_ints() should never fail when 'rxbuf1'
	   and 'rxbuf2' are supplied (see implementation of _sort_ints() in
	   S4Vectors/src/sort_utils.c for the details). We perform this check
	   nonetheless just to be on the safe side in case the implementation
	   of sort_ints() changes in the future. */
	if (ret < 0)
		error("SparseArray internal error in "
		      "subassign_NULL_by_OPBuf():\n"
		      "    sort_ints() returned an error");

	SEXP ans_nzvals = PROTECT(allocVector(TYPEOF(vals), ans_nzcount));
	SEXP ans_nzoffs = PROTECT(NEW_INTEGER(ans_nzcount));
	int *ans_nzoffs_p = INTEGER(ans_nzoffs);
	ans_nzcount = 0;
/*
	for (int i = 0; i < dim0; i++) {
		int k1 = idx0_to_k_map[i];
		if (k1 == -1)
			continue;
		R_xlen_t Loff = GET_OPBUF_LOFF(opbuf, k1);
		copy_Rvector_elt_FUN(vals, Loff,
				     ans_nzvals, (R_xlen_t) ans_nzcount);
		ans_nzoffs_p[ans_nzcount] = i;
		ans_nzcount++;
	}
*/
	/* Walk on the (idx0,Loff) pairs in ascending 'idx0' order. */
	for (int k0 = 0; k0 < opbuf->nelt; k0++) {
		int k = idx0_order_buf[k0];
		int idx0 = opbuf->idx0s[k];
		if (idx0_to_k_map[idx0] != k)
			continue;
		R_xlen_t Loff = GET_OPBUF_LOFF(opbuf, k);
		copy_Rvector_elt_FUN(vals, Loff,
				     ans_nzvals, (R_xlen_t) ans_nzcount);
		ans_nzoffs_p[ans_nzcount] = idx0;
		ans_nzcount++;
	}
	SEXP ans = zip_leaf(ans_nzvals, ans_nzoffs, 1);
	UNPROTECT(2);
	return ans;
}

/* Returns -1 if subassignment is a no-op. */
static int compute_subassignment_nzcount(SEXP leaf, int dim0,
		const OPBuf *opbuf, SEXP vals,
		RVectorEltIsZero_FUNType Rvector_elt_is_zero_FUN,
		SameRVectorVals_FUNType same_Rvector_vals_FUN,
		int *idx0_to_k_map)
{
	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
	//print_idx0_to_k_map(idx0_to_k_map, dim0);
	int out_nzcount = 0;
	int k2 = 0, nzoff = INTEGER(nzoffs)[0];
	int is_noop = 1;
	for (int i = 0; i < dim0; i++) {
		int k1 = idx0_to_k_map[i];
		if (i != nzoff) {
			if (k1 == -1)
				continue;
			R_xlen_t Loff = GET_OPBUF_LOFF(opbuf, k1);
			int is_zero = Rvector_elt_is_zero_FUN(vals, Loff);
			if (is_zero) {
				idx0_to_k_map[i] = -1;
				continue;
			}
			out_nzcount++;
			is_noop = 0;
			continue;
		}
		if (k1 == -1) {
			out_nzcount++;
		} else {
			R_xlen_t Loff = GET_OPBUF_LOFF(opbuf, k1);
			int is_zero = Rvector_elt_is_zero_FUN(vals, Loff);
			if (is_zero) {
				is_noop = 0;
			} else {
				out_nzcount++;
				R_xlen_t Loff = GET_OPBUF_LOFF(opbuf, k1);
				if (!same_Rvector_vals_FUN(nzvals, k2,
							   vals, Loff))
				{
					is_noop = 0;
				}
			}
		}
		/* Move to next nzoffs[]. */
		k2++;
		nzoff = k2 < nzcount ? INTEGER(nzoffs)[k2] : -1;
	}
	if (is_noop && out_nzcount != nzcount)  /* sanity check */
		error("SparseArray internal error in "
		      "compute_subassignment_nzcount():\n"
		      "    leaf subassignment is a no-op "
		      "but nzcount(out_leaf) != nzcount(in_leaf)");
	return is_noop ? -1 : out_nzcount;
}

static void do_subassign_nonNULL_leaf_by_OPBuf(SEXP leaf, int dim0,
		const OPBuf *opbuf, SEXP vals,
		SEXP ans_nzvals, SEXP ans_nzoffs,
		RVectorEltIsZero_FUNType Rvector_elt_is_zero_FUN,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN,
		const int *idx0_to_k_map)
{
	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
	int *ans_nzoffs_p = INTEGER(ans_nzoffs);
	int ans_nzcount = 0;
	int k2 = 0, nzoff = INTEGER(nzoffs)[0];
	for (int i = 0; i < dim0; i++) {
		int k1 = idx0_to_k_map[i];
		if (i != nzoff) {
			if (k1 == -1)
				continue;
			R_xlen_t Loff = GET_OPBUF_LOFF(opbuf, k1);
			copy_Rvector_elt_FUN(vals, Loff,
					ans_nzvals, (R_xlen_t) ans_nzcount);
			ans_nzoffs_p[ans_nzcount] = i;
			ans_nzcount++;
			continue;
		}
		if (k1 == -1) {
			copy_Rvector_elt_FUN(nzvals, k2,
					ans_nzvals, (R_xlen_t) ans_nzcount);
			ans_nzoffs_p[ans_nzcount] = i;
			ans_nzcount++;
		} else {
			R_xlen_t Loff = GET_OPBUF_LOFF(opbuf, k1);
			int is_zero = Rvector_elt_is_zero_FUN(vals, Loff);
			if (!is_zero) {
				copy_Rvector_elt_FUN(vals, Loff,
					ans_nzvals, (R_xlen_t) ans_nzcount);
				ans_nzoffs_p[ans_nzcount] = i;
				ans_nzcount++;
			}
		}
		/* Move to next nzoffs[]. */
		k2++;
		nzoff = k2 < nzcount ? INTEGER(nzoffs)[k2] : -1;
	}
	if (ans_nzcount != LENGTH(ans_nzoffs))  /* sanity check */
		error("SparseArray internal error in "
		      "do_subassign_nonNULL_leaf_by_OPBuf():\n"
		      "    ans_nzcount != LENGTH(ans_nzoffs)");
	return;
}

/* 'leaf' cannot be R_NilValue. */
static SEXP subassign_nonNULL_leaf_by_OPBuf(SEXP leaf, int dim0,
		const OPBuf *opbuf, SEXP vals,
		RVectorEltIsZero_FUNType fun1,
		SameRVectorVals_FUNType fun2,
		CopyRVectorElt_FUNType fun3,
		int *idx0_to_k_map)
{
	int ans_nzcount = compute_subassignment_nzcount(leaf, dim0,
				opbuf, vals, fun1, fun2, idx0_to_k_map);
	if (ans_nzcount == -1)  /* no-op */
		return leaf;
	if (ans_nzcount == 0)
		return R_NilValue;
	SEXP ans_nzvals = PROTECT(allocVector(TYPEOF(vals), ans_nzcount));
	SEXP ans_nzoffs = PROTECT(NEW_INTEGER(ans_nzcount));
	do_subassign_nonNULL_leaf_by_OPBuf(leaf, dim0,
				opbuf, vals, ans_nzvals, ans_nzoffs,
				fun1, fun3, idx0_to_k_map);
	SEXP ans = zip_leaf(ans_nzvals, ans_nzoffs, 1);
	UNPROTECT(2);
	return ans;
}

static SEXP subassign_leaf_by_OPBuf(SEXP leaf, int dim0,
		const OPBuf *opbuf, SEXP vals,
		RVectorEltIsZero_FUNType fun1,
		SameRVectorVals_FUNType fun2,
		CopyRVectorElt_FUNType fun3,
		int *idx0_order_buf, unsigned short int *rxbuf1, int *rxbuf2,
		int *idx0_to_k_map)
{
	SEXP ans;
	/* PROTECT(ans) not necessary because reset_idx0_to_k_map() won't
	   trigger R's garbage collector. */
	if (leaf == R_NilValue) {
		ans = subassign_NULL_by_OPBuf(dim0, opbuf, vals,
					fun1, fun3,
					idx0_order_buf, rxbuf1, rxbuf2,
					idx0_to_k_map);
	} else {
		init_idx0_to_k_map(idx0_to_k_map, opbuf->idx0s, opbuf->nelt);
		ans = subassign_nonNULL_leaf_by_OPBuf(leaf, dim0, opbuf, vals,
					fun1, fun2, fun3, idx0_to_k_map);
	}
	reset_idx0_to_k_map(idx0_to_k_map, opbuf->idx0s, opbuf->nelt);
	return ans;
}


/****************************************************************************
 * C_subassign_SVT_by_Lindex()
 */

/* 'Lidx0' is trusted to be a non-NA value >= 0 and < 'dimcumprod[ndim - 1]'.
   Returns NULL if we didn't land anywhere. */
static OPBufTree *find_host_node_for_Lidx0(OPBufTree *opbuf_tree,
		R_xlen_t Lidx0,
		const int *dim, int ndim,
		const R_xlen_t *dimcumprod, int *idx0)
{
	for (int along = ndim - 1; along >= 1; along--) {
		R_xlen_t p = dimcumprod[along - 1];
		int i = Lidx0 / p;  /* always >= 0 and < 'dim[along]' */
		Lidx0 %= p;
		if (opbuf_tree->node_type == NULL_NODE)
			_alloc_OPBufTree_children(opbuf_tree, dim[along]);
		opbuf_tree = get_OPBufTree_child(opbuf_tree, i);
	}
	/* At this point:
	   - 'Lidx0' is guaranteed to be < 'dimcumprod[0]' (note that
	     'dimcumprod[0]' should always be = 'dim[0]' and <= INT_MAX);
	   - 'opbuf_tree' is guaranteed to be a node of type NULL_NODE or
	     LEAF_NODE. */
	*idx0 = (int) Lidx0;
	return opbuf_tree;
}

static int build_OPBufTree_from_Lindex1(OPBufTree *opbuf_tree, SEXP Lindex,
		const int *x_dim, int x_ndim,
		const R_xlen_t *dimcumprod)
{
	int max_outleaf_len = 0;
	int in_len = LENGTH(Lindex);
	R_xlen_t x_len = dimcumprod[x_ndim - 1];
	/* Walk along 'Lindex'. */
	for (int Loff = 0; Loff < in_len; Loff++) {
		R_xlen_t Lidx0;
		int ret = extract_long_idx0(Lindex, (R_xlen_t) Loff, x_len,
					    &Lidx0);
		if (ret < 0)
			return ret;
		int idx0;
		OPBufTree *host_node = find_host_node_for_Lidx0(
						opbuf_tree, Lidx0,
						x_dim, x_ndim,
						dimcumprod, &idx0);
		ret = _append_idx0Loff_to_host_node(host_node, idx0, Loff);
		if (ret < 0)
			return ret;
		if (ret > max_outleaf_len)
			max_outleaf_len = ret;
	}
	return max_outleaf_len;
}

static int build_OPBufTree_from_Lindex2(OPBufTree *opbuf_tree, SEXP Lindex,
		const int *x_dim, int x_ndim,
		const R_xlen_t *dimcumprod)
{
	error("build_OPBufTree_from_Lindex2() not ready yet");
	return 0;
}

static int build_OPBufTree_from_Lindex(OPBufTree *opbuf_tree, SEXP Lindex,
		const int *x_dim, int x_ndim,
		const R_xlen_t *dimcumprod)
{
	/* _free_OPBufTree(opbuf_tree) resets 'opbuf_tree->node_type'
	   to NULL_NODE. */
	_free_OPBufTree(opbuf_tree);
	return XLENGTH(Lindex) <= (R_xlen_t) INT_MAX ?
		build_OPBufTree_from_Lindex1(opbuf_tree, Lindex,
				x_dim, x_ndim, dimcumprod) :
		build_OPBufTree_from_Lindex2(opbuf_tree, Lindex,
				x_dim, x_ndim, dimcumprod);
}

/* Recursive tree traversal of 'opbuf_tree'. */
static SEXP REC_subassign_SVT_by_OPBufTree(OPBufTree *opbuf_tree,
		SEXP SVT, const int *dim, int ndim, SEXP vals,
		RVectorEltIsZero_FUNType fun1,
		SameRVectorVals_FUNType fun2,
		CopyRVectorElt_FUNType fun3,
		int *idx0_order_buf, unsigned short int *rxbuf1, int *rxbuf2,
		int *idx0_to_k_map, int pardim)
{
	if (opbuf_tree->node_type == NULL_NODE)
		return SVT;

	if (ndim == 1) {
		/* Both 'opbuf_tree' and 'SVT' are leaves. */
		OPBuf *opbuf = get_OPBufTree_leaf(opbuf_tree);
		SEXP ans = subassign_leaf_by_OPBuf(SVT, dim[0], opbuf, vals,
					fun1, fun2, fun3,
					idx0_order_buf, rxbuf1, rxbuf2,
					idx0_to_k_map);
		/* PROTECT not really necessary since neither _free_OPBufTree()
		   won't trigger R's garbage collector but this could change
		   someday so we'd better not take any risk. */
		PROTECT(ans);
		_free_OPBufTree(opbuf_tree);
		UNPROTECT(1);
		return ans;
	}

	/* Both 'opbuf_tree' and 'SVT' are inner nodes. */
	int n = get_OPBufTree_nchildren(opbuf_tree);  /* = dim[ndim - 1] */
	SEXP ans = PROTECT(NEW_LIST(n));
	int is_empty = 1;
	for (int i = 0; i < n; i++) {
		OPBufTree *child = get_OPBufTree_child(opbuf_tree, i);
		SEXP subSVT = SVT == R_NilValue ? R_NilValue
						: VECTOR_ELT(SVT, i);
		SEXP ans_elt = REC_subassign_SVT_by_OPBufTree(child,
					subSVT, dim, ndim - 1, vals,
					fun1, fun2, fun3,
					idx0_order_buf, rxbuf1, rxbuf2,
					idx0_to_k_map, pardim);
		if (ans_elt != R_NilValue) {
			PROTECT(ans_elt);
			SET_VECTOR_ELT(ans, i, ans_elt);
			UNPROTECT(1);
			is_empty = 0;
		}
	}
	UNPROTECT(1);
	return is_empty ? R_NilValue : ans;
}

/* --- .Call ENTRY POINT ---
   'Lindex' must be a numeric vector (integer or double), possibly a long one.
   NAs are not allowed (they'll trigger an error).
   'vals' must be a vector (atomic or list) of type 'x_type'. */
SEXP C_subassign_SVT_by_Lindex(
		SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP x_na_background,
		SEXP Lindex, SEXP vals)
{
	SEXPTYPE Rtype = _get_and_check_Rtype_from_Rstring(x_type,
				"C_subassign_SVT_by_Lindex", "x_type");
	if (TYPEOF(vals) != Rtype)
		error("SparseArray internal error in "
		      "C_subassign_SVT_by_Lindex():\n"
		      "    SVT_SparseArray object and 'vals' "
		      "must have the same type");

	int x_has_NAbg = _get_and_check_na_background(x_na_background,
				"C_subassign_SVT_by_Lindex", "x_na_background");

	if (!(IS_INTEGER(Lindex) || IS_NUMERIC(Lindex)))
		error("'Lindex' must be an integer or numeric vector");

	int x_ndim = LENGTH(x_dim);
	R_xlen_t nvals = XLENGTH(vals);
	if (XLENGTH(Lindex) != nvals)
		error("length(Lindex) != length(vals)");
	if (nvals == 0)
		return x_SVT;  /* no-op */

	RVectorEltIsZero_FUNType fun1;
	if (x_has_NAbg) {
		fun1 = select_Rvector_elt_is_NA_FUN(Rtype);
	} else {
		fun1 = select_Rvector_elt_is_zero_FUN(Rtype);
	}
	SameRVectorVals_FUNType fun2 = select_same_Rvector_vals_FUN(Rtype);
	CopyRVectorElt_FUNType fun3 = _select_copy_Rvector_elt_FUN(Rtype);

	int x_dim0 = INTEGER(x_dim)[0];
	if (x_ndim == 1)
		return subassign_leaf_by_Lindex(
				x_SVT, x_dim0, x_has_NAbg,
				Lindex, vals);

	/* 1st pass: Build the OPBufTree. */
	//clock_t t0 = clock();
	OPBufTree *opbuf_tree = _get_global_opbuf_tree();
	R_xlen_t *dimcumprod = (R_xlen_t *) R_alloc(x_ndim, sizeof(R_xlen_t));
	R_xlen_t p = 1;
	for (int along = 0; along < x_ndim; along++) {
		p *= INTEGER(x_dim)[along];
		dimcumprod[along] = p;
	}
	int max_outleaf_len =
		build_OPBufTree_from_Lindex(opbuf_tree, Lindex,
				INTEGER(x_dim), x_ndim, dimcumprod);
	if (max_outleaf_len < 0) {
		UNPROTECT(1);
		_bad_Lindex_error(max_outleaf_len);
	}

	//double dt = (1.0 * clock() - t0) * 1000.0 / CLOCKS_PER_SEC;
	//printf("1st pass: %2.3f ms\n", dt);

	//printf("max_outleaf_len = %d\n", max_outleaf_len);
	//_print_OPBufTree(opbuf_tree, 1);

	/* 2nd pass: Subset SVT by OPBufTree. */
	//t0 = clock();
	int *idx0_to_k_map = (int *) R_alloc(x_dim0, sizeof(int));
	for (int i = 0; i < x_dim0; i++)
		idx0_to_k_map[i] = -1;
	/* Three buffers needed by sort_ints(). */
	int *idx0_order_buf = (int *) R_alloc(max_outleaf_len, sizeof(int));
	unsigned short int *rxbuf1 = (unsigned short int *)
			R_alloc(max_outleaf_len, sizeof(unsigned short int));
	int *rxbuf2 = (int *) R_alloc(max_outleaf_len, sizeof(int));
	/* Get 1-based rank of biggest dimension (ignoring the 1st dim).
	   Parallel execution will be along that dimension. */
	int pardim = which_max(INTEGER(x_dim) + 1, x_ndim - 1) + 2;
	SEXP ans = REC_subassign_SVT_by_OPBufTree(opbuf_tree,
				x_SVT, INTEGER(x_dim), x_ndim, vals,
				fun1, fun2, fun3,
				idx0_order_buf, rxbuf1, rxbuf2,
				idx0_to_k_map, pardim);
	//dt = (1.0 * clock() - t0) * 1000.0 / CLOCKS_PER_SEC;
	//printf("2nd pass: %2.3f ms\n", dt);
	return ans;
}


/****************************************************************************
 * C_subassign_SVT_by_Mindex()
 */

static void check_Mindex_dim(SEXP Mindex, R_xlen_t nvals, int ndim,
		const char *what1, const char *what2, const char *what3)
{
	SEXP Mindex_dim = GET_DIM(Mindex);
	if (Mindex_dim == R_NilValue || LENGTH(Mindex_dim) != 2)
		error("'%s' must be a matrix", what1);
	if (!IS_INTEGER(Mindex))
		error("'%s' must be an integer matrix", what1);
	if (INTEGER(Mindex_dim)[0] != nvals)
		error("nrow(%s) != %s", what1, what2);
	if (INTEGER(Mindex_dim)[1] != ndim)
		error("ncol(%s) != %s", what1, what3);
	return;
}

/* --- .Call ENTRY POINT --- */
SEXP C_subassign_SVT_by_Mindex(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		SEXP Mindex, SEXP vals)
{
	SEXPTYPE Rtype = _get_and_check_Rtype_from_Rstring(x_type,
					"C_subassign_SVT_by_Mindex", "x_type");
	if (TYPEOF(vals) != Rtype)
		error("SparseArray internal error in "
		      "C_subassign_SVT_by_Mindex():\n"
		      "    SVT_SparseArray object and 'vals' "
		      "must have the same type");

	int x_ndim = LENGTH(x_dim);
	R_xlen_t nvals = XLENGTH(vals);
	check_Mindex_dim(Mindex, nvals, x_ndim,
			 "Mindex", "length(vals)", "length(dim(x))");
	if (nvals == 0)
		return x_SVT;  /* no-op */

	//RVectorEltIsZero_FUNType fun1 = select_Rvector_elt_is_zero_FUN(Rtype);
	//SameRVectorVals_FUNType fun2 = select_same_Rvector_vals_FUN(Rtype);
	//CopyRVectorElt_FUNType fun3 = _select_copy_Rvector_elt_FUN(Rtype);

	int x_dim0 = INTEGER(x_dim)[0];
	if (x_ndim == 1)
		return subassign_leaf_by_Lindex(x_SVT, x_dim0, 0, Mindex, vals);

	/* 1st pass: Build the OPBufTree. */
	error("C_subassign_SVT_by_Mindex() not ready yet");

	/* 2nd pass: Subset SVT by OPBufTree. */

	return R_NilValue;
}


/****************************************************************************
 * make_SVT_node()
 */

static SEXP shallow_copy_list(SEXP x)
{
	int x_len, i;
	SEXP ans;

	if (!isVectorList(x))  // IS_LIST() is broken
		error("SparseArray internal error in shallow_copy_list():\n"
		      "    'x' is not a list");
	x_len = LENGTH(x);
	ans = PROTECT(NEW_LIST(x_len));
	for (i = 0; i < x_len; i++)
		SET_VECTOR_ELT(ans, i, VECTOR_ELT(x, i));
	UNPROTECT(1);
	return ans;
}

/* 'SVT' must be R_NilValue or a list of length 'd' ('d' cannot be 0).
   Always returns a list of length 'd'. Can be a newly allocated list
   or 'SVT' itself. */
static inline SEXP make_SVT_node(SEXP SVT, int d, SEXP SVT0)
{
	if (d == 0)
		error("SparseArray internal error in make_SVT_node():\n"
		      "    d == 0");
	if (SVT == R_NilValue)
		return NEW_LIST(d);
	if (!isVectorList(SVT) || LENGTH(SVT) != d)
		error("SparseArray internal error in make_SVT_node():\n"
		      "    'SVT' is not R_NilValue or a list of length 'd'");
	/* Shallow copy **only** if 'SVT' == corresponding node in
	   original 'SVT0'. */
	if (SVT == SVT0)
		return shallow_copy_list(SVT);
	return SVT;
}


/****************************************************************************
 * C_subassign_SVT_with_short_Rvector()
 */

typedef struct left_bufs_t {
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;
	SEXP Rvector;
	int *offs;
	SEXP precomputed_leaf;
	int full_replacement;
} LeftBufs;

/* 'short_Rvector' must have a length >= 1.
   'dim0' must be a multiple of 'short_Rvector' length. */
static SEXP precompute_leaf_from_short_Rvector(
		int dim0, SEXP index0, SEXP short_Rvector,
		LeftBufs *left_bufs)
{
	left_bufs->full_replacement = 1;
	SEXP left_Rvector = left_bufs->Rvector;
	int short_len = LENGTH(short_Rvector);
	if (index0 == R_NilValue) {
		if (short_len == dim0) {
			left_Rvector = short_Rvector;
		} else {
			/* Copy a recycled version of 'short_Rvector'
			   to 'left_bufs->Rvector'. 'left_bufs->Rvector' is
			   of length 'dim0'. */
			for (int i1 = 0; i1 < dim0; i1++) {
				left_bufs->copy_Rvector_elt_FUN(short_Rvector,
						i1 % short_len,
						left_Rvector, i1);
			}
		}
	} else {
		for (int i1 = 0; i1 < dim0; i1++)
			left_bufs->offs[i1] = 0;
		/* Recycle and subassign 'short_Rvector' into 'left_Rvector'. */
		int d2 = LENGTH(index0);
		for (int i2 = 0; i2 < d2; i2++) {
			int coord = INTEGER(index0)[i2];
			if (INVALID_COORD(coord, dim0))
				error("subscript contains "
				      "out-of-bound indices or NAs");
			int i1 = coord - 1;
			left_bufs->copy_Rvector_elt_FUN(short_Rvector,
						i2 % short_len,
						left_Rvector, i1);
			left_bufs->offs[i1] = 1;
		}
		for (int i1 = 0; i1 < dim0; i1++) {
			if (left_bufs->offs[i1] == 0) {
				left_bufs->full_replacement = 0;
				break;
			}
		}
	}
	//printf("full_replacement=%d\n", left_bufs->full_replacement);
	return _make_leaf_from_Rsubvec(left_Rvector, 0, dim0,
				       left_bufs->offs,
				       left_bufs->full_replacement);
}

/* 'short_Rvector' must have a length >= 1.
   'dim0' must be a multiple of 'short_Rvector' length. */
static LeftBufs init_left_bufs(int dim0, SEXP index0, SEXP short_Rvector)
{
	LeftBufs left_bufs;
	SEXPTYPE Rtype;
	R_xlen_t short_len;
	SEXP leaf;

	Rtype = TYPEOF(short_Rvector);
	left_bufs.copy_Rvector_elt_FUN = _select_copy_Rvector_elt_FUN(Rtype);
	if (left_bufs.copy_Rvector_elt_FUN == NULL)
		error("SparseArray internal error in init_left_bufs():\n"
		      "    short Rvector has invalid type");

	short_len = XLENGTH(short_Rvector);
	if (short_len == 0 || LENGTH(index0) % short_len != 0)
		error("SparseArray internal error in init_left_bufs():\n"
		      "    invalid short Rvector length");

	left_bufs.offs = (int *) R_alloc(dim0, sizeof(int));
	left_bufs.Rvector = PROTECT(_new_Rvector0(Rtype, dim0));
	leaf = PROTECT(
		precompute_leaf_from_short_Rvector(
					dim0, index0, short_Rvector,
					&left_bufs)
	);
	left_bufs.precomputed_leaf = leaf;
	UNPROTECT(2);
	return left_bufs;
}

/* 'index0' must be either R_NilValue or an integer vector.
   'short_Rvector' must have a length >= 1. */
static SEXP subassign_leaf_with_short_Rvector(SEXP leaf, int dim0,
		SEXP index0, SEXP short_Rvector,
		LeftBufs *left_bufs)
{
	if (left_bufs->full_replacement || leaf == R_NilValue)
		return left_bufs->precomputed_leaf;

	SEXP left_Rvector = left_bufs->Rvector;
	_expand_leaf(leaf, left_Rvector, 0);
	int short_len = LENGTH(short_Rvector);
	int d2 = LENGTH(index0);
	for (int i2 = 0; i2 < d2; i2++) {
		int coord = INTEGER(index0)[i2];
		if (INVALID_COORD(coord, dim0))
			error("subscript contains "
			      "out-of-bound indices or NAs");
		int i1 = coord - 1;
		/* Virtual recycling of 'short_Rvector'. */
		left_bufs->copy_Rvector_elt_FUN(
				short_Rvector, i2 % short_len,
				left_Rvector, i1);
	}
	SEXP ans = PROTECT(_make_leaf_from_Rsubvec(left_Rvector, 0, dim0,
						   left_bufs->offs, 0));
	if (ans != R_NilValue) {
		/* Remove nonzeros introduced in 'left_bufs->Rvector'. */
		SEXP ans_nzoffs = get_leaf_nzoffs(ans);
		_set_selected_Rsubvec_elts_to_zero(left_Rvector, 0,
					     INTEGER(ans_nzoffs),
					     LENGTH(ans_nzoffs));
	}
	UNPROTECT(1);
	return ans;
}

/* Recursive. 'ndim' must be >= 2. */
static SEXP REC_subassign_SVT_with_short_Rvector(SEXP SVT, SEXP SVT0,
		const int *dim, int ndim, SEXP Nindex,
		SEXP short_Rvector, LeftBufs *left_bufs)
{
	SEXP subSVT0 = R_NilValue;
	int d1 = dim[ndim - 1];
	SEXP Nindex_elt = VECTOR_ELT(Nindex, ndim - 1);
	int d2 = Nindex_elt == R_NilValue ? d1 : LENGTH(Nindex_elt);
	//printf("ndim = %d: d2 = %d\n", ndim, d2);
	for (int i2 = 0; i2 < d2; i2++) {
		int i1;
		if (Nindex_elt == R_NilValue) {
			i1 = i2;
		} else {
			int coord = INTEGER(Nindex_elt)[i2];
			if (INVALID_COORD(coord, d1))
				error("subscript contains "
				      "out-of-bound indices or NAs");
			i1 = coord - 1;
		}
		//printf("ndim = %d: i1 = %d i2 = %d\n", ndim, i1, i2);
		SEXP subSVT = VECTOR_ELT(SVT, i1);
		if (ndim == 2) {
			subSVT = PROTECT(
				subassign_leaf_with_short_Rvector(
					subSVT, dim[0],
					VECTOR_ELT(Nindex, 0), short_Rvector,
					left_bufs)
			);
		} else {
			if (SVT0 != R_NilValue)
				subSVT0 = VECTOR_ELT(SVT0, i1);
			subSVT = PROTECT(
				make_SVT_node(subSVT, dim[ndim - 2], subSVT0)
			);
			subSVT = PROTECT(
				REC_subassign_SVT_with_short_Rvector(
					subSVT, subSVT0,
					dim, ndim - 1, Nindex,
					short_Rvector, left_bufs)
			);
		}
		SET_VECTOR_ELT(SVT, i1, subSVT);
		UNPROTECT(ndim == 2 ? 1 : 2);
	}
	int is_empty = 1;
	for (int i1 = 0; i1 < d1; i1++) {
		if (VECTOR_ELT(SVT, i1) != R_NilValue) {
			is_empty = 0;
			break;
		}
	}
	return is_empty ? R_NilValue : SVT;
}

/* --- .Call ENTRY POINT ---
   'Nindex' must be an N-index, that is, a list of integer vectors (or NULLs),
   one along each dimension in the array. */
SEXP C_subassign_SVT_with_short_Rvector(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		SEXP Nindex, SEXP Rvector)
{
	SEXPTYPE Rtype = _get_and_check_Rtype_from_Rstring(x_type,
				"C_subassign_SVT_with_short_Rvector", "x_type");
	if (TYPEOF(Rvector) != Rtype)
		error("SparseArray internal error in "
		      "C_subassign_SVT_with_short_Rvector():\n"
		      "    SVT_SparseArray object and 'Rvector' "
		      "must have the same type");

	const int *dim = INTEGER(x_dim);
	int ndim = LENGTH(x_dim);
	for (int along = 0; along < ndim; along++)
		if (dim[along] == 0)
			return x_SVT;  /* no-op */

	int dim0 = dim[0];
	SEXP index0 = VECTOR_ELT(Nindex, 0);

	LeftBufs left_bufs = init_left_bufs(dim0, index0, Rvector);
	PROTECT(left_bufs.Rvector);
	PROTECT(left_bufs.precomputed_leaf);

	if (ndim == 1) {
		SEXP ans = subassign_leaf_with_short_Rvector(
					x_SVT, dim0,
					index0, Rvector, &left_bufs);
		UNPROTECT(2);
		return ans;
	}

	SEXP ans = PROTECT(make_SVT_node(x_SVT, dim[ndim - 1], x_SVT));
	ans = REC_subassign_SVT_with_short_Rvector(ans, x_SVT,
					dim, ndim, Nindex,
					Rvector, &left_bufs);
	UNPROTECT(3);
	return ans;
}


/****************************************************************************
 * C_subassign_SVT_with_Rarray() and C_subassign_SVT_with_SVT()
 */

/* --- .Call ENTRY POINT ---
   The left and right arrays ('x' and 'Rarray') must have the same number
   of dimensions.
   'Nindex' must be an N-index, that is, a list of integer vectors (or NULLs),
   one along each dimension in the arrays. */
SEXP C_subassign_SVT_with_Rarray(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		SEXP Nindex, SEXP Rarray)
{
	error("not ready yet");
	return R_NilValue;
}

/* --- .Call ENTRY POINT ---
   The left and right arrays ('x' and 'v') must have the same number
   of dimensions.
  'Nindex' must be an N-index, that is, a list of integer vectors (or NULLs),
   one along each dimension in the arrays. */
SEXP C_subassign_SVT_with_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		SEXP Nindex, SEXP v_dim, SEXP v_type, SEXP v_SVT)
{
	error("not ready yet");
	return R_NilValue;
}


/****************************************************************************
 * ABANDONNED CODE
 *
 * The code below was an early attempt at solving the
 * C_subassign_SVT_with_short_Rvector() problem with a non-recursive
 * implementation. When I realized it was not going to work, I switched
 * to the REC_subassign_SVT_with_short_Rvector() solution (which is
 * recursive). I'm keeping the code below for now because the NindexIterator
 * thing works great (even though SVT_SparseArray subassignment is not a
 * good use case for it) and I might need it at some point for other things.
 */

#if 0

typedef struct Nindex_iterator_t {
	int ndim;
	const int *dim;
	SEXP Nindex;
	int margin;
	int *selection_dim;       /* of length 'ndim - margin' */
	int *selection_midx_buf;  /* of length 'ndim - margin' */
	long long selection_len;
	long long counter;
	int *coords0_buf;         /* of length 'ndim - margin' */
} NindexIterator;

static long long init_NindexIterator(NindexIterator *Nindex_iter,
		const int *dim, int ndim, SEXP Nindex, int margin)
{
	long long selection_len;
	int along, sd;
	SEXP Nindex_elt;

	if (!isVectorList(Nindex) || LENGTH(Nindex) != ndim)
		error("incorrect number of subscripts");
	Nindex_iter->ndim = ndim;
	Nindex_iter->dim = dim;
	Nindex_iter->Nindex = Nindex;
	Nindex_iter->margin = margin;
	Nindex_iter->selection_dim =
		(int *) R_alloc(ndim - margin, sizeof(int));
	Nindex_iter->selection_midx_buf =
		(int *) R_alloc(ndim - margin, sizeof(int));
	selection_len = 1;
	for (along = 0; along < ndim; along++) {
		Nindex_elt = VECTOR_ELT(Nindex, along);
		if (Nindex_elt == R_NilValue) {
			sd = dim[along];
		} else if (IS_INTEGER(Nindex_elt)) {
			sd = LENGTH(Nindex_elt);
		} else {
			error("subscripts must be integer vectors");
		}
		selection_len *= sd;
		if (along < margin)
			continue;
		Nindex_iter->selection_dim[along - margin] = sd;
		Nindex_iter->selection_midx_buf[along - margin] = 0;
	}
	Nindex_iter->selection_len = selection_len;
	Nindex_iter->counter = -1;
	Nindex_iter->coords0_buf =
		(int *) R_alloc(ndim - margin, sizeof(int));
	return selection_len;
}

static inline int next_midx(int ndim, const int *max_idx_plus_one,
			    int *midx_buf)
{
	int along, i;

	for (along = 0; along < ndim; along++) {
		i = midx_buf[along] + 1;
		if (i < max_idx_plus_one[along]) {
			midx_buf[along] = i;
			break;
		}
		midx_buf[along] = 0;
        }
        return along;
}

/* Returns:
       1 = if the array coords before the move was not the last one in the
	   array selection and the move to the next one was successful;
       0 = if the array coords before the move was the last one in the
	   array selection and so the move to the next one was not possible;
     < 0 = if error
   Typical use:
       while (ret = next_coords0(&Nindex_iter)) {
           if (ret < 0) {
               an error occured
           }
           handle current array coords
       }
 */
static inline int next_coords0(NindexIterator *Nindex_iter)
{
	int moved_along, along, *coords0_p, coord;
	const int *midx_p;
	SEXP Nindex_elt;

	if (Nindex_iter->selection_len == 0)
		return 0;
	if (Nindex_iter->counter == -1) {
		moved_along = Nindex_iter->ndim;
	} else {
		/* Update 'Nindex_iter->selection_midx_buf'. */
		moved_along = Nindex_iter->margin +
			next_midx(Nindex_iter->ndim - Nindex_iter->margin,
				  Nindex_iter->selection_dim,
				  Nindex_iter->selection_midx_buf);
		if (moved_along == Nindex_iter->ndim)
			return 0;
	}
	Nindex_iter->counter++;
	//printf("Nindex_iter->counter=%lld\n", Nindex_iter->counter);
	//printf("moved_along=%d\n", moved_along);

	/* Update 'Nindex_iter->coords0_buf'. */
	midx_p = Nindex_iter->selection_midx_buf;
	coords0_p = Nindex_iter->coords0_buf;
	for (along = Nindex_iter->margin; along < Nindex_iter->ndim; along++) {
		if (along > moved_along)
                        break;
		Nindex_elt = VECTOR_ELT(Nindex_iter->Nindex, along);
		if (Nindex_elt == R_NilValue) {
			*coords0_p = *midx_p;
		} else {
			coord = INTEGER(Nindex_elt)[*midx_p];
			if (INVALID_COORD(coord, Nindex_iter->dim[along]))
				error("subscript contains "
				      "out-of-bound indices or NAs");
			*coords0_p = coord - 1;
		}
		midx_p++;
		coords0_p++;
	}
	printf("coords0: ");
	coords0_p = Nindex_iter->coords0_buf;
	for (along = Nindex_iter->margin; along < Nindex_iter->ndim; along++) {
		printf(" %3d", *coords0_p);
		coords0_p++;
	}
	printf("\n");
	return 1;
}

/* 'Nindex_iter' declared and initialized with init_NindexIterator()
   in the caller (must be called with 'margin' set to 1). */
static int subassign_SVT_with_leaf(SEXP SVT, SEXP SVT0,
		NindexIterator *Nindex_iter, SEXP Rleaf)
{
	int ret, i;
	SEXP leaf_parent, leaf;

	while ((ret = next_coords0(Nindex_iter))) {
		if (ret < 0) {
			error("an error occured");
		}
		ret = descend_to_bottom_by_coords0(SVT, SVT0,
				Nindex_iter->dim, Nindex_iter->ndim,
				Nindex_iter->coords0_buf,
				&leaf_parent, &i, &leaf);
		if (ret < 0)
			return -1;
	}
	return 0;
}

#endif

