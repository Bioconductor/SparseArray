/****************************************************************************
 ****************************************************************************
 **									   **
 **        Single-bracket subsetting (`[`) of a SparseArray object         **
 **									   **
 ****************************************************************************
 ****************************************************************************/
#include "SparseArray_subsetting.h"

#include "argcheck_utils.h"
#include "OPBufTree.h"
#include "thread_control.h"  /* for which_max() */
#include "leaf_utils.h"

#include <limits.h>  /* for INT_MAX */
#include <string.h>  /* for memcpy() */
//#include <time.h>


/* There's no NA value for types "raw" or "list" so in this case we set to
   zero, that is, to Rbyte0 for "raw" and to NULL for "list". This mimics
   what 'as.raw(11:15)[c(3L, NA)]' and 'as.list(11:15)[c(3L, NA)]' do.
   TODO (maybe): Move this to Rvector_utils.c. */
static inline void set_Rvector_elt_to_NA(SEXP Rvector, R_xlen_t i)
{
	SEXPTYPE Rtype = TYPEOF(Rvector);
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		INTEGER(Rvector)[i] = NA_INTEGER;
		return;
	    case REALSXP:
		REAL(Rvector)[i] = NA_REAL;
		return;
	    case CPLXSXP: {
		Rcomplex *z = COMPLEX(Rvector) + i;
		z->r = z->i = NA_REAL;
		return;
	    }
	    case RAWSXP:
		RAW(Rvector)[i] = Rbyte0;
		return;
	    case STRSXP:
		SET_STRING_ELT(Rvector, i, NA_STRING);
		return;
	    case VECSXP:
		SET_VECTOR_ELT(Rvector, i, R_NilValue);
		return;
	}
	error("SparseArray internal error in "
	      "set_Rvector_elt_to_NA():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
}


/****************************************************************************
 * subset_NULL_by_Lindex()
 * subset_leaf_by_Lindex()
 * subset_leaf_by_OPBuf()
 * subset_leaf_as_sparse()
 */

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

/* Returns a value >= 0 and < 'nzcount' if success, or -1 if failure. */
static inline int bsearch_idx0_to_k2(int idx0, const int *nzoffs_p, int nzcount)
{
	/* Compare with first offset. */
	int k1 = 0;
	int nzoff = nzoffs_p[k1];
	if (idx0 < nzoff)
		return -1;
	if (idx0 == nzoff)
		return k1;

	/* Compare with last offset. */
	int k2 = nzcount - 1;
	nzoff = nzoffs_p[k2];
	if (idx0 > nzoff)
		return -1;
	if (idx0 == nzoff)
		return k2;

	/* Binary search.
	   Seems that using >> 1 instead of / 2 is faster, even when compiling
	   with 'gcc -O2' (one would hope that the optimizer is able to do that
	   kind of optimization). */
	int k;
	while ((k = (k1 + k2) >> 1) != k1) {
		nzoff = nzoffs_p[k];
		if (idx0 == nzoff)
			return k;
		if (idx0 > nzoff)
			k1 = k;
		else
			k2 = k;
	}
	return -1;
}

/* Mapping 'idx0' to k thru the lookup table is always faster than using a
   binary search. However building the lookup table has a small cost so is
   only worth it if we're going to use it to map more than one 'idx0' value.
   TODO: Right now we use a cutoff value of 10 but this needs to be refined.
   The overhead of building the lookup table depends on the value of 'nzcount'
   so the cutoff value should be a function of 'nzcount'. */
#define	MAP_IDX0_TO_K(idx0, nzoffs_p, nzcount) \
	use_lookup_table ? lookup_table[(idx0)] \
			 : bsearch_idx0_to_k2((idx0), (nzoffs_p), (nzcount))

static int subset_NULL_by_Lindex(int dim0, SEXP Lindex, SEXP ans)
{
	int n = LENGTH(Lindex);

	/* We only care about NAs or NaNs in 'Lindex'. */
	for (int k1 = 0; k1 < n; k1++) {
		int idx0;
		int ret = extract_idx0(Lindex, k1, dim0, &idx0);
		if (ret == SUBSCRIPT_ELT_IS_NA) {
			/* 'Lindex[k1]' is NA or NaN. */
			set_Rvector_elt_to_NA(ans, (R_xlen_t) k1);
			continue;
		}
		if (ret < 0)
			return ret;
	}
	return 0;
}

/* Used for linear subsetting of a 1D SVT_SparseArray object.
   'Lindex' must be a numeric vector. It cannot be a long one! It is expected
   to contain 1-based indices that are >= 1 and <= 'dim0' (in the case of a
   1D SVT_SparseArray object, 'dim0' is also the length of the object).
   NA indices are ok. */
static int subset_leaf_by_Lindex(SEXP leaf, int dim0, SEXP Lindex, SEXP ans,
		CopyRVectorEltFUN copy_Rvector_elt_FUN)
{
	if (leaf == R_NilValue)
		return subset_NULL_by_Lindex(dim0, Lindex, ans);

	int n = LENGTH(Lindex);
	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
	const int *nzoffs_p = INTEGER(nzoffs);

	/* Walk on 'Lindex'. */
	for (int k1 = 0; k1 < n; k1++) {
		int idx0;
		int ret = extract_idx0(Lindex, k1, dim0, &idx0);
		if (ret == SUBSCRIPT_ELT_IS_NA) {
			/* 'Lindex[k1]' is NA or NaN. */
			set_Rvector_elt_to_NA(ans, (R_xlen_t) k1);
			continue;
		}
		if (ret < 0)
			return ret;
		int k2 = bsearch_idx0_to_k2(idx0, nzoffs_p, nzcount);
		if (k2 >= 0)
			copy_Rvector_elt_FUN(nzvals, (R_xlen_t) k2,
					     ans, (R_xlen_t) k1);
	}
	return 0;
}

static void subset_leaf_by_OPBuf(SEXP leaf, const OPBuf *opbuf, SEXP ans,
		CopyRVectorEltFUN copy_Rvector_elt_FUN, int *lookup_table)
{
	if (leaf == R_NilValue || opbuf->nelt == 0)
		return;

	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
	const int *nzoffs_p = INTEGER(nzoffs);
	const int *idx0s = opbuf->idx0s;
	const int *Loffs = opbuf->Loffs;
	const R_xlen_t *xLoffs = opbuf->xLoffs;

	/* See comment preceding MAP_IDX0_TO_K() definition above. */
	int use_lookup_table = opbuf->nelt > 10;
	if (use_lookup_table)
		build_lookup_table(lookup_table, nzoffs_p, nzcount);
	/* Walk on the (idx0,Loff) pairs. */
	if (Loffs != NULL) {
		for (int k1 = 0; k1 < opbuf->nelt; k1++) {
			int idx0 = idx0s[k1];
			int k2 = MAP_IDX0_TO_K(idx0, nzoffs_p, nzcount);
			if (k2 >= 0)
				copy_Rvector_elt_FUN(nzvals, (R_xlen_t) k2,
						     ans, (R_xlen_t) Loffs[k1]);
		}
	} else {
		for (int k1 = 0; k1 < opbuf->nelt; k1++) {
			int idx0 = idx0s[k1];
			int k2 = MAP_IDX0_TO_K(idx0, nzoffs_p, nzcount);
			if (k2 >= 0)
				copy_Rvector_elt_FUN(nzvals, (R_xlen_t) k2,
						     ans, xLoffs[k1]);
		}
	}
	if (use_lookup_table)
		reset_lookup_table(lookup_table, nzoffs_p, nzcount);
	return;
}

/* 'subscript' must be a numeric vector. It cannot be a long one! It is
   expected to contain 1-based indices that are >= 1 and <= 'sv->len'.
   NA indices will trigger an error.
   'sv_selection' and 'out_nzoffs' must be arrays that are long enough
   to hold at least 'LENGTH(subscript)' ints. */
static int subset_SV(const SparseVec *sv, SEXP subscript,
		     int *sv_selection, int *out_nzoffs, int *lookup_table)
{
	int out_nzcount = 0;
	int n = LENGTH(subscript);
	if (n == 0)
		return out_nzcount;

	int sv_nzcount = get_SV_nzcount(sv);

	/* See comment preceding MAP_IDX0_TO_K() definition above. */
	int use_lookup_table = n > 10;
	if (use_lookup_table)
		build_lookup_table(lookup_table, sv->nzoffs, sv_nzcount);
	/* Walk on 'subscript'. */
	int idx0 = 0;  // only for -Wmaybe-uninitialized
	for (int i1 = 0; i1 < n; i1++) {
		int ret = extract_idx0(subscript, i1, sv->len, &idx0);
		if (ret < 0)
			_bad_Nindex_error(ret, 1);
		int k2 = MAP_IDX0_TO_K(idx0, sv->nzoffs, sv_nzcount);
		if (k2 >= 0) {
			sv_selection[out_nzcount] = k2;
			out_nzoffs[out_nzcount] = i1;
			out_nzcount++;
		}
	}
	if (use_lookup_table)
		reset_lookup_table(lookup_table, sv->nzoffs, sv_nzcount);
	return out_nzcount;
}

/* Takes a non-NULL leaf (standard or lacunar), and returns a leaf
   that can be NULL, standard, or lacunar.
   'subscript' must be NULL or a numeric vector. It cannot be a long one!
   It is expected to contain 1-based indices that are >= 1 and <= 'dim0'.
   NA indices will trigger an error. */
static SEXP subset_leaf_as_sparse(SEXP leaf, int dim0, SEXP subscript,
		int *selection_buf, int *nzoffs_buf, int *lookup_table)
{
	if (subscript == R_NilValue)
		return leaf;

	SEXP leaf_nzvals = get_leaf_nzvals(leaf);
	SparseVec sv = leaf2SV(leaf, TYPEOF(leaf_nzvals), dim0, 0);
	int ans_nzcount = subset_SV(&sv, subscript,
				    selection_buf, nzoffs_buf, lookup_table);
	if (ans_nzcount == 0)
		return R_NilValue;

	SEXP ans_nzoffs = PROTECT(NEW_INTEGER(ans_nzcount));
	memcpy(INTEGER(ans_nzoffs), nzoffs_buf, sizeof(int) * ans_nzcount);
	if (leaf_nzvals == R_NilValue) {  /* input leaf is lacunar */
		SEXP ans = _make_lacunar_leaf(ans_nzoffs);
		UNPROTECT(1);
		return ans;
	}
	/* input leaf is standard */
	if (LACUNAR_MODE_IS_ON) {
		int all_ones = _all_selected_Rsubvec_elts_equal_one(
						leaf_nzvals, 0,
						selection_buf, ans_nzcount);
		if (all_ones) {
			SEXP ans = _make_lacunar_leaf(ans_nzoffs);
			UNPROTECT(1);
			return ans;
		}
	}
	SEXP ans_nzvals = PROTECT(
		_subset_Rsubvec(leaf_nzvals, 0, selection_buf, ans_nzcount)
	);
	SEXP ans = zip_leaf(ans_nzvals, ans_nzoffs, 0);
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * build_OPBufTree_from_Lindex()
 * build_OPBufTree_from_Mindex()
 *
 * Both are used by .Call ENTRY POINTs C_subset_SVT_by_[L|M]index() below
 * in this file.
 */

/* 'Lidx0' is trusted to be a non-NA value >= 0 and < 'dimcumprod[ndim - 1]'.
   'SVT' is trusted to not be NULL.
   Returns NULL if we didn't land anywhere. */
static OPBufTree *find_host_node_for_Lidx0(OPBufTree *opbuf_tree,
		R_xlen_t Lidx0,
		SEXP SVT, const int *dim, int ndim,
		const R_xlen_t *dimcumprod, int *idx0)
{
	for (int along = ndim - 1; along >= 1; along--) {
		R_xlen_t p = dimcumprod[along - 1];
		int i = Lidx0 / p;  /* always >= 0 and < 'dim[along]' */
		SVT = VECTOR_ELT(SVT, i);
		if (SVT == R_NilValue)
			return NULL;
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

/* 'SVT' is trusted to not be NULL.
   Returns NULL if we didn't land anywhere. */
static OPBufTree *find_host_node_for_Mindex_row(OPBufTree *opbuf_tree,
		SEXP Mindex, R_xlen_t Moff, int out_len,
		SEXP SVT, const int *dim, int ndim,
		int *idx0, int *ret_code)
{
	for (int along = ndim - 1; along >= 1; along--, Moff -= out_len) {
		int i, d = dim[along];  /* = LENGTH(SVT) */
		*ret_code = extract_idx0(Mindex, Moff, d, &i);
		if (*ret_code < 0)
			return NULL;
		SVT = VECTOR_ELT(SVT, i);
		if (SVT == R_NilValue)
			return NULL;
		if (opbuf_tree->node_type == NULL_NODE)
			_alloc_OPBufTree_children(opbuf_tree, d);
		opbuf_tree = get_OPBufTree_child(opbuf_tree, i);
	}
	/* At this point 'opbuf_tree' is guaranteed to be a node of type
	   NULL_NODE or LEAF_NODE. */
	*ret_code = extract_idx0(Mindex, Moff, dim[0], idx0);
	return opbuf_tree;
}

/* To use on an 'Lindex' that has a length <= INT_MAX.
   Based on **unsafe** _append_idx0Loff_to_host_node().
   See comments for _append_idx0Loff_to_host_node() and
   append_idx0Loff_to_OPBuf() in OPBufTree.c for the details.
   Returns a negative value in case of error. */
static int build_OPBufTree_from_Lindex1(OPBufTree *opbuf_tree, SEXP Lindex,
		SEXP x_SVT, const int *x_dim, int x_ndim, SEXP ans,
		const R_xlen_t *dimcumprod)
{
	int max_outleaf_len = 0;
	int out_len = LENGTH(Lindex);  /* = LENGTH(ans) */
	R_xlen_t x_len = dimcumprod[x_ndim - 1];
	/* Walk along 'Lindex' (and 'ans').
	   Direction of the walk doesn't matter. */
	for (int Loff = 0; Loff < out_len; Loff++) {
		R_xlen_t Lidx0;
		int ret = extract_long_idx0(Lindex, (R_xlen_t) Loff, x_len,
					    &Lidx0);
		if (ret == SUBSCRIPT_ELT_IS_NA) {
			/* 'Lindex[Loff]' is NA or NaN. */
			set_Rvector_elt_to_NA(ans, (R_xlen_t) Loff);
			continue;
		}
		if (ret < 0)
			return ret;
		if (x_SVT == R_NilValue)
			continue;
		int idx0;
		OPBufTree *host_node = find_host_node_for_Lidx0(
						opbuf_tree, Lidx0,
						x_SVT, x_dim, x_ndim,
						dimcumprod, &idx0);
		if (host_node == NULL)  /* we didn't land anywhere */
			continue;
		/* Turn 'host_node' into a leaf node (LEAF_NODE) if it is
		   not one already, and append the (idx0,Loff) pair to it.
		   _append_idx0Loff_to_host_node() will return the new
		   length of the leaf node, that is, the new nb of pairs
		   in it (this will always be > 0), or a negative value in
		   case of error. */
		ret = _append_idx0Loff_to_host_node(host_node, idx0, Loff);
		if (ret < 0)
			return ret;
		if (ret > max_outleaf_len)
			max_outleaf_len = ret;
	}
	return max_outleaf_len;
}

/* To use on an 'Lindex' that has a length > INT_MAX (long vector).
   Returns a negative value in case of error. */
static int build_OPBufTree_from_Lindex2(OPBufTree *opbuf_tree, SEXP Lindex,
		SEXP x_SVT, const int *x_dim, int x_ndim, SEXP ans,
		const R_xlen_t *dimcumprod)
{
	int max_outleaf_len = 0;
	R_xlen_t out_len = XLENGTH(Lindex);  /* = XLENGTH(ans) */
	R_xlen_t x_len = dimcumprod[x_ndim - 1];
	/* Walk along 'Lindex' (and 'ans').
	   Direction of the walk doesn't matter. */
	//for (R_xlen_t Loff = 0; Loff < out_len; Loff++) {
	for (R_xlen_t Loff = out_len - 1; Loff >= 0; Loff--) {
		R_xlen_t Lidx0;
		int ret = extract_long_idx0(Lindex, Loff, x_len, &Lidx0);
		if (ret == SUBSCRIPT_ELT_IS_NA) {
			/* 'Lindex[Loff]' is NA or NaN. */
			set_Rvector_elt_to_NA(ans, Loff);
			continue;
		}
		if (ret < 0)
			return ret;
		if (x_SVT == R_NilValue)
			continue;
		int idx0;
		OPBufTree *host_node = find_host_node_for_Lidx0(
						opbuf_tree, Lidx0,
						x_SVT, x_dim, x_ndim,
						dimcumprod, &idx0);
		if (host_node == NULL)  /* we didn't land anywhere */
			continue;
		/* Turn 'host_node' into a leaf node (LEAF_NODE) if it is
		   not one already, and append the (idx0,Loff) pair to it.
		   _append_idx0Loff_to_host_node() will return the new
		   length of the leaf node, that is, the new nb of pairs
		   in it (this will always be > 0), or a negative value in
		   case of error.
		   Note that the only possible error code at the moment is
		   MAX_OPBUF_LEN_REACHED (defined in OPBufTree.h), which
		   indicates that 'host_node' is already at its max length
		   (INT_MAX). This can only happen if 'Lindex' is a long
		   vector and more than INT_MAX indices in it hit the same
		   leaf in 'SVT'. A rather crazy and unlikely situation! */
		ret = _append_idx0xLoff_to_host_node(host_node, idx0, Loff);
		if (ret < 0)
			return ret;
		if (ret > max_outleaf_len)
			max_outleaf_len = ret;
	}
	return max_outleaf_len;
}

/* Returns the length of the longest leaf in the output tree (i.e. the
   resulting OPBufTree), or < 0 in case of an error. Therefore a returned
   value of 0 means that no leaves in the input tree ('x_SVT') are touched by
   the subsetting operation (i.e. no leaves in 'x_SVT' are "hit" by 'Lindex').

   IMPORTANT: The output tree will have the same morphology as the input tree
   except that not all nodes in the latter are necessarily mapped to a node
   in the former. In particular the output tree will only have leaves that
   correspond to leaves in the input tree that are touched by the subsetting
   operation.
   In other words, the output tree will be a pruned version of the input tree,
   unless all the leaves in the latter are touched by the subsetting operation,
   in which case the two trees will have identical morphologies. This means
   that all the nodes in the resulting OPBufTree will have a corresponding
   node in the input SVT. This is the reason why the recursive tree traversal
   in REC_subset_SVT_by_OPBufTree() below is guided by 'opbuf_tree' and not
   by 'SVT'. */
static int build_OPBufTree_from_Lindex(OPBufTree *opbuf_tree, SEXP Lindex,
		SEXP x_SVT, const int *x_dim, int x_ndim, SEXP ans,
		const R_xlen_t *dimcumprod)
{
	/* _free_OPBufTree(opbuf_tree) resets 'opbuf_tree->node_type'
	   to NULL_NODE. */
	_free_OPBufTree(opbuf_tree);
	return XLENGTH(Lindex) <= (R_xlen_t) INT_MAX ?
		build_OPBufTree_from_Lindex1(opbuf_tree, Lindex,
				x_SVT, x_dim, x_ndim, ans, dimcumprod) :
		build_OPBufTree_from_Lindex2(opbuf_tree, Lindex,
				x_SVT, x_dim, x_ndim, ans, dimcumprod);
}

static int build_OPBufTree_from_Mindex(OPBufTree *opbuf_tree, SEXP Mindex,
		SEXP x_SVT, const int *x_dim, int x_ndim, SEXP ans)
{
	/* _free_OPBufTree(opbuf_tree) resets 'opbuf_tree->node_type'
	   to NULL_NODE. */
	_free_OPBufTree(opbuf_tree);
	int max_outleaf_len = 0;
	int out_len = LENGTH(ans);  /* = nrow(Mindex) */
	R_xlen_t Moff = (R_xlen_t) out_len * (x_ndim - 1);
	/* Walk along 'ans'. Direction of the walk doesn't matter. */
	for (int Loff = 0; Loff < out_len; Loff++, Moff++) {
		int idx0, ret;
		OPBufTree *host_node = find_host_node_for_Mindex_row(
						opbuf_tree,
						Mindex, Moff, out_len,
						x_SVT, x_dim, x_ndim,
						&idx0, &ret);
		if (ret < 0)
			return ret;
		if (host_node == NULL)  /* we didn't land anywhere */
			continue;
		/* Turn 'host_node' into a leaf node (LEAF_NODE) if it is
		   not one already, and append the (idx0,Loff) pair to it.
		   _append_idx0Loff_to_host_node() will return the new
		   length of the leaf node, that is, the new nb of pairs
		   in it (this will always be > 0), or a negative value in
		   case of error. */
		ret = _append_idx0Loff_to_host_node(host_node, idx0, Loff);
		if (ret < 0)
			return ret;
		if (ret > max_outleaf_len)
			max_outleaf_len = ret;
	}
	return max_outleaf_len;
}


/****************************************************************************
 * C_subset_SVT_by_[L|M]index()
 */

static SEXP new_Rvector(SEXPTYPE Rtype, R_xlen_t len, int na_background)
{
	if (na_background)
		return _new_RvectorNA(Rtype, len);
	return _new_Rvector0(Rtype, len);
}

/* Recursive tree traversal must be guided by 'opbuf_tree', not by 'SVT'. See
   IMPORTANT comment for build_OPBufTree_from_Lindex() above. */
static void REC_subset_SVT_by_OPBufTree(OPBufTree *opbuf_tree,
		SEXP SVT, const int *dim, int ndim, SEXP ans,
		CopyRVectorEltFUN fun, int *lookup_table, int pardim)
{
	if (opbuf_tree->node_type == NULL_NODE)
		return;

	if (ndim == 1) {
		/* Both 'opbuf_tree' and 'SVT' are leaves. */
		OPBuf *opbuf = get_OPBufTree_leaf(opbuf_tree);
		subset_leaf_by_OPBuf(SVT, opbuf, ans, fun, lookup_table);
		_free_OPBufTree(opbuf_tree);
		return;
	}

	/* Both 'opbuf_tree' and 'SVT' are inner nodes. */
	int n = get_OPBufTree_nchildren(opbuf_tree);  /* same as dim[ndim - 1]
							 or LENGTH(SVT) */
	/* Parallel execution along the biggest dimension only.
	   WARNING: We need to delcare 'lookup_table' as private if we
	   want to do this otherwise various threads will concurrently write
	   to it and step on each others feet. */
	//#pragma omp parallel for schedule(static) if(ndim == pardim)
	for (int i = 0; i < n; i++) {
		OPBufTree *child = get_OPBufTree_child(opbuf_tree, i);
		SEXP subSVT = VECTOR_ELT(SVT, i);
		REC_subset_SVT_by_OPBufTree(child,
				subSVT, dim, ndim - 1, ans,
				fun, lookup_table, pardim);
	}
	_free_OPBufTree(opbuf_tree);
	return;
}

/* --- .Call ENTRY POINT ---
   'Lindex' must be a numeric vector (integer or double), possibly a long one.
   NA indices are accepted. */
SEXP C_subset_SVT_by_Lindex(
		SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP x_na_background,
		SEXP Lindex)
{
	SEXPTYPE Rtype = _get_and_check_Rtype_from_Rstring(x_type,
				"C_subset_SVT_by_Lindex", "x_type");
	CopyRVectorEltFUN fun = _select_copy_Rvector_elt_FUN(Rtype);
	int x_has_NAbg = _get_and_check_na_background(x_na_background,
				"C_subset_SVT_by_Lindex", "x_na_background");

	if (!(IS_INTEGER(Lindex) || IS_NUMERIC(Lindex)))
		error("'Lindex' must be an integer or numeric vector");

	int x_ndim = LENGTH(x_dim);
	int x_dim0 = INTEGER(x_dim)[0];
	R_xlen_t ans_len = XLENGTH(Lindex);
	SEXP ans = PROTECT(new_Rvector(Rtype, ans_len, x_has_NAbg));

	if (x_ndim == 1) {
		int ret = subset_leaf_by_Lindex(x_SVT, x_dim0, Lindex, ans,
						fun);
		UNPROTECT(1);
		if (ret < 0)
			_bad_Lindex_error(ret);
		return ans;
	}

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
				x_SVT, INTEGER(x_dim), x_ndim, ans,
				dimcumprod);
	if (max_outleaf_len < 0) {
		UNPROTECT(1);
		_bad_Lindex_error(max_outleaf_len);
	}
	if (x_SVT == R_NilValue) {
		UNPROTECT(1);
		return ans;
	}
	//double dt = (1.0 * clock() - t0) * 1000.0 / CLOCKS_PER_SEC;
	//printf("1st pass: %2.3f ms\n", dt);

	//printf("max_outleaf_len = %d\n", max_outleaf_len);

	/* 2nd pass: Subset SVT by OPBufTree. */
	if (max_outleaf_len > 0) {
		//clock_t t0 = clock();
		int *lookup_table = (int *) R_alloc(x_dim0, sizeof(int));
		for (int i = 0; i < x_dim0; i++)
			lookup_table[i] = -1;
		/* Get 1-based rank of biggest dimension (ignoring the 1st
		   dim). Parallel execution will be along that dimension. */
		int pardim = which_max(INTEGER(x_dim) + 1, x_ndim - 1) + 2;
		REC_subset_SVT_by_OPBufTree(opbuf_tree,
				x_SVT, INTEGER(x_dim), x_ndim, ans,
				fun, lookup_table, pardim);
		//double dt = (1.0 * clock() - t0) * 1000.0 / CLOCKS_PER_SEC;
		//printf("2nd pass: %2.3f ms\n", dt);
	}

	UNPROTECT(1);
	return ans;
}

static int check_Mindex(SEXP Mindex, int ndim,
		const char *what1, const char *what2)
{
	SEXP Mindex_dim = GET_DIM(Mindex);
	if (Mindex_dim == R_NilValue || LENGTH(Mindex_dim) != 2)
		error("'%s' must be a matrix", what1);
	if (!(IS_INTEGER(Mindex) || IS_NUMERIC(Mindex)))
		error("'%s' must be an integer matrix", what1);
	if (INTEGER(Mindex_dim)[1] != ndim)
		error("ncol(%s) != %s", what1, what2);
	return INTEGER(Mindex_dim)[0];
}

/* --- .Call ENTRY POINT ---
   'Mindex' must be a numeric matrix (integer or double) with one column per
   dimension in the array to subset. NAs in the matrix are forbidden at the
   moment (they'll trigger an error), except in the 1D case. */
SEXP C_subset_SVT_by_Mindex(
		SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP x_na_background,
		SEXP Mindex)
{
	SEXPTYPE Rtype = _get_and_check_Rtype_from_Rstring(x_type,
				"C_subset_SVT_by_Mindex", "x_type");
	CopyRVectorEltFUN fun = _select_copy_Rvector_elt_FUN(Rtype);
	int x_has_NAbg = _get_and_check_na_background(x_na_background,
				"C_subset_SVT_by_Mindex", "x_na_background");

	int x_ndim = LENGTH(x_dim);
	int x_dim0 = INTEGER(x_dim)[0];
	int ans_len = check_Mindex(Mindex, x_ndim, "Mindex", "length(dim(x))");
	SEXP ans = PROTECT(new_Rvector(Rtype, (R_xlen_t) ans_len, x_has_NAbg));

	if (x_ndim == 1) {
		int ret = subset_leaf_by_Lindex(x_SVT, x_dim0, Mindex, ans,
						fun);
		UNPROTECT(1);
		if (ret < 0)
			_bad_Mindex_error(ret);
		return ans;
	}
	if (x_SVT == R_NilValue) {
		UNPROTECT(1);
		return ans;
	}

	/* 1st pass: Build OPBufTree. */
	OPBufTree *opbuf_tree = _get_global_opbuf_tree();
	int max_outleaf_len =
		build_OPBufTree_from_Mindex(opbuf_tree, Mindex,
					    x_SVT, INTEGER(x_dim), x_ndim, ans);
	if (max_outleaf_len < 0) {
		UNPROTECT(1);
		_bad_Mindex_error(max_outleaf_len);
	}

	/* 2nd pass: Subset SVT by OPBufTree. */
	if (max_outleaf_len > 0) {
		int *lookup_table = (int *) R_alloc(x_dim0, sizeof(int));
		for (int i = 0; i < x_dim0; i++)
			lookup_table[i] = -1;
		/* Get 1-based rank of biggest dimension (ignoring the 1st
		   dim). Parallel execution will be along that dimension. */
		int pardim = which_max(INTEGER(x_dim) + 1, x_ndim - 1) + 2;
		REC_subset_SVT_by_OPBufTree(opbuf_tree,
				x_SVT, INTEGER(x_dim), x_ndim, ans,
				fun, lookup_table, pardim);
	}

	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_subset_SVT_by_Nindex()
 */

static SEXP compute_subset_dim(SEXP Nindex, SEXP x_dim)
{
	int ndim = LENGTH(x_dim);
	if (!isVectorList(Nindex) || LENGTH(Nindex) != ndim)
		error("'Nindex' must be a list with one list "
		      "element along each dimension in 'x'");

	SEXP ans_dim = PROTECT(duplicate(x_dim));
	for (int along = 0; along < ndim; along++) {
		SEXP subscript = VECTOR_ELT(Nindex, along);
		if (subscript == R_NilValue)
			continue;
		if (!(IS_INTEGER(subscript) || IS_NUMERIC(subscript))) {
			UNPROTECT(1);
			_bad_Nindex_error(BAD_SUBSCRIPT_TYPE, along + 1);
		}
		R_xlen_t d = XLENGTH(subscript);
		if (d > (R_xlen_t) INT_MAX) {
			UNPROTECT(1);
			_bad_Nindex_error(SUBSCRIPT_IS_TOO_LONG, along + 1);
		}
		INTEGER(ans_dim)[along] = (int) d;
	}
	UNPROTECT(1);
	return ans_dim;
}

/* Recursive tree traversal.
   Returns R_NilValue or a list of length 'ans_dim[ndim - 1]'. */
static SEXP REC_subset_SVT_by_Nindex(SEXP SVT, SEXP Nindex,
		const int *x_dim, const int *ans_dim, int ndim,
		int *selection_buf, int *nzoffs_buf, int *lookup_table)
{
	if (SVT == R_NilValue)
		return R_NilValue;

	/* compute_subset_dim() already checked that 'subscript' is either
	   NULL or a numeric vector that is not a long vector.
	   If not NULL, 'subscript' is expected to contain 1-based indices
	   that are >= 1 and <= 'x_dim[ndim - 1]'. NA indices will trigger
	   an error. */
	SEXP subscript = VECTOR_ELT(Nindex, ndim - 1);

	if (ndim == 1) {
		/* 'SVT' is a leaf (i.e. 1D SVT). */
		return subset_leaf_as_sparse(SVT, x_dim[0], subscript,
				selection_buf, nzoffs_buf, lookup_table);
	}

	/* 'SVT' is a regular node (list). */
	int SVT_len = LENGTH(SVT);        /* same as 'x_dim[ndim - 1]' */
	int ans_len = ans_dim[ndim - 1];  /* same as 'LENGTH(subscript)'
					     if 'subscript' is not NULL */
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	int is_empty = 1;
	for (int i = 0; i < ans_len; i++) {
		int idx0 = i;
		if (subscript != R_NilValue) {
			int ret = extract_idx0(subscript, i, SVT_len, &idx0);
			if (ret < 0)
				_bad_Nindex_error(ret, ndim);
		}
		SEXP subSVT = VECTOR_ELT(SVT, idx0);
		SEXP ans_elt = REC_subset_SVT_by_Nindex(subSVT, Nindex,
					 x_dim, ans_dim, ndim - 1,
					 selection_buf, nzoffs_buf,
					 lookup_table);
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
   'Nindex' must be an N-index, that is, a list of numeric vectors (or NULLs),
   one along each dimension in the array to subset. Note that, strictly
   speaking, the vectors in an N-index are expected to be integer vectors,
   but C_subset_SVT_by_Nindex() can handle subscripts of type "double".
   NAs in the subscripts are forbidden (they'll trigger an error).  */
SEXP C_subset_SVT_by_Nindex(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP Nindex)
{
	/* Returned value ignored. */
	_get_and_check_Rtype_from_Rstring(x_type,
					  "C_subset_SVT_by_Nindex", "x_type");

	SEXP ans_dim = PROTECT(compute_subset_dim(Nindex, x_dim));
	int ans_dim0 = INTEGER(ans_dim)[0];
	int *selection_buf = (int *) R_alloc(ans_dim0, sizeof(int));
	int *nzoffs_buf = (int *) R_alloc(ans_dim0, sizeof(int));
	int x_dim0 = INTEGER(x_dim)[0];
	int *lookup_table = (int *) R_alloc(x_dim0, sizeof(int));
	for (int i = 0; i < x_dim0; i++)
		lookup_table[i] = -1;
	SEXP ans_SVT = REC_subset_SVT_by_Nindex(x_SVT, Nindex,
				 INTEGER(x_dim),
				 INTEGER(ans_dim), LENGTH(ans_dim),
				 selection_buf, nzoffs_buf, lookup_table);
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

