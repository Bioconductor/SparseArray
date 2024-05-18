/****************************************************************************
 ****************************************************************************
 **									   **
 **        Single-bracket subsetting (`[`) of a SparseArray object         **
 **									   **
 ****************************************************************************
 ****************************************************************************/
#include "SparseArray_subsetting.h"

#include "OPBufTree.h"
#include "thread_control.h"  /* for which_max() */
#include "Rvector_utils.h"
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
 * extract_long_idx0() and extract_idx0()
 */

/* In addition to 0, extract_long_idx0() and extract_idx0() can return one
   of the five values defined below. Negative values are always considered
   errors. Note that this could change for SUBSCRIPT_ELT_IS_BEYOND_MAX if we
   decided to strictly mimic subsetting in R base, in which case we would need
   to set its value to a positive one. SUBSCRIPT_ELT_IS_NA is set to a positive
   value because this is sometimes considered an error (by subset_SV() and
   REC_subset_SVT_by_Nindex()) and sometimes not (by subset_leaf_by_Lindex()
   and build_OPBufTree_from_Lindex()).
   Be aware that MAX_OPBUF_LEN_REACHED is set to -1 (see OPBufTree.h) so
   do **not** use that value. */
#define	BAD_SUBSCRIPT_TYPE             -2
#define	SUBSCRIPT_IS_TOO_LONG          -3  /* returned by extract_idx0() only */
#define	SUBSCRIPT_ELT_IS_LESS_THAN_ONE -4
#define	SUBSCRIPT_ELT_IS_BEYOND_MAX    -5
#define	SUBSCRIPT_ELT_IS_NA             1

/* 'subscript' must be a numeric vector, possibly a long one. It is expected
   to contain 1-based indices that are >= 1 and <= 'max'.
   Returns 0 or one of the values defined above. Will set 'idx0' only if
   returning 0. In other words, caller **must** ignore 'idx0' when a non-zero
   value is returned. */
static inline int extract_long_idx0(SEXP subscript, R_xlen_t i,
				    R_xlen_t max, R_xlen_t *idx0)
{
	if (IS_INTEGER(subscript)) {
		int idx = INTEGER(subscript)[i];
		if (idx == NA_INTEGER)
			return SUBSCRIPT_ELT_IS_NA;
		if (idx < 1)
			return SUBSCRIPT_ELT_IS_LESS_THAN_ONE;
		if ((R_xlen_t) idx > max)
			return SUBSCRIPT_ELT_IS_BEYOND_MAX;
		idx--;	/* from 1-based to 0-based */
		*idx0 = (R_xlen_t) idx;
		return 0;
	}
	if (IS_NUMERIC(subscript)) {
		double idx = REAL(subscript)[i];
		/* ISNAN(): True for *both* NA and NaN. See <R_ext/Arith.h> */
		if (ISNAN(idx))
			return SUBSCRIPT_ELT_IS_NA;
		if (idx < 1.0)
			return SUBSCRIPT_ELT_IS_LESS_THAN_ONE;
		if (idx > (double) max)
			return SUBSCRIPT_ELT_IS_BEYOND_MAX;
		idx -= 1.0;  /* from 1-based to 0-based */
		*idx0 = (R_xlen_t) idx;
		return 0;
	}
	return BAD_SUBSCRIPT_TYPE;
}

/* Like extract_long_idx0(), except that 'subscript' cannot be a long vector
   nor can it contain values > INT_MAX ('max' must be supplied as an 'int'). */
static inline int extract_idx0(SEXP subscript, int i, int max, int *idx0)
{
	if (XLENGTH(subscript) > INT_MAX)
		return SUBSCRIPT_IS_TOO_LONG;
	R_xlen_t lidx0 = 0;
	int ret = extract_long_idx0(subscript, (R_xlen_t) i,
				    (R_xlen_t) max, &lidx0);
	if (ret != 0)
		return ret;
	*idx0 = (int) lidx0;
	return 0;
}

/* 'ret_code' must be a **non-zero** code returned by extract_idx0() above.
   All of them are considered errors regardless of their sign. */
static void bad_Nindex_error(int ret_code, int along1)
{
	if (ret_code == BAD_SUBSCRIPT_TYPE)
		error("'Nindex[[%d]]' is not a numeric vector (or a NULL)",
		      along1);
	if (ret_code == SUBSCRIPT_IS_TOO_LONG)
		error("'Nindex[[%d]]' is a long vector", along1);
	if (ret_code == SUBSCRIPT_ELT_IS_NA)
		error("'Nindex[[%d]]' contains NAs", along1);
	error("'Nindex[[%d]]' contains out-of-bound indices", along1);
}


/****************************************************************************
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
static inline int bsearch_idx0_to_k2(int idx0, const int *nzoffs, int nzcount)
{
	/* Compare with first offset. */
	int k1 = 0;
	int nzoff = nzoffs[k1];
	if (idx0 < nzoff)
		return -1;
	if (idx0 == nzoff)
		return k1;

	/* Compare with last offset. */
	int k2 = nzcount - 1;
	nzoff = nzoffs[k2];
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
		nzoff = nzoffs[k];
		if (idx0 == nzoff)
			return k;
		if (idx0 > nzoff)
			k1 = k;
		else
			k2 = k;
	}
	return -1;
}

static inline void copy_nzval_elt(SEXP nzvals, int k,
		SEXP out, R_xlen_t out_offset,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	if (nzvals == R_NilValue) {
		/* lacunar leaf */
		_set_Rsubvec_elts_to_one(out, out_offset, (R_xlen_t) 1);
	} else {
		/* standard leaf */
		copy_Rvector_elt_FUN(nzvals, (R_xlen_t) k, out, out_offset);
	}
	return;
}

/* 'Lindex' must be a numeric vector. It should not be a long one!
   It is expected to contain 1-based indices that are >= 1 and <= 'dim0'.
   NA indices are ok. */
static int subset_leaf_by_Lindex(SEXP leaf, int dim0, SEXP Lindex,
		SEXP ans, CopyRVectorElt_FUNType fun, int *lookup_table)
{
	int n = LENGTH(Lindex);
	if (leaf == R_NilValue || n == 0)
		return 0;

	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);

	/* We use bsearch_idx0_to_k2() instead of the lookup table if 'Lindex'
	   is short. This avoids the overhead of building the lookup table and
	   so should be slightly more efficient.
	   TODO: Right now we use a cutoff value of 10 but this needs to be
	   refined. The overhead of building the lookup table depends on the
	   value of 'nzcount' so the cutoff value should be a function
	   of 'nzcount'. */
	int use_lookup_table = n > 10;
	if (use_lookup_table)
		build_lookup_table(lookup_table, INTEGER(nzoffs), nzcount);
	/* Walk on 'Lindex'. */
	for (int k1 = 0; k1 < n; k1++) {
		int idx0;
		int ret = extract_idx0(Lindex, k1, dim0, &idx0);
		if (ret < 0)
			return ret;
		if (ret == SUBSCRIPT_ELT_IS_NA) {
			/* 'Lindex[k1]' is NA or NaN. */
			set_Rvector_elt_to_NA(ans, (R_xlen_t) idx0);
			continue;
		}
		int k2;
		if (use_lookup_table) {
			k2 = lookup_table[idx0];
		} else {
			k2 = bsearch_idx0_to_k2(idx0, INTEGER(nzoffs), nzcount);
		}
		if (k2 >= 0)
			copy_nzval_elt(nzvals, k2, ans, (R_xlen_t) k1, fun);
	}
	if (use_lookup_table)
		reset_lookup_table(lookup_table, INTEGER(nzoffs), nzcount);
	return 0;
}

static void subset_leaf_by_OPBuf(SEXP leaf, const OPBuf *opbuf, SEXP ans,
		CopyRVectorElt_FUNType fun, int *lookup_table)
{
	if (leaf == R_NilValue || opbuf->nelt == 0)
		return;

	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
	const int *from_offs = opbuf->soffs;
	const R_xlen_t *to_offs = opbuf->loffs;

	/* We use bsearch_idx0_to_k2() instead of the lookup table if
	   the number of (from_off,to_off) pairs is small. This avoids the
	   overhead of building the lookup table and so should be slightly
	   more efficient.
	   TODO: Right now we use a cutoff value of 10 but this needs to be
	   refined. The overhead of building the lookup table depends on the
	   value of 'nzcount' so the cutoff value should be a function
	   of 'nzcount'. */
	int use_lookup_table = opbuf->nelt > 10;
	if (use_lookup_table)
		build_lookup_table(lookup_table, INTEGER(nzoffs), nzcount);
	/* Walk on the (from_off,to_off) pairs. */
	for (int k1 = 0; k1 < opbuf->nelt; k1++) {
		int from_off = from_offs[k1];
		int k2;
		if (use_lookup_table) {
			k2 = lookup_table[from_off];
		} else {
			k2 = bsearch_idx0_to_k2(from_off,
						INTEGER(nzoffs), nzcount);
		}
		if (k2 >= 0)
			copy_nzval_elt(nzvals, k2,
				       ans, to_offs[k1], fun);
	}
	if (use_lookup_table)
		reset_lookup_table(lookup_table, INTEGER(nzoffs), nzcount);
	return;
}

/* 'subscript' must be a numeric vector. It cannot be a long one! It is
   expected to contain 1-based indices that are >= 1 and <= 'sv->len'.
   NA indices will trigger an error.
   'sv_selection' and 'out_nzoffs' must be arrays that are long enough to
   hold at least 'LENGTH(subscript)' ints. */
static int subset_SV(const SparseVec *sv, SEXP subscript,
		     int *sv_selection, int *out_nzoffs,
		     int *lookup_table)
{
	int out_nzcount = 0;
	int n = LENGTH(subscript);
	if (n == 0)
		return out_nzcount;
	int sv_nzcount = get_SV_nzcount(sv);
	build_lookup_table(lookup_table, sv->nzoffs, sv_nzcount);
	/* Walk on 'subscript'. */
	for (int i1 = 0; i1 < n; i1++) {
		int idx0;
		int ret = extract_idx0(subscript, i1, sv->len, &idx0);
		if (ret != 0)
			bad_Nindex_error(ret, 1);
		//k2 = bsearch_idx0_to_k2(idx0, sv->nzoffs, sv_nzcount);
		int k2 = lookup_table[idx0];
		if (k2 >= 0) {
			sv_selection[out_nzcount] = k2;
			out_nzoffs[out_nzcount] = i1;
			out_nzcount++;
		}
	}
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
	SparseVec sv = leaf2SV(leaf, TYPEOF(leaf_nzvals), dim0);
	int ans_nzcount = subset_SV(&sv, subscript,
				    selection_buf, nzoffs_buf, lookup_table);
	if (ans_nzcount == 0)
		return R_NilValue;

	SEXP ans_nzoffs = PROTECT(NEW_INTEGER(ans_nzcount));
	memcpy(INTEGER(ans_nzoffs), nzoffs_buf, sizeof(int) * ans_nzcount);
	if (leaf_nzvals == R_NilValue) {
		/* Leaf to subset is lacunar --> subsetting preserves that. */
		SEXP ans = _make_lacunar_leaf(ans_nzoffs);
		UNPROTECT(1);
		return ans;
	}
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
	SEXP ans = zip_leaf(ans_nzvals, ans_nzoffs);
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * C_subset_SVT_by_[M|L]index()
 */

static int check_Mindex_dim(SEXP Mindex, int ndim,
		const char *what1, const char *what2)
{
	SEXP Mindex_dim = GET_DIM(Mindex);
	if (Mindex_dim == R_NilValue || LENGTH(Mindex_dim) != 2)
		error("'%s' must be a matrix", what1);
	if (!IS_INTEGER(Mindex))
		error("'%s' must be an integer matrix", what1);
	if (INTEGER(Mindex_dim)[1] != ndim)
		error("ncol(%s) != %s", what1, what2);
	return INTEGER(Mindex_dim)[0];
}

/* Three possible outcomes:
   (1) We land on a leaf node (OPBuf), in which case we append the loff/soff
       pair to the leaf node and return its new length, that is, the new nb
       of (loff,soff) pairs contained in it. This will always be > 0.
   (2) We didn't land anywhere, in which case we return 0.
   (3) We encountered an error, in which case we return a negative value. */
static int add_offset_pair_to_OPBufTree(R_xlen_t loff, R_xlen_t idx0,
		SEXP SVT, const int *dim, const R_xlen_t *dimcumprod, int ndim,
		OPBufTree *opbuf_tree)
{
	for (int along = ndim - 1; along >= 1; along--) {
		R_xlen_t p = dimcumprod[along - 1];
		int i = idx0 / p;  /* always >= 0 and < 'dim[along]' */
		SVT = VECTOR_ELT(SVT, i);
		if (SVT == R_NilValue)
			return 0;
		idx0 %= p;
		if (opbuf_tree->node_type == NULL_NODE)
			_alloc_OPBufTree_children(opbuf_tree, dim[along]);
		opbuf_tree = get_OPBufTree_child(opbuf_tree, i);
	}
	/* 'idx0' is guaranteed to be < dimcumprod[0] = dim[0] <= INT_MAX. */
	if (opbuf_tree->node_type == NULL_NODE)
		_alloc_OPBufTree_leaf(opbuf_tree);
	OPBuf *opbuf = get_OPBufTree_leaf(opbuf_tree);
	/* If 'opbuf->nelt' is INT_MAX then _append_to_OPBuf() will return
	   MAX_OPBUF_LEN_REACHED (negative value, see OPBufTree.h). Note that
	   this will only happen if 'Lindex' is a long vector and more than
	   INT_MAX in it hit the same leaf in 'SVT'. A rather crazy and
	   unlikely situation! */
	return _append_to_OPBuf(opbuf, loff, (int) idx0);
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
	int max_outleaf_len = 0;
	R_xlen_t out_len = XLENGTH(Lindex);
	/* Walk along 'Lindex'. */
	int ret = 0;
	for (R_xlen_t loff = 0; loff < out_len; loff++) {
		R_xlen_t idx0;
		ret = extract_long_idx0(Lindex, loff,
					dimcumprod[x_ndim - 1], &idx0);
		if (ret < 0)
			return ret;
		if (ret == SUBSCRIPT_ELT_IS_NA) {
			/* 'Lindex[loff]' is NA or NaN. */
			set_Rvector_elt_to_NA(ans, loff);
			continue;
		}
		ret = add_offset_pair_to_OPBufTree(loff, idx0,
					x_SVT, x_dim, dimcumprod, x_ndim,
					opbuf_tree);
		if (ret < 0)
			return ret;
		if (ret > max_outleaf_len)
			max_outleaf_len = ret;
	}
	return max_outleaf_len;
}

/* Recursive tree traversal must be guided by 'opbuf_tree', not by 'SVT'. See
   long comment above why. */
static void REC_subset_SVT_by_OPBufTree(OPBufTree *opbuf_tree,
		SEXP SVT, const int *dim, int ndim, SEXP ans,
		CopyRVectorElt_FUNType fun, int *lookup_table, int pardim)
{
	if (opbuf_tree->node_type == NULL_NODE)
		return;

	if (ndim == 1) {
		/* 'opbuf_tree' and 'SVT' are leaves. */
		OPBuf *opbuf = get_OPBufTree_leaf(opbuf_tree);
		subset_leaf_by_OPBuf(SVT, opbuf, ans, fun, lookup_table);
		_free_OPBufTree(opbuf_tree);
		return;
	}

	/* 'opbuf_tree' and 'SVT' are inner nodes. */
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

/* --- .Call ENTRY POINT --- */
SEXP C_subset_SVT_by_Mindex(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP Mindex)
{
	SEXPTYPE Rtype = _get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("SparseArray internal error in "
		      "C_subset_SVT_by_Mindex():\n"
		      "    SVT_SparseArray object has invalid type");
	int x_ndim = LENGTH(x_dim);
	int ans_len = check_Mindex_dim(Mindex, x_ndim,
				       "Mindex", "length(dim(x))");
	SEXP ans = PROTECT(_new_Rvector0(Rtype, (R_xlen_t) ans_len));

	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_subset_SVT_by_Lindex(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP Lindex)
{
	SEXPTYPE Rtype = _get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("SparseArray internal error in "
		      "C_subset_SVT_by_Lindex():\n"
		      "    SVT_SparseArray object has invalid type");

	int x_ndim = LENGTH(x_dim);
	if (!IS_INTEGER(Lindex) && !IS_NUMERIC(Lindex))
		error("'Lindex' must be an integer or numeric vector");

	R_xlen_t ans_len = XLENGTH(Lindex);
	SEXP ans = PROTECT(_new_Rvector0(Rtype, ans_len));
	if (x_SVT == R_NilValue) {
		UNPROTECT(1);
		return ans;
	}

	if (x_ndim == 1)
		error("x_ndim == 1 not ready yet");

	/* 1st pass: Build OPBufTree. */
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
		if (max_outleaf_len == BAD_SUBSCRIPT_TYPE)
			error("'Lindex' must be a numeric vector");
		if (max_outleaf_len == SUBSCRIPT_ELT_IS_LESS_THAN_ONE ||
		    max_outleaf_len == SUBSCRIPT_ELT_IS_BEYOND_MAX)
			error("'Lindex' contains out-of-bound indices");
		if (max_outleaf_len == MAX_OPBUF_LEN_REACHED)
			error("too many indices in 'Lindex' hit the same "
			      "leaf in the Sparse Vector Tree representation");
		error("SparseArray internal error in "
		      "C_subset_SVT_by_Lindex():\n"
		      "    unexpected error code %d", max_outleaf_len);
	}
	//double dt = (1.0 * clock() - t0) * 1000.0 / CLOCKS_PER_SEC;
	//printf("1st pass: %2.3f ms\n", dt);

	//printf("max_outleaf_len = %d\n", max_outleaf_len);

	/* 2nd pass: Subset SVT by OPBufTree. */
	if (max_outleaf_len > 0) {
		//clock_t t0 = clock();
		CopyRVectorElt_FUNType fun =
			_select_copy_Rvector_elt_FUN(Rtype);
		int x_dim0 = INTEGER(x_dim)[0];
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
			bad_Nindex_error(BAD_SUBSCRIPT_TYPE, along + 1);
		}
		R_xlen_t d = XLENGTH(subscript);
		if (d > INT_MAX) {
			UNPROTECT(1);
			bad_Nindex_error(SUBSCRIPT_IS_TOO_LONG, along + 1);
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
		int idx0;
		if (subscript == R_NilValue) {
			idx0 = i;
		} else {
			int ret = extract_idx0(subscript, i, SVT_len, &idx0);
			if (ret != 0)
				bad_Nindex_error(ret, ndim);
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
   'Nindex' must be an N-index, that is, a list of integer vectors (or NULLs),
   one along each dimension in the array. */
SEXP C_subset_SVT_by_Nindex(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP Nindex)
{
	SEXPTYPE Rtype = _get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("SparseArray internal error in "
		      "C_subset_SVT_by_Nindex():\n"
		      "    SVT_SparseArray object has invalid type");

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

