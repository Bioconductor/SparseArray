/****************************************************************************
 ****************************************************************************
 **									   **
 **        Single-bracket subsetting (`[`) of a SparseArray object         **
 **									   **
 ****************************************************************************
 ****************************************************************************/
#include "SparseArray_subsetting.h"

#include "OPBufTree.h"
#include "Rvector_utils.h"
#include "leaf_utils.h"

#include <limits.h>  /* for INT_MAX */
#include <string.h>  /* for memcpy() */


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
 * subset_leaf1()
 * subset_leaf2()
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

static void subset_leaf1(SEXP leaf,
		const int *from_offs, const R_xlen_t *to_offs, int noffs,
		SEXP ans,
		int *lookup_table, CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	if (leaf == R_NilValue || noffs == 0)
		return;
	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
	/* TODO: Try to use map_i2_to_k2_with_bsearch() instead of the lookup
	   table if the vector of from_off/to_off pairs is "short" (e.g.
	   noffs < 5). This would avoid the overhead of building the lookup
	   table and so might be slightly more efficient. */
	build_lookup_table(lookup_table, INTEGER(nzoffs), nzcount);
	/* Walk on the vector of from_off/to_off pairs. */
	for (int k1 = 0; k1 < noffs; k1++) {
		int from_off = from_offs[k1];
		//k2 = map_i2_to_k2_with_bsearch(from_off,
		//			INTEGER(nzoffs), nzcount);
		//k2 = map_i2_to_k2_with_lookup_table(from_off, lookup_table);
		int k2 = lookup_table[from_off];
		if (k2 >= 0) {
			if (nzvals == R_NilValue) {
				/* lacunar leaf */
				_set_Rsubvec_elts_to_one(ans,
					to_offs[k1], (R_xlen_t) 1);
			} else {
				/* standard leaf */
				copy_Rvector_elt_FUN(nzvals, (R_xlen_t) k2,
						     ans, to_offs[k1]);
			}
		}
	}
	reset_lookup_table(lookup_table, INTEGER(nzoffs), nzcount);
	return;
}

static inline int get_i2(const int *idx, int i1, int dim0)
{
	int i2;

	i2 = idx[i1];
	if (i2 == NA_INTEGER) {
		UNPROTECT(1);
		error("'Nindex' cannot contain NAs");
	}
	if (i2 < 1 || i2 > dim0) {
		UNPROTECT(1);
		error("'Nindex' contains out-of-bound "
		      "indices");
	}
	return --i2;
}

/* Takes a non-NULL leaf (standard or lacunar), and returns a leaf that can
   be NULL, standard, or lacunar. */
static SEXP subset_leaf2(SEXP leaf, int dim0, SEXP idx,
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
		int i2 = get_i2(INTEGER(idx), i1, dim0);
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

/* MAX_OPBUF_LEN_REACHED is set to -1 in OPBufTree.h so make sure to use
   a different negative value for INVALID_LINEAR_INDEX_VALUE. */
#define	INVALID_LINEAR_INDEX_VALUE -2
#define	LINEAR_INDEX_VALUE_IS_NA    1

static inline int get_Lindex_elt(SEXP Lindex, R_xlen_t i, R_xlen_t *out)
{
	if (IS_INTEGER(Lindex)) {
		int Lindex_elt = INTEGER(Lindex)[i];
		if (Lindex_elt == NA_INTEGER)
			return LINEAR_INDEX_VALUE_IS_NA;
		if (Lindex_elt < 1)
			return INVALID_LINEAR_INDEX_VALUE;
		*out = (R_xlen_t) Lindex_elt;
		return 0;
	}
	double Lindex_elt = REAL(Lindex)[i];
	/* ISNAN(): True for *both* NA and NaN. See <R_ext/Arith.h> */
	if (ISNAN(Lindex_elt))
		return LINEAR_INDEX_VALUE_IS_NA;
	if (Lindex_elt < 1 || Lindex_elt >= 1.00 + R_XLEN_T_MAX)
		return INVALID_LINEAR_INDEX_VALUE;
	*out = (R_xlen_t) Lindex_elt;
	return 0;
}

/* Three possible outcomes:
   (1) We land on a leaf node (OPBuf), in which case we append the loff/soff
       pair to the leaf node and return its new length, that is, the new nb
       of loff/soff pairs contained in it. This will always be > 0.
   (2) We didn't land anywhere, in which case we return 0.
   (3) We encountered an error, in which case we return a negative value. */
static int attach_Lindex_elt_to_OPBufTree(R_xlen_t Lindex_elt, R_xlen_t loff,
		SEXP SVT, const int *dim, const R_xlen_t *dimcumprod, int ndim,
		OPBufTree *opbuf_tree)
{
	R_xlen_t idx0 = Lindex_elt - 1;  /* 0-based */
	if (idx0 >= dimcumprod[ndim - 1])
		return INVALID_LINEAR_INDEX_VALUE;
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
	_free_OPBufTree(opbuf_tree);  /* reset 'opbuf_tree' to NULL_NODE */
	int max_outleaf_len = 0;
	R_xlen_t out_len = XLENGTH(Lindex);
	/* Walk along 'Lindex'. */
	int ret = 0;
	for (R_xlen_t loff = 0; loff < out_len; loff++) {
		R_xlen_t Lindex_elt;
		ret = get_Lindex_elt(Lindex, loff, &Lindex_elt);
		if (ret < 0)
			return ret;
		if (ret == LINEAR_INDEX_VALUE_IS_NA) {
			/* Lindex[loff] is NA or NaN. */
			set_Rvector_elt_to_NA(ans, loff);
			continue;
		}
		ret = attach_Lindex_elt_to_OPBufTree(Lindex_elt, loff,
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
		int *lookup_table, CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	if (opbuf_tree->node_type == NULL_NODE)
		return;

	if (ndim == 1) {
		/* 'opbuf_tree' and 'SVT' are leaves. */
		OPBuf *opbuf = get_OPBufTree_leaf(opbuf_tree);
		subset_leaf1(SVT,
			     opbuf->soffs, opbuf->loffs, opbuf->nelt, ans,
			     lookup_table, copy_Rvector_elt_FUN);
		_free_OPBufTree(opbuf_tree);
		return;
	}

	/* 'opbuf_tree' and 'SVT' are inner nodes. */
	int n = get_OPBufTree_nchildren(opbuf_tree);  /* same as dim[ndim - 1]
							 or LENGTH(SVT) */
	for (int i = 0; i < n; i++) {
		OPBufTree *child = get_OPBufTree_child(opbuf_tree, i);
		SEXP subSVT = VECTOR_ELT(SVT, i);
		REC_subset_SVT_by_OPBufTree(child,
				subSVT, dim, ndim - 1, ans,
				lookup_table, copy_Rvector_elt_FUN);
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
	R_xlen_t *dimcumprod = (R_xlen_t *) R_alloc(x_ndim, sizeof(R_xlen_t));
	R_xlen_t p = 1;
	for (int along = 0; along < x_ndim; along++) {
		p *= INTEGER(x_dim)[along];
		dimcumprod[along] = p;
	}
	OPBufTree *opbuf_tree = _get_global_opbuf_tree();
	int max_outleaf_len =
		build_OPBufTree_from_Lindex(opbuf_tree, Lindex,
				x_SVT, INTEGER(x_dim), x_ndim, ans,
				dimcumprod);
	if (max_outleaf_len < 0) {
		UNPROTECT(1);
		if (max_outleaf_len == INVALID_LINEAR_INDEX_VALUE)
			error("'Lindex' contains invalid linear indices");
		if (max_outleaf_len == MAX_OPBUF_LEN_REACHED)
			error("too many indices in 'Lindex' hit the same "
			      "leaf in the Sparse Vector Tree representation");
		error("SparseArray internal error in "
		      "C_subset_SVT_by_Lindex():\n"
		      "    unexpected error code %d", max_outleaf_len);
	}

	/* 2nd pass: Subset SVT by OPBufTree. */
	if (max_outleaf_len > 0) {
		int x_dim0 = INTEGER(x_dim)[0];
		int *lookup_table = (int *) R_alloc(x_dim0, sizeof(int));
		for (int i = 0; i < x_dim0; i++)
			lookup_table[i] = -1;
		CopyRVectorElt_FUNType fun =
			_select_copy_Rvector_elt_FUN(Rtype);
		REC_subset_SVT_by_OPBufTree(opbuf_tree,
				x_SVT, INTEGER(x_dim), x_ndim, ans,
				lookup_table, fun);
	}

	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_subset_SVT_by_Nindex()
 */

static SEXP compute_subset_dim(SEXP Nindex, SEXP x_dim)
{
	int ndim, along;
	SEXP ans_dim, Nindex_elt;
	R_xlen_t d;

	ndim = LENGTH(x_dim);
	if (!isVectorList(Nindex) || LENGTH(Nindex) != ndim)
		error("'Nindex' must be a list with one list "
		      "element along each dimension in 'x'");

	ans_dim = PROTECT(duplicate(x_dim));

	for (along = 0; along < ndim; along++) {
		Nindex_elt = VECTOR_ELT(Nindex, along);
		if (Nindex_elt == R_NilValue)
			continue;
		if (!IS_INTEGER(Nindex_elt)) {
			UNPROTECT(1);
			error("each list element in 'Nindex' must "
			      "be either NULL or an integer vector");
		}
		d = XLENGTH(Nindex_elt);
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

/* Recursive tree traversal.
   Returns R_NilValue or a list of length 'ans_dim[ndim - 1]'. */
static SEXP REC_subset_SVT_by_Nindex(SEXP SVT, SEXP Nindex,
		const int *x_dim, const int *ans_dim, int ndim,
		int *i1_buf, int *k2_buf, int *lookup_table)
{
	if (SVT == R_NilValue)
		return R_NilValue;

	SEXP idx = VECTOR_ELT(Nindex, ndim - 1);

	if (ndim == 1) {
		/* 'SVT' is a leaf. */
		return subset_leaf2(SVT, x_dim[0], idx,
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
		SEXP subSVT = VECTOR_ELT(SVT, i2);
		SEXP ans_elt = REC_subset_SVT_by_Nindex(subSVT, Nindex,
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
	int x_dim0 = INTEGER(x_dim)[0];
	int *lookup_table = (int *) R_alloc(x_dim0, sizeof(int));
	int ans_dim0 = INTEGER(ans_dim)[0];
	int *i1_buf = (int *) R_alloc(ans_dim0, sizeof(int));
	int *k2_buf = (int *) R_alloc(ans_dim0, sizeof(int));
	for (int i = 0; i < x_dim0; i++)
		lookup_table[i] = -1;
	SEXP ans_SVT = REC_subset_SVT_by_Nindex(x_SVT, Nindex,
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

