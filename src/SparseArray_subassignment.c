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


/* Copied from S4Arrays/src/array_selection.h */
#define INVALID_COORD(coord, maxcoord) \
	((coord) == NA_INTEGER || (coord) < 1 || (coord) > (maxcoord))

/* Maybe move this to src/OPBufTree.c */
static OPBuf R_alloc_OPBuf(int buflen)
{
	OPBuf opbuf;
	opbuf.buflen = buflen;
	opbuf.idx0s = (int *) R_alloc(buflen, sizeof(int));
	opbuf.Loffs = (int *) R_alloc(buflen, sizeof(int));
	opbuf.xLoffs = NULL;
	return opbuf;
}

static R_xlen_t *alloc_and_compute_cumprod(const int *x, int x_len)
{
	R_xlen_t *cumprod = (R_xlen_t *) R_alloc(x_len, sizeof(R_xlen_t));
	R_xlen_t prod = 1;
	for (int i = 0; i < x_len; i++) {
		prod *= x[i];
		cumprod[i] = prod;
	}
	return cumprod;
}


/****************************************************************************
 * subassign_leaf_by_OPBuf()
 */

static SEXP subassign_leaf_by_OPBuf(
		SEXP leaf, const OPBuf *opbuf, SEXP Rvector,
		OPBuf *sorted_opbuf,
		int *order_buf, unsigned short int *rxbuf1, int *rxbuf2,
		SparseVec *buf_sv)
{
	if (opbuf->xLoffs != NULL && sorted_opbuf->xLoffs == NULL)
		sorted_opbuf->xLoffs = (R_xlen_t *)
			R_alloc(sorted_opbuf->buflen, sizeof(R_xlen_t));
	_sort_and_remove_dups_OPBuf(opbuf, sorted_opbuf,
				    order_buf, rxbuf1, rxbuf2);
	if (opbuf->Loffs != NULL)
		return _subassign_leaf_with_Rvector_subset(leaf,
				sorted_opbuf->idx0s, sorted_opbuf->nelt,
				Rvector, sorted_opbuf->Loffs, buf_sv);
	if (opbuf->xLoffs != NULL)
		return _subassign_leaf_with_Rvector_xsubset(leaf,
				sorted_opbuf->idx0s, sorted_opbuf->nelt,
				Rvector, sorted_opbuf->xLoffs, buf_sv);
	error("SparseArray internal error in "
	      "subassign_leaf_by_OPBuf()\n"
	      "    'sorted_opbuf->Loffs' and 'sorted_opbuf->xLoffs' are NULL");
	return R_NilValue;  /* will never reach this */
}


/****************************************************************************
 * build_OPBufTree_from_Lindex()
 * build_OPBufTree_from_Mindex()
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

static OPBufTree *find_host_node_for_Mindex_row(OPBufTree *opbuf_tree,
		SEXP Mindex, int Mnrow, R_xlen_t Moff,
		const int *dim, int ndim,
		int *idx0, int *ret_code)
{
	for (int along = ndim - 1; along >= 1; along--, Moff -= Mnrow) {
		int i, d = dim[along];
		*ret_code = extract_idx0(Mindex, Moff, d, &i);
		if (*ret_code < 0)
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
   Returns a negative value in case of error. */
static int build_OPBufTree_from_Lindex1(OPBufTree *opbuf_tree, SEXP Lindex,
		const int *x_dim, int x_ndim,
		const R_xlen_t *dimcumprod)
{
	int max_opbuf_nelt = 0;
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
		if (ret > max_opbuf_nelt)
			max_opbuf_nelt = ret;
	}
	return max_opbuf_nelt;
}

static int build_OPBufTree_from_Lindex2(OPBufTree *opbuf_tree, SEXP Lindex,
		const int *x_dim, int x_ndim,
		const R_xlen_t *dimcumprod)
{
	int max_opbuf_nelt = 0;
	R_xlen_t in_len = XLENGTH(Lindex);
	R_xlen_t x_len = dimcumprod[x_ndim - 1];
	/* Walk along 'Lindex'. */
	for (R_xlen_t Loff = 0; Loff < in_len; Loff++) {
		R_xlen_t Lidx0;
		int ret = extract_long_idx0(Lindex, Loff, x_len, &Lidx0);
		if (ret < 0)
			return ret;
		int idx0;
		OPBufTree *host_node = find_host_node_for_Lidx0(
						opbuf_tree, Lidx0,
						x_dim, x_ndim,
						dimcumprod, &idx0);
		ret = _append_idx0xLoff_to_host_node(host_node, idx0, Loff);
		if (ret < 0)
			return ret;
		if (ret > max_opbuf_nelt)
			max_opbuf_nelt = ret;
	}
	return max_opbuf_nelt;
}

static int build_OPBufTree_from_Lindex(OPBufTree *opbuf_tree,
		SEXP Lindex, const int *x_dim, int x_ndim)
{
	/* _free_OPBufTree(opbuf_tree) resets 'opbuf_tree->node_type'
	   to NULL_NODE. */
	_free_OPBufTree(opbuf_tree);
	R_xlen_t *dimcumprod = alloc_and_compute_cumprod(x_dim, x_ndim);
	return XLENGTH(Lindex) <= (R_xlen_t) INT_MAX ?
		build_OPBufTree_from_Lindex1(opbuf_tree, Lindex,
				x_dim, x_ndim, dimcumprod) :
		build_OPBufTree_from_Lindex2(opbuf_tree, Lindex,
				x_dim, x_ndim, dimcumprod);
}

static int build_OPBufTree_from_Mindex(OPBufTree *opbuf_tree,
		SEXP Mindex, int Mnrow, const int *x_dim, int x_ndim)
{
	/* _free_OPBufTree(opbuf_tree) resets 'opbuf_tree->node_type'
	   to NULL_NODE. */
	_free_OPBufTree(opbuf_tree);
	int max_opbuf_nelt = 0;
	R_xlen_t Moff = (R_xlen_t) Mnrow * (x_ndim - 1);
	/* Walk along the rows of 'Mindex'. */
	for (int Loff = 0; Loff < Mnrow; Loff++, Moff++) {
		int idx0, ret;
		OPBufTree *host_node = find_host_node_for_Mindex_row(
						opbuf_tree,
						Mindex, Mnrow, Moff,
						x_dim, x_ndim,
						&idx0, &ret);
		if (ret < 0)
			return ret;
		ret = _append_idx0Loff_to_host_node(host_node, idx0, Loff);
		if (ret < 0)
			return ret;
		if (ret > max_opbuf_nelt)
			max_opbuf_nelt = ret;
	}
	return max_opbuf_nelt;
}


/****************************************************************************
 * subassign_SVT_by_OPBufTree()
 */

/* Recursive tree traversal of 'opbuf_tree'. */
static SEXP REC_subassign_SVT_by_OPBufTree(OPBufTree *opbuf_tree,
		SEXP SVT, int ndim, SEXP vals,
		OPBuf *sorted_opbuf,
		int *order_buf, unsigned short int *rxbuf1, int *rxbuf2,
		SparseVec *buf_sv, int pardim)
{
	if (opbuf_tree->node_type == NULL_NODE)
		return SVT;

	if (ndim == 1) {
		/* Both 'opbuf_tree' and 'SVT' are leaves. */
		OPBuf *opbuf = get_OPBufTree_leaf(opbuf_tree);
		SEXP ans = subassign_leaf_by_OPBuf(SVT, opbuf, vals,
					sorted_opbuf,
					order_buf, rxbuf1, rxbuf2, buf_sv);
		/* PROTECT not really necessary since _free_OPBufTree()
		   won't trigger R's garbage collector but this could change
		   someday so we'd better not take any risk. */
		PROTECT(ans);
		_free_OPBufTree(opbuf_tree);
		UNPROTECT(1);
		return ans;
	}

	/* Both 'opbuf_tree' and 'SVT' are inner nodes.
	   'n' is their outermost dimension. */
	int n = get_OPBufTree_nchildren(opbuf_tree);
	SEXP ans = PROTECT(NEW_LIST(n));
	int is_empty = 1;
	for (int i = 0; i < n; i++) {
		OPBufTree *child = get_OPBufTree_child(opbuf_tree, i);
		SEXP subSVT = SVT == R_NilValue ? R_NilValue
						: VECTOR_ELT(SVT, i);
		SEXP ans_elt = REC_subassign_SVT_by_OPBufTree(child,
					subSVT, ndim - 1, vals,
					sorted_opbuf,
					order_buf, rxbuf1, rxbuf2,
					buf_sv, pardim);
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

static SEXP subassign_SVT_by_OPBufTree(SEXP SVT, const int *dim, int ndim,
		SEXPTYPE Rtype, int na_background,
		OPBufTree *opbuf_tree, int max_opbuf_nelt, SEXP vals)
{
	int dim0 = dim[0];

	SparseVec buf_sv = _alloc_buf_SparseVec(Rtype, dim0, na_background, 0);
	if (IS_STRSXP_OR_VECSXP(buf_sv.Rtype))
		PROTECT(buf_sv.nzvals);

	int buflen = max_opbuf_nelt < dim0 ? max_opbuf_nelt : dim0;
	OPBuf sorted_opbuf = R_alloc_OPBuf(buflen);

	int *order_buf = (int *) R_alloc(max_opbuf_nelt, sizeof(int));

	unsigned short int *rxbuf1 = (unsigned short int *)
			R_alloc(max_opbuf_nelt, sizeof(unsigned short int));

	int *rxbuf2 = (int *) R_alloc(max_opbuf_nelt, sizeof(int));

	/* Get 1-based rank of biggest dimension (ignoring the 1st dim).
	   Parallel execution will be along that dimension. */
	int pardim = which_max(dim + 1, ndim - 1) + 2;

	SEXP ans = REC_subassign_SVT_by_OPBufTree(opbuf_tree,
				 SVT, ndim, vals,
				 &sorted_opbuf,
				 order_buf, rxbuf1, rxbuf2,
				 &buf_sv, pardim);

	if (IS_STRSXP_OR_VECSXP(buf_sv.Rtype))
		UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_subassign_SVT_by_Lindex()
 */

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

	int ndim = LENGTH(x_dim);
	R_xlen_t nvals = XLENGTH(vals);
	if (XLENGTH(Lindex) != nvals)
		error("length(Lindex) != length(vals)");
	if (nvals == 0)
		return x_SVT;  /* no-op */

	/* --- STEP 1: Build the OPBufTree --- */

	OPBufTree *opbuf_tree = _get_global_opbuf_tree();
	int max_opbuf_nelt = build_OPBufTree_from_Lindex(opbuf_tree, Lindex,
							 INTEGER(x_dim), ndim);
	if (max_opbuf_nelt < 0)
		_bad_Lindex_error(max_opbuf_nelt);

	//printf("max_opbuf_nelt = %d\n", max_opbuf_nelt);
	//_print_OPBufTree(opbuf_tree, 1);

	/* --- STEP 2: Subassign SVT by OPBufTree --- */

	return subassign_SVT_by_OPBufTree(x_SVT, INTEGER(x_dim), ndim,
			Rtype, x_has_NAbg,
			opbuf_tree, max_opbuf_nelt, vals);
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
	if (!(IS_INTEGER(Mindex) || IS_NUMERIC(Mindex)))
		error("'%s' must be a numeric matrix", what1);
	if (INTEGER(Mindex_dim)[0] != nvals)
		error("nrow(%s) != %s", what1, what2);
	if (INTEGER(Mindex_dim)[1] != ndim)
		error("ncol(%s) != %s", what1, what3);
	return;
}

/* --- .Call ENTRY POINT --- */
SEXP C_subassign_SVT_by_Mindex(
		SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP x_na_background,
		SEXP Mindex, SEXP vals)
{
	SEXPTYPE Rtype = _get_and_check_Rtype_from_Rstring(x_type,
				"C_subassign_SVT_by_Mindex", "x_type");
	if (TYPEOF(vals) != Rtype)
		error("SparseArray internal error in "
		      "C_subassign_SVT_by_Mindex():\n"
		      "    SVT_SparseArray object and 'vals' "
		      "must have the same type");

	int x_has_NAbg = _get_and_check_na_background(x_na_background,
				"C_subassign_SVT_by_Mindex", "x_na_background");

	int ndim = LENGTH(x_dim);
	R_xlen_t nvals = XLENGTH(vals);
	check_Mindex_dim(Mindex, nvals, ndim,
			 "Mindex", "length(vals)", "length(dim(x))");
	if (nvals == 0)
		return x_SVT;  /* no-op */

	/* --- STEP 1: Build the OPBufTree --- */

	OPBufTree *opbuf_tree = _get_global_opbuf_tree();
	int max_opbuf_nelt = build_OPBufTree_from_Mindex(opbuf_tree,
					     Mindex, (int) nvals,
					     INTEGER(x_dim), ndim);
	if (max_opbuf_nelt < 0)
		_bad_Mindex_error(max_opbuf_nelt);

	/* --- STEP 2: Subassign SVT by OPBufTree --- */

	return subassign_SVT_by_OPBufTree(x_SVT, INTEGER(x_dim), ndim,
			Rtype, x_has_NAbg,
			opbuf_tree, max_opbuf_nelt, vals);
}


/****************************************************************************
 * Some helpers shared between C_subassign_SVT_with_short_Rvector(),
 * C_subassign_SVT_with_Rarray(), and C_subassign_SVT_with_SVT()
 */

static int check_offs(SEXP offs, int d)
{
	int n = LENGTH(offs);
	const int *offs_p = INTEGER(offs);
	int prev_off = -1;
	for (int i = 0; i < n; i++) {
		int off = offs_p[i];
		if (off == NA_INTEGER)
			error("subscripts contain NAs");
		if (off < 0 || off >= d)
			error("subscripts contain out-of-bound indices");
		if (off <= prev_off)
			error("SparseArray internal error in check_offs():\n"
			      "    subscripts are not strictly sorted");
		prev_off = off;
	}
	return n;
}

static int check_Noffs(SEXP Noffs, const int *dim, int ndim, const int *arr_dim)
{
	if (LENGTH(Noffs) != ndim)
		error("SparseArray internal error in check_Noffs():\n"
		      "    'Noffs' must have one list element per "
		      "dimension in the array to subassign");
	for (int along = 0; along < ndim; along++) {
		if (dim[along] == 0)
			return 1;  /* subassignment is a no-op */
		SEXP offs = VECTOR_ELT(Noffs, along);
		int doas;  /* dim of array selection */
		if (offs == R_NilValue) {
			doas = dim[along];
		} else if (IS_INTEGER(offs)) {
			doas = check_offs(offs, dim[along]);
		} else {
			error("subscripts must be integer vectors");
		}
		if (doas == 0)
			return 1;  /* subassignment is a no-op */
		if (arr_dim != NULL && arr_dim[along] != doas)
			error("SparseArray internal error in check_Noffs():\n"
			      "    dimensions of right array don't "
			      "match dimensions of array selection");
	}
	return 0;
}

/* Offsets were already checked upfront by check_Noffs() above. */
static int get_off(SEXP offs, int i)
{
	return offs == R_NilValue ? i : INTEGER(offs)[i];
}

static SEXP new_SVT(int d, SEXP SVT0)
{
	if (d == 0)
		error("SparseArray internal error in new_SVT():\n"
		      "    d == 0");
	SEXP SVT = PROTECT(NEW_LIST(d));
	if (SVT0 != R_NilValue) {
		if (!isVectorList(SVT0))  // IS_LIST() is broken
			error("SparseArray internal error in new_SVT():\n"
			      "    'SVT0' is not a list");
		if (LENGTH(SVT0) != d)
			error("SparseArray internal error in new_SVT():\n"
			      "    'LENGTH(SVT0) != d'");
		/* Shallow copy. */
		for (int i = 0; i < d; i++)
			SET_VECTOR_ELT(SVT, i, VECTOR_ELT(SVT0, i));
	}
	UNPROTECT(1);
	return SVT;
}

static SEXP post_process_SVT(SEXP SVT, SEXP SVT0)
{
	int SVT_len = LENGTH(SVT), is_empty = 1;
	for (int i = 0; i < SVT_len; i++) {
		if (VECTOR_ELT(SVT, i) != R_NilValue) {
			is_empty = 0;
			break;
		}
	}
	if (is_empty)
		return R_NilValue;
	if (SVT0 == R_NilValue)
		return SVT;
	int is_noop = 1;
	for (int i = 0; i < SVT_len; i++) {
		if (VECTOR_ELT(SVT, i) != VECTOR_ELT(SVT0, i)) {
			is_noop = 0;
			break;
		}
	}
	return is_noop ? SVT0 : SVT;
}


/****************************************************************************
 * C_subassign_SVT_with_short_Rvector()
 */

static int *is_full_replacement(SEXP Noffs, const int *dim, int ndim)
{
	int *fully = (int *) R_alloc(ndim, sizeof(int));
	int ok = 1;
	for (int along = 0; along < ndim; along++) {
		SEXP offs = VECTOR_ELT(Noffs, along);
		ok = ok && (offs == R_NilValue || LENGTH(offs) == dim[along]);
		fully[along] = ok;
	}
	return fully;
}

/* Recursive. */
static SEXP REC_subassign_SVT_with_short_Rvector(SEXP SVT,
		const int *dim, int ndim, SEXP Noffs, SEXP short_Rvector,
		const int *fully, SparseVec *buf_sv)
{
	SEXP offs = VECTOR_ELT(Noffs, ndim - 1);
	if (ndim == 1)
		return _subassign_leaf_with_Rvector_block(SVT, offs,
				offs == R_NilValue ? dim[0] : LENGTH(offs),
				short_Rvector, 0, buf_sv);
	SEXP ans_shared_elt = R_NilValue;
	int use_shared = fully[ndim - 2];
	if (use_shared) {
		SEXP subSVT = SVT == R_NilValue ? R_NilValue :
					VECTOR_ELT(SVT, get_off(offs, 0));
		ans_shared_elt =
			REC_subassign_SVT_with_short_Rvector(subSVT,
					dim, ndim - 1, Noffs, short_Rvector,
					fully, buf_sv);
		if (ans_shared_elt != R_NilValue) {
			PROTECT(ans_shared_elt);
		} else if (fully[ndim - 1]) {
			return R_NilValue;
		}
	}
	int d1 = dim[ndim - 1];
	int d2 = offs == R_NilValue ? d1 : LENGTH(offs);
	SEXP ans = PROTECT(new_SVT(d1, SVT));
	for (int i2 = 0; i2 < d2; i2++) {
		int i1 = get_off(offs, i2);
		if (use_shared) {
			SET_VECTOR_ELT(ans, i1, ans_shared_elt);
			continue;
		}
		SEXP subSVT = VECTOR_ELT(ans, i1);
		SEXP ans_elt = PROTECT(
			REC_subassign_SVT_with_short_Rvector(subSVT,
					dim, ndim - 1, Noffs, short_Rvector,
					fully, buf_sv)
		);
		SET_VECTOR_ELT(ans, i1, ans_elt);
		UNPROTECT(1);
	}
	if (use_shared && ans_shared_elt != R_NilValue)
		UNPROTECT(1);
	ans = post_process_SVT(ans, SVT);
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT ---
   'Noffs' must be a list of integer vectors (or NULLs), one along each
   dimension in the arrays. Each non-NULL list element must contain valid
   offsets (i.e. zero-based indices) along the corresponding dimension in 'x'.
   The offsets must be sorted in **strictly** ascending order. */
SEXP C_subassign_SVT_with_short_Rvector(
		SEXP x_dim, SEXP x_type, SEXP x_SVT,
		SEXP Noffs, SEXP Rvector)
{
	SEXPTYPE x_Rtype = _get_and_check_Rtype_from_Rstring(x_type,
			     "C_subassign_SVT_with_short_Rvector", "x_type");
	if (TYPEOF(Rvector) != x_Rtype)
		error("SparseArray internal error in "
		      "C_subassign_SVT_with_short_Rvector():\n"
		      "    SVT_SparseArray object and 'Rvector' "
		      "must have the same type");

	const int *dim = INTEGER(x_dim);
	int ndim = LENGTH(x_dim);
	if (check_Noffs(Noffs, dim, ndim, NULL))
		return x_SVT;  /* no-op */

	const int *fully = is_full_replacement(Noffs, dim, ndim);

	/* _subassign_leaf_with_Rvector_block(), the workhorse behind
	   REC_subassign_SVT_with_short_Rvector(), does not always need
	   'buf_sv.nzvals' and is able to allocate it the first time it
	   needs it (if it ever needs it). So we set the 'nzoffs_only'
	   argument to 1 in our _alloc_buf_SparseVec() call below. */
	SparseVec buf_sv = _alloc_buf_SparseVec(x_Rtype, dim[0], 0, 1);
	if (IS_STRSXP_OR_VECSXP(buf_sv.Rtype))
		buf_sv.nzvals = PROTECT(
			allocVector(buf_sv.Rtype, (R_xlen_t) buf_sv.len)
		);
	SEXP ans = REC_subassign_SVT_with_short_Rvector(x_SVT, dim, ndim,
						Noffs, Rvector,
						fully, &buf_sv);
	if (IS_STRSXP_OR_VECSXP(buf_sv.Rtype))
		UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_subassign_SVT_with_Rarray()
 */

static SEXP check_Rarray(SEXP Rarray, int ndim, SEXPTYPE x_Rtype)
{
	SEXP Rarray_dim = GET_DIM(Rarray);
	if (Rarray_dim == R_NilValue)
		error("SparseArray internal error in check_Rarray():\n"
		      "    'Rarray' must be an array");
	if (LENGTH(Rarray_dim) != ndim)
		error("SparseArray internal error in check_Rarray():\n"
		      "    SVT_SparseArray object and 'Rarray' "
		      "must have the same number of dimensions");
	if (TYPEOF(Rarray) != x_Rtype)
		error("SparseArray internal error in check_Rarray():\n"
		      "    SVT_SparseArray object and 'Rarray' "
		      "must have the same type");
	return Rarray_dim;
}

/* Recursive. */
static SEXP REC_subassign_SVT_with_Rsubarr(SEXP SVT,
		const int *dim, int ndim, SEXP Noffs,
		SEXP Rarray, R_xlen_t arr_offset, const R_xlen_t *subarr_lens,
		SparseVec *buf_sv)
{
	SEXP offs = VECTOR_ELT(Noffs, ndim - 1);
	if (ndim == 1)
		return _subassign_leaf_with_Rvector_block(SVT,
					offs, subarr_lens[0],
					Rarray, arr_offset, buf_sv);
	int d1 = dim[ndim - 1];
	int d2 = offs == R_NilValue ? d1 : LENGTH(offs);
	SEXP ans = PROTECT(new_SVT(d1, SVT));
	R_xlen_t offset_inc = subarr_lens[ndim - 2];
	for (int i2 = 0; i2 < d2; i2++, arr_offset += offset_inc) {
		int i1 = get_off(offs, i2);
		SEXP subSVT = VECTOR_ELT(ans, i1);
		SEXP ans_elt = PROTECT(
			REC_subassign_SVT_with_Rsubarr(subSVT,
					dim, ndim - 1, Noffs,
					Rarray, arr_offset, subarr_lens,
					buf_sv)
		);
		SET_VECTOR_ELT(ans, i1, ans_elt);
		UNPROTECT(1);
	}
	ans = post_process_SVT(ans, SVT);
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT ---
   The left and right arrays ('x' and 'Rarray') must have the same number
   of dimensions.
   'Noffs' must be a list of integer vectors (or NULLs), one along each
   dimension in the arrays. Each non-NULL list element must contain valid
   offsets (i.e. zero-based indices) along the corresponding dimension in 'x'.
   The offsets must be sorted in **strictly** ascending order. */
SEXP C_subassign_SVT_with_Rarray(
		SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP x_na_background,
		SEXP Noffs, SEXP Rarray)
{
	SEXPTYPE x_Rtype = _get_and_check_Rtype_from_Rstring(x_type,
			     "C_subassign_SVT_with_Rarray", "x_type");
	int x_has_NAbg = _get_and_check_na_background(x_na_background,
			     "C_subassign_SVT_with_Rarray", "x_na_background");

	int ndim = LENGTH(x_dim);
	SEXP Rarray_dim = check_Rarray(Rarray, ndim, x_Rtype);

	const int *dim = INTEGER(x_dim);
	const int *arr_dim = INTEGER(Rarray_dim);
	if (check_Noffs(Noffs, dim, ndim, arr_dim))
		return x_SVT;  /* no-op */

	/* _subassign_leaf_with_Rvector_block(), the workhorse behind
	   REC_subassign_SVT_with_Rsubarr(), does not always need
	   'buf_sv.nzvals' and is able to allocate it the first time it
	   needs it (if it ever needs it). So we set the 'nzoffs_only'
	   argument to 1 in our _alloc_buf_SparseVec() call below. */
	SparseVec buf_sv = _alloc_buf_SparseVec(x_Rtype, dim[0], x_has_NAbg, 1);
	if (IS_STRSXP_OR_VECSXP(buf_sv.Rtype))
		buf_sv.nzvals = PROTECT(
			allocVector(buf_sv.Rtype, (R_xlen_t) buf_sv.len)
		);
	R_xlen_t *subarr_lens = alloc_and_compute_cumprod(arr_dim, ndim);
	SEXP ans = REC_subassign_SVT_with_Rsubarr(x_SVT, dim, ndim, Noffs,
						  Rarray, 0, subarr_lens,
						  &buf_sv);
	if (IS_STRSXP_OR_VECSXP(buf_sv.Rtype))
		UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_subassign_SVT_with_SVT()
 */

/* Recursive. */
static SEXP REC_subassign_SVT1_with_SVT2(
		SEXP SVT1, const int *dim1, int ndim, SEXP Noffs,
		SEXP SVT2, const int *dim2, SparseVec *buf_sv)
{
	if (SVT1 == R_NilValue && SVT2 == R_NilValue)
		return R_NilValue;
	SEXP offs = VECTOR_ELT(Noffs, ndim - 1);
	if (ndim == 1)
		return _subassign_leaf_with_leaf(SVT1, offs, dim2[0],
						 SVT2, buf_sv);
	int d1 = dim1[ndim - 1];
	int d2 = dim2[ndim - 1];
	SEXP ans = PROTECT(new_SVT(d1, SVT1));
	for (int i2 = 0; i2 < d2; i2++) {
		int i1 = get_off(offs, i2);
		SEXP subSVT1 = VECTOR_ELT(ans, i1);
		SEXP subSVT2 = SVT2 == R_NilValue ? R_NilValue :
						    VECTOR_ELT(SVT2, i2);
		SEXP ans_elt = PROTECT(
			REC_subassign_SVT1_with_SVT2(subSVT1,
					dim1, ndim - 1, Noffs,
					subSVT2, dim2, buf_sv)
		);
		SET_VECTOR_ELT(ans, i1, ans_elt);
		UNPROTECT(1);
	}
	ans = post_process_SVT(ans, SVT1);
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT ---
   The left and right arrays ('x' and 'y') must have the same number
   of dimensions.
   'Noffs' must be a list of integer vectors (or NULLs), one along each
   dimension in the arrays. Each non-NULL list element must contain valid
   offsets (i.e. zero-based indices) along the corresponding dimension in 'x'.
   The offsets must be sorted in **strictly** ascending order. */
SEXP C_subassign_SVT_with_SVT(
		SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP x_na_background,
		SEXP Noffs,
		SEXP y_dim, SEXP y_type, SEXP y_SVT, SEXP y_na_background)
{
	SEXPTYPE x_Rtype = _get_and_check_Rtype_from_Rstring(x_type,
			     "C_subassign_SVT_with_SVT", "x_type");
	int x_has_NAbg = _get_and_check_na_background(x_na_background,
			     "C_subassign_SVT_with_SVT", "x_na_background");
	SEXPTYPE y_Rtype = _get_and_check_Rtype_from_Rstring(y_type,
			     "C_subassign_SVT_with_SVT", "y_type");
	int y_has_NAbg = _get_and_check_na_background(y_na_background,
			     "C_subassign_SVT_with_SVT", "y_na_background");
	if (x_Rtype != y_Rtype)
		error("SparseArray internal error in "
		      "C_subassign_SVT_with_SVT():\n"
		      "    x_Rtype != y_Rtype");
	if (x_has_NAbg != y_has_NAbg)
		error("SparseArray internal error in "
		      "C_subassign_SVT_with_SVT():\n"
		      "    x_has_NAbg != y_has_NAbg");

	int ndim = LENGTH(x_dim);
	if (LENGTH(y_dim) != ndim)
		error("SparseArray internal error in "
		      "C_subassign_SVT_with_SVT():\n"
		      "    LENGTH(x_dim) != LENGTH(y_dim)");

	if (check_Noffs(Noffs, INTEGER(x_dim), ndim, INTEGER(y_dim)))
		return x_SVT;  /* no-op */

	SparseVec buf_sv = _alloc_buf_SparseVec(x_Rtype, INTEGER(x_dim)[0],
						x_has_NAbg, 0);
	if (IS_STRSXP_OR_VECSXP(buf_sv.Rtype))
		PROTECT(buf_sv.nzvals);
	SEXP ans = REC_subassign_SVT1_with_SVT2(x_SVT,
						INTEGER(x_dim), ndim, Noffs,
						y_SVT, INTEGER(y_dim),
						&buf_sv);
	if (IS_STRSXP_OR_VECSXP(buf_sv.Rtype))
		UNPROTECT(1);
	return ans;
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
	int along, doas;
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
			doas = dim[along];
		} else if (IS_INTEGER(Nindex_elt)) {
			doas = LENGTH(Nindex_elt);
		} else {
			error("subscripts must be integer vectors");
		}
		selection_len *= doas;
		if (along < margin)
			continue;
		Nindex_iter->selection_dim[along - margin] = doas;
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

