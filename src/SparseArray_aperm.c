/****************************************************************************
 *                   Permutation of a SparseArray object                    *
 ****************************************************************************/
#include "SparseArray_aperm.h"

#include "Rvector_utils.h"
#include "leaf_utils.h"

#include <string.h>  /* for memset() */


/****************************************************************************
 * C_transpose_2D_SVT()
 *
 * TODO: Do we still need this? C_aperm_SVT() below reimplements this for
 * the multidimensional case and with no significant/measurable differences
 * in terms of performance for the 2D case compared to C_transpose_2D_SVT().
 */

static void collect_stats_on_rows(SEXP SVT, int nrow, int ncol,
		int *nzcount_buf, int *onecount_buf)
{
	memset(nzcount_buf, 0, sizeof(int) * nrow);
	if (onecount_buf != NULL)
		memset(onecount_buf, 0, sizeof(int) * nrow);
	for (int j = 0; j < ncol; j++) {
		SEXP leaf = VECTOR_ELT(SVT, j);
		if (leaf == R_NilValue)
			continue;
		SEXP nzvals, nzoffs;
		int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
		const int *nzoffs_p = INTEGER(nzoffs);
		for (int k = 0; k < nzcount; k++, nzoffs_p++) {
			nzcount_buf[*nzoffs_p]++;
			if (onecount_buf == NULL)
				continue;
			if (nzvals == R_NilValue ||  /* lacunar leaf */
			    _all_Rsubvec_elts_equal_one(nzvals, 0, 1))
				onecount_buf[*nzoffs_p]++;
		}
	}
	return;
}

static void **set_quick_out_nzvals_p(SEXP out_SVT, SEXPTYPE Rtype)
{
	int out_SVT_len = LENGTH(out_SVT);
	switch (Rtype) {
	    case LGLSXP: case INTSXP: {
		int **vals_p = (int **) R_alloc(out_SVT_len, sizeof(int *));
		int **p = vals_p;
		for (int i = 0; i < out_SVT_len; i++, p++) {
			SEXP leaf = VECTOR_ELT(out_SVT, i);
			if (leaf == R_NilValue)
				continue;
			SEXP nzvals = get_leaf_nzvals(leaf);
			*p = nzvals == R_NilValue ? NULL : INTEGER(nzvals);
		}
		return (void **) vals_p;
	    }
	    case REALSXP: {
		double **vals_p = (double **)
			R_alloc(out_SVT_len, sizeof(double *));
		double **p = vals_p;
		for (int i = 0; i < out_SVT_len; i++, p++) {
			SEXP leaf = VECTOR_ELT(out_SVT, i);
			if (leaf == R_NilValue)
				continue;
			SEXP nzvals = get_leaf_nzvals(leaf);
			*p = nzvals == R_NilValue ? NULL : REAL(nzvals);
		}
		return (void **) vals_p;
	    }
	    case CPLXSXP: {
		Rcomplex **vals_p = (Rcomplex **)
			R_alloc(out_SVT_len, sizeof(Rcomplex *));
		Rcomplex **p = vals_p;
		for (int i = 0; i < out_SVT_len; i++, p++) {
			SEXP leaf = VECTOR_ELT(out_SVT, i);
			if (leaf == R_NilValue)
				continue;
			SEXP nzvals = get_leaf_nzvals(leaf);
			*p = nzvals == R_NilValue ? NULL : COMPLEX(nzvals);
		}
		return (void **) vals_p;
	    }
	    case RAWSXP: {
		Rbyte **vals_p = (Rbyte **)
			R_alloc(out_SVT_len, sizeof(Rbyte *));
		Rbyte **p = vals_p;
		for (int i = 0; i < out_SVT_len; i++, p++) {
			SEXP leaf = VECTOR_ELT(out_SVT, i);
			if (leaf == R_NilValue)
				continue;
			SEXP nzvals = get_leaf_nzvals(leaf);
			*p = nzvals == R_NilValue ? NULL : RAW(nzvals);
		}
		return (void **) vals_p;
	    }
	}
	/* STRSXP and VECSXP cases. */
	return NULL;
}

typedef void (*TransposeCol_FUNType)(int col_idx, SEXP col,
		void **quick_out_nzvals_p, int **quick_out_nzoffs_p,
		SEXP out_SVT, int *nzcount_buf);

/* Ignores 'out_SVT' and 'nzcount_buf'. */
static void transpose_INTEGER_col(int col_idx, SEXP col,
		void **quick_out_nzvals_p, int **quick_out_nzoffs_p,
		SEXP out_SVT, int *nzcount_buf)
{
	SEXP col_nzvals, col_nzoffs;
	int col_nzcount = unzip_leaf(col, &col_nzvals, &col_nzoffs);
	const int *col_nzvals_p = NULL;  /* -Wmaybe-uninitialized */
	int v;
	if (col_nzvals != R_NilValue) {
		col_nzvals_p = INTEGER(col_nzvals);
	} else {
		v = int1;
	}
	int **out_nzvals_p = (int **) quick_out_nzvals_p;
	const int *col_nzoffs_p = INTEGER(col_nzoffs);
	for (int k = 0; k < col_nzcount; k++) {
		int row_idx = *col_nzoffs_p;
		int **p = out_nzvals_p + row_idx;
		if (*p != NULL) {
			if (col_nzvals != R_NilValue)
				v = col_nzvals_p[k];
			*((*p)++) = v;
		}
		*(quick_out_nzoffs_p[row_idx]++) = col_idx;
		col_nzoffs_p++;
	}
	return;
}

/* Ignores 'out_SVT' and 'nzcount_buf'. */
static void transpose_NUMERIC_col(int col_idx, SEXP col,
		void **quick_out_nzvals_p, int **quick_out_nzoffs_p,
		SEXP out_SVT, int *nzcount_buf)
{
	SEXP col_nzvals, col_nzoffs;
	int col_nzcount = unzip_leaf(col, &col_nzvals, &col_nzoffs);
	const double *col_nzvals_p = NULL;  /* -Wmaybe-uninitialized */
	double v;
	if (col_nzvals != R_NilValue) {
		col_nzvals_p = REAL(col_nzvals);
	} else {
		v = double1;
	}
	double **out_nzvals_p = (double **) quick_out_nzvals_p;
	const int *col_nzoffs_p = INTEGER(col_nzoffs);
	for (int k = 0; k < col_nzcount; k++) {
		int row_idx = *col_nzoffs_p;
		double **p = out_nzvals_p + row_idx;
		if (*p != NULL) {
			if (col_nzvals != R_NilValue)
				v = col_nzvals_p[k];
			*((*p)++) = v;
		}
		*(quick_out_nzoffs_p[row_idx]++) = col_idx;
		col_nzoffs_p++;
	}
	return;
}

/* Ignores 'out_SVT' and 'nzcount_buf'. */
static void transpose_COMPLEX_col(int col_idx, SEXP col,
		void **quick_out_nzvals_p, int **quick_out_nzoffs_p,
		SEXP out_SVT, int *nzcount_buf)
{
	SEXP col_nzvals, col_nzoffs;
	int col_nzcount = unzip_leaf(col, &col_nzvals, &col_nzoffs);
	const Rcomplex *col_nzvals_p = NULL;  /* -Wmaybe-uninitialized */
	Rcomplex v;
	if (col_nzvals != R_NilValue) {
		col_nzvals_p = COMPLEX(col_nzvals);
	} else {
		v = Rcomplex1;
	}
	Rcomplex **out_nzvals_p = (Rcomplex **) quick_out_nzvals_p;
	const int *col_nzoffs_p = INTEGER(col_nzoffs);
	for (int k = 0; k < col_nzcount; k++) {
		int row_idx = *col_nzoffs_p;
		Rcomplex **p = out_nzvals_p + row_idx;
		if (*p != NULL) {
			if (col_nzvals != R_NilValue)
				v = col_nzvals_p[k];
			*((*p)++) = v;
		}
		*(quick_out_nzoffs_p[row_idx]++) = col_idx;
		col_nzoffs_p++;
	}
	return;
}

/* Ignores 'out_SVT' and 'nzcount_buf'. */
static void transpose_RAW_col(int col_idx, SEXP col,
		void **quick_out_nzvals_p, int **quick_out_nzoffs_p,
		SEXP out_SVT, int *nzcount_buf)
{
	SEXP col_nzvals, col_nzoffs;
	int col_nzcount = unzip_leaf(col, &col_nzvals, &col_nzoffs);
	const Rbyte *col_nzvals_p = NULL;  /* -Wmaybe-uninitialized */
	Rbyte v;
	if (col_nzvals != R_NilValue) {
		col_nzvals_p = RAW(col_nzvals);
	} else {
		v = Rbyte1;
	}
	Rbyte **out_nzvals_p = (Rbyte **) quick_out_nzvals_p;
	const int *col_nzoffs_p = INTEGER(col_nzoffs);
	for (int k = 0; k < col_nzcount; k++) {
		int row_idx = *col_nzoffs_p;
		Rbyte **p = out_nzvals_p + row_idx;
		if (*p != NULL) {
			if (col_nzvals != R_NilValue)
				v = col_nzvals_p[k];
			*((*p)++) = v;
		}
		*(quick_out_nzoffs_p[row_idx]++) = col_idx;
		col_nzoffs_p++;
	}
	return;
}

/* Ignores 'quick_out_nzvals_p'. */
static void transpose_CHARACTER_col(int col_idx, SEXP col,
		void **quick_out_nzvals_p, int **quick_out_nzoffs_p,
		SEXP out_SVT, int *nzcount_buf)
{
	SEXP col_nzvals, col_nzoffs;
	int col_nzcount = unzip_leaf(col, &col_nzvals, &col_nzoffs);
	const int *col_nzoffs_p = INTEGER(col_nzoffs);
	for (int k = 0; k < col_nzcount; k++) {
		int row_idx = *col_nzoffs_p;
		SEXP out_leaf = VECTOR_ELT(out_SVT, row_idx);
		_copy_CHARACTER_elt(col_nzvals, (R_xlen_t) k,
			get_leaf_nzvals(out_leaf),
			(R_xlen_t) nzcount_buf[row_idx]++);
		*(quick_out_nzoffs_p[row_idx]++) = col_idx;
		col_nzoffs_p++;
	}
	return;
}

/* Ignores 'quick_out_nzvals_p'. */
static void transpose_LIST_col(int col_idx, SEXP col,
		void **quick_out_nzvals_p, int **quick_out_nzoffs_p,
		SEXP out_SVT, int *nzcount_buf)
{
	SEXP col_nzvals, col_nzoffs;
	int col_nzcount = unzip_leaf(col, &col_nzvals, &col_nzoffs);
	const int *col_nzoffs_p = INTEGER(col_nzoffs);
	for (int k = 0; k < col_nzcount; k++) {
		int row_idx = *col_nzoffs_p;
		SEXP out_leaf = VECTOR_ELT(out_SVT, row_idx);
		_copy_LIST_elt(col_nzvals, (R_xlen_t) k,
			get_leaf_nzvals(out_leaf),
			(R_xlen_t) nzcount_buf[row_idx]++);
		*(quick_out_nzoffs_p[row_idx]++) = col_idx;
		col_nzoffs_p++;
	}
	return;
}

static TransposeCol_FUNType select_transpose_col_FUN(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return transpose_INTEGER_col;
	    case REALSXP:             return transpose_NUMERIC_col;
	    case CPLXSXP:             return transpose_COMPLEX_col;
	    case RAWSXP:              return transpose_RAW_col;
	    case STRSXP:              return transpose_CHARACTER_col;
	    case VECSXP:              return transpose_LIST_col;
	}
	return NULL;
}

static SEXP transpose_2D_SVT(SEXP SVT, int nrow, int ncol, SEXPTYPE Rtype,
		int *nzcount_buf, int *onecount_buf)
{
	TransposeCol_FUNType transpose_col_FUN =
			select_transpose_col_FUN(Rtype);
	if (transpose_col_FUN == NULL)
		error("SparseArray internal error in "
		      "transpose_2D_SVT():\n"
		      "    SVT_SparseMatrix object has invalid type");

	/* 1st pass: Collect number of nonzero elements (and ones) for each
	   row in the input matrix. */
	collect_stats_on_rows(SVT, nrow, ncol, nzcount_buf, onecount_buf);

	/* 2nd pass: Build the transposed SVT. */
	SEXP ans = PROTECT(NEW_LIST(nrow));
	int **quick_out_nzoffs_p = (int **) R_alloc(nrow, sizeof(int *));
	for (int i = 0; i < nrow; i++) {
		int nzcount = nzcount_buf[i];
		if (nzcount == 0)
			continue;
		SEXP ans_nzoffs = PROTECT(NEW_INTEGER(nzcount));
		quick_out_nzoffs_p[i] = INTEGER(ans_nzoffs);
		SEXP ans_nzvals;
		if (onecount_buf != NULL && onecount_buf[i] == nzcount) {
			/* Create a lacunar leaf. Will turn it into a
			   standard leaf later if LACUNAR_MODE_IS_ON
			   is set to 0. */
			ans_nzvals = R_NilValue;
		} else {
			ans_nzvals = PROTECT(allocVector(Rtype, nzcount));
		}
		SEXP ans_elt = PROTECT(zip_leaf(ans_nzvals, ans_nzoffs));
		SET_VECTOR_ELT(ans, i, ans_elt);
		UNPROTECT(ans_nzvals == R_NilValue ? 2 : 3);
	}
	void **quick_out_nzvals_p = set_quick_out_nzvals_p(ans, Rtype);

	memset(nzcount_buf, 0, sizeof(int) * nrow);
	for (int j = 0; j < ncol; j++) {
		SEXP leaf = VECTOR_ELT(SVT, j);
		if (leaf == R_NilValue)
			continue;
		transpose_col_FUN(j, leaf,
			quick_out_nzvals_p, quick_out_nzoffs_p,
			ans, nzcount_buf);
	}

	if (onecount_buf != NULL && !LACUNAR_MODE_IS_ON) {
		/* Turn lacunar leaves into standard leaves. */
		for (int i = 0; i < nrow; i++) {
			SEXP ans_leaf = VECTOR_ELT(ans, i);
			if (ans_leaf == R_NilValue)
				continue;
			SEXP nzvals, nzoffs;
			int nzcount = unzip_leaf(ans_leaf, &nzvals, &nzoffs);
			if (nzvals != R_NilValue)
				continue;
			nzvals = PROTECT(_new_Rvector1(Rtype, nzcount));
			replace_leaf_nzvals(ans_leaf, nzvals);
			UNPROTECT(1);
		}
	}

	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_transpose_2D_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT)
{
	SEXPTYPE x_Rtype = _get_Rtype_from_Rstring(x_type);
	if (x_Rtype == 0)
		error("SparseArray internal error in "
		      "C_transpose_2D_SVT():\n"
		      "    SVT_SparseMatrix object has invalid type");

	if (LENGTH(x_dim) != 2)
		error("object to transpose must have exactly 2 dimensions");

	if (x_SVT == R_NilValue)
		return x_SVT;

	int x_nrow = INTEGER(x_dim)[0];
	int x_ncol = INTEGER(x_dim)[1];
	int *nzcount_buf = (int *) R_alloc(x_nrow, sizeof(int));
	int *onecount_buf = NULL;
	if (x_Rtype != STRSXP && x_Rtype != VECSXP)
		onecount_buf = (int *) R_alloc(x_nrow, sizeof(int));
	return transpose_2D_SVT(x_SVT, x_nrow, x_ncol, x_Rtype,
				nzcount_buf, onecount_buf);
}


/****************************************************************************
 * check_perm() and compute_ans_dim_and_perm_margins()
 *
 * The permutation of dimensions is described by an integer vector 'perm' of
 * length 'ndim' where 'ndim' is the number of dimensions of the array to
 * permute. The inner and outer margins are the number of inner and outer
 * dimensions that are not touched by the permutation. For example:
 *
 *     perm             inner_margin  outer_margin
 *     -------------------------------------------
 *     c(2,1)                      0             0
 *     c(2,1,3)                    0             1
 *     c(3,1,2)                    0             0
 *     c(1,2,3,4)                  4             4
 *     c(1,5,3,4,2,6,7)            1             2
 *
 * Note that normally 'inner_margin' + 'outer_margin' <= 'ndim' - 2, except
 * when the permutation is the identity in which case both 'inner_margin'
 * and 'outer_margin' are equal to 'ndim'.
 */

static void check_perm(SEXP perm, int ndim)
{
	if (!IS_INTEGER(perm))
		error("'perm' must be an integer vector");
	if (LENGTH(perm) != ndim)
		error("'length(perm)' not equal to number "
		      "of dimensions of array to permute");
	int *p_is_taken = (int *) R_alloc(ndim, sizeof(int));
	memset(p_is_taken, 0, sizeof(int) * ndim);
	for (int along = 0; along < ndim; along++) {
		int p = INTEGER(perm)[along];
		if (p == NA_INTEGER || p < 1 || p > ndim)
			error("invalid 'perm' argument");
		p--;
		if (p_is_taken[p])
			error("'perm' cannot contain duplicates");
		p_is_taken[p] = 1;
	}
	return;
}

static void compute_ans_dim_and_perm_margins(
		const int *dim, int ndim, const int *perm,
		int *ans_dim, int *inner_margin, int *outer_margin)
{
	int along;
	*inner_margin = ndim;
	for (along = 0; along < ndim; along++) {
		int p = perm[along] - 1;
		ans_dim[along] = dim[p];
		if (*inner_margin == ndim && p != along)
			*inner_margin = along;
	}
	/* Compute the outer margin. */
	for (along = ndim - 1; along >= 0; along--) {
		int p = perm[along] - 1;
		if (p != along)
			break;
	}
	*outer_margin = ndim - (along + 1);
	//printf("inner_margin = %d -- outer_margin = %d\n",
	//       *inner_margin, *outer_margin);
	return;
}


/****************************************************************************
 * Buffers used by aperm_with_new_SVT_leaves()
 */

static void *alloc_quick_out_nzvals_p(unsigned long long int n, SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return R_alloc(n, sizeof(int *));
	    case REALSXP:             return R_alloc(n, sizeof(double *));
	    case CPLXSXP:             return R_alloc(n, sizeof(Rcomplex *));
	    case RAWSXP:              return R_alloc(n, sizeof(Rbyte *));
	    case STRSXP: case VECSXP: return R_alloc(n, sizeof(SEXP));
	}
	error("SparseArray internal error in alloc_quick_out_nzvals_p():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return NULL;  /* will never reach this */
}

static inline void init_quick_out_nzvals_p(void *p, SEXPTYPE Rtype, SEXP vals)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: *((int      **) p) = DATAPTR(vals); break;
	    case REALSXP:             *((double   **) p) = DATAPTR(vals); break;
	    case CPLXSXP:             *((Rcomplex **) p) = DATAPTR(vals); break;
	    case RAWSXP:              *((Rbyte    **) p) = DATAPTR(vals); break;
	    case STRSXP: case VECSXP: *((SEXP      *) p) = vals; break;
	}
	return;
}

static inline void *inc_quick_out_nzvals_p(void *p, SEXPTYPE Rtype,
					 unsigned long long int inc)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return ((int      **) p) + inc;
	    case REALSXP:             return ((double   **) p) + inc;
	    case CPLXSXP:             return ((Rcomplex **) p) + inc;
	    case RAWSXP:              return ((Rbyte    **) p) + inc;
	    case STRSXP: case VECSXP: return ((SEXP      *) p) + inc;
	}
	return NULL;
}

typedef struct aperm0_bufs_t {
	int *nzcount_buf;  // an array of dimensions tail(dim(x)[perm], n=-1)
	unsigned long long int nzcount_buf_len;
	unsigned long long int *nzcount_buf_incs;
	unsigned long long int *outer_incs;
	void *quick_out_nzvals_p;
	int **quick_out_nzoffs_p;
} Aperm0Bufs;

/* Will always be called with 'ndim' >= 2. */
static void init_A0Bufs(Aperm0Bufs *A0Bufs, const int *dim, int ndim,
			const int *perm, SEXPTYPE Rtype)
{
	unsigned long long int *nzcount_buf_incs;
	unsigned long long int *outer_incs;
	unsigned long long int nzcount_buf_len;
	int along, p;

	nzcount_buf_incs = (unsigned long long int *)
		R_alloc(ndim - 1, sizeof(unsigned long long int));
	outer_incs = (unsigned long long int *)
		R_alloc(ndim, sizeof(unsigned long long int));
	outer_incs[perm[0] - 1] = 0;
	nzcount_buf_len = 1;
	for (along = 1; along < ndim; along++) {
		p = perm[along] - 1;
		nzcount_buf_incs[along - 1] = nzcount_buf_len;
		outer_incs[p] = nzcount_buf_len;
		nzcount_buf_len *= dim[p];
	}
/*
	printf("nzcount_buf_len = %llu\n", nzcount_buf_len);
	for (along = 0; along < ndim - 1; along++)
		printf("nzcount_buf_incs[%d] = %llu\n",
			along, nzcount_buf_incs[along]);
	for (along = 0; along < ndim; along++)
		printf("outer_incs[%d] = %llu\n",
			along, outer_incs[along]);
*/
	A0Bufs->nzcount_buf = (int *)
		R_alloc(nzcount_buf_len, sizeof(int));
	A0Bufs->nzcount_buf_len = nzcount_buf_len;
	A0Bufs->nzcount_buf_incs = nzcount_buf_incs;
	A0Bufs->outer_incs = outer_incs;
	A0Bufs->quick_out_nzvals_p =
		alloc_quick_out_nzvals_p(nzcount_buf_len, Rtype);
	A0Bufs->quick_out_nzoffs_p = (int **)
		R_alloc(nzcount_buf_len, sizeof(int *));
	return;
}


/****************************************************************************
 * alloc_aperm0_SVT_bufs()
 */

/* Depending on whether the 'perm' vector has an inner margin or not,
   we allocate the 'coords0_buf' buffer or the 'A0Bufs' buffers.
   The former is used by aperm_with_same_SVT_leaves() and the latter
   are used by aperm_with_new_SVT_leaves(). */
static int *alloc_aperm0_SVT_bufs(const int *dim, int ndim, SEXPTYPE Rtype,
		const int *perm, int inner_margin, Aperm0Bufs *A0Bufs)
{
	if (perm[0] == 1) {
		/* 'ndim - inner_margin' is guaranteed to be >= 2. */
		return (int *) R_alloc(ndim - inner_margin, sizeof(int));
	}
	/* 'ndim' is guaranteed to be >= 2. */
	init_A0Bufs(A0Bufs, dim, ndim, perm, Rtype);
	return NULL;
}


/****************************************************************************
 * EASY CASE: aperm_with_same_SVT_leaves()
 */

static void graft_subSVT_onto_ans(SEXP subSVT, const int *perm,
				  const int *ans_dim, int ans_ndim,
				  int inner_margin, const int *coords0_buf,
				  SEXP ans)
{
	for (int along = ans_ndim - 2; along >= inner_margin; along--) {
		int i = coords0_buf[perm[along + 1] - 1 - inner_margin];
		SEXP ans_elt = VECTOR_ELT(ans, i);
		if (ans_elt == R_NilValue) {
			ans_elt = PROTECT(NEW_LIST(ans_dim[along]));
			SET_VECTOR_ELT(ans, i, ans_elt);
			UNPROTECT(1);
		}
		ans = ans_elt;
	}
	int i = coords0_buf[perm[inner_margin] - 1 - inner_margin];
	/* Sanity check. */
	if (VECTOR_ELT(ans, i) != R_NilValue)
		error("SparseArray internal error in "
		      "graft_subSVT_onto_ans():\n"
		      "    graft spot is already taken");
	SET_VECTOR_ELT(ans, i, subSVT);
	return;
}

/* Recursive.
   'SVT' cannot be R_NilValue.
   Length of 'perm' array is 'ans_ndim'.
   Length of 'coords0_buf' array is 'ans_ndim' - 'inner_margin'. */
static void REC_aperm_with_same_SVT_leaves(
		SEXP SVT, const int *dim, int ndim, const int *perm,
		const int *ans_dim, int ans_ndim,
		int inner_margin, int *coords0_buf, SEXP ans)
{
	int SVT_len = LENGTH(SVT);
	for (int i = 0; i < SVT_len; i++) {
		SEXP subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue)
			continue;
		coords0_buf[ndim - inner_margin - 1] = i;
		if (ndim > inner_margin + 1) {
			REC_aperm_with_same_SVT_leaves(
					subSVT, dim, ndim - 1, perm,
					ans_dim, ans_ndim,
					inner_margin, coords0_buf, ans);
		} else {
			graft_subSVT_onto_ans(
					subSVT, perm,
					ans_dim, ans_ndim,
					inner_margin, coords0_buf, ans);
		}
	}
	return;
}

static SEXP aperm_with_same_SVT_leaves(
		SEXP SVT, const int *dim, int ndim,
		const int *perm, const int *ans_dim,
		int inner_margin, int *coords0_buf)
{
	int ans_len;
	SEXP ans;

	if (SVT == R_NilValue)
		return R_NilValue;

	ans_len = ans_dim[ndim - 1];
	ans = PROTECT(NEW_LIST(ans_len));
	REC_aperm_with_same_SVT_leaves(SVT, dim, ndim, perm,
				       ans_dim, ndim,
				       inner_margin, coords0_buf, ans);
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * HARD CASE: aperm_with_new_SVT_leaves()
 */

static inline void count_leaf_nzvals(SEXP leaf,
		unsigned long long int outer_inc0,
		unsigned long long int outer_offset0,
		int *nzcount_buf)
{
	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
	const int *nzoffs_p = INTEGER(nzoffs);
	for (int k = 0; k < nzcount; k++) {
		unsigned long long int outer_idx = outer_offset0 +
						   outer_inc0 * *nzoffs_p;
		nzcount_buf[outer_idx]++;
		nzoffs_p++;
	}
	return;
}

/* Recursive. */
static void REC_count_SVT_nzvals(SEXP SVT, int ndim,
		const unsigned long long int *outer_incs,
		unsigned long long int outer_offset0,
		int *nzcount_buf)
{
	unsigned long long int outer_inc;
	int SVT_len, i;
	SEXP subSVT;

	if (SVT == R_NilValue)
		return;
	outer_inc = outer_incs[ndim - 1];
	if (ndim == 1) {
		count_leaf_nzvals(SVT, outer_inc, outer_offset0, nzcount_buf);
		return;
	}
	SVT_len = LENGTH(SVT);
	for (i = 0; i < SVT_len; i++, outer_offset0 += outer_inc) {
		subSVT = VECTOR_ELT(SVT, i);
		REC_count_SVT_nzvals(subSVT, ndim - 1,
				     outer_incs, outer_offset0,
				     nzcount_buf);
	}
	return;
}

/* Recursive. */
static SEXP REC_grow_tree_and_alloc_leaves(const int *dim, int ndim,
		SEXPTYPE Rtype,
		const unsigned long long int *nzcount_buf_incs,
		const int *nzcount_buf,
		void *quick_out_nzvals_p, int **quick_out_nzoffs_p)
{
	if (ndim == 1) {
		if (*nzcount_buf == 0)
			return R_NilValue;
		SEXP nzvals, nzoffs;
		SEXP ans = PROTECT(
			_alloc_and_unzip_leaf(Rtype, *nzcount_buf,
					      &nzvals, &nzoffs)
		);
		*quick_out_nzoffs_p = INTEGER(nzoffs);
		init_quick_out_nzvals_p(quick_out_nzvals_p, Rtype, nzvals);
		UNPROTECT(1);
		return ans;
	}
	int ans_len = dim[ndim - 1];
	unsigned long long int buf_inc = nzcount_buf_incs[ndim - 2];
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	int is_empty = 1;
	for (int i = 0; i < ans_len; i++) {
		SEXP ans_elt = REC_grow_tree_and_alloc_leaves(dim, ndim - 1,
					Rtype, nzcount_buf_incs, nzcount_buf,
					quick_out_nzvals_p, quick_out_nzoffs_p);
		if (ans_elt != R_NilValue) {
			PROTECT(ans_elt);
			SET_VECTOR_ELT(ans, i, ans_elt);
			UNPROTECT(1);
			is_empty = 0;
		}
		nzcount_buf += buf_inc;
		quick_out_nzoffs_p += buf_inc;
		quick_out_nzvals_p = inc_quick_out_nzvals_p(quick_out_nzvals_p,
							    Rtype, buf_inc);
	}
	UNPROTECT(1);
	return is_empty ? R_NilValue : ans;
}

#define	FUNDEF_spray_ans_with_nzvals(type)				\
	(SEXP nzvals, const int *nzoffs,				\
	 unsigned long long int outer_inc0,				\
	 unsigned long long int outer_offset0,				\
	 int *nzcount_buf,						\
	 void *quick_out_nzvals_p, int **quick_out_nzoffs_p,		\
	 int inner_idx)							\
{									\
	int n = LENGTH(nzvals);						\
	const type *nzvals_p = DATAPTR(nzvals);				\
	type **out_nzvals_p = quick_out_nzvals_p;			\
	for (int k = 0; k < n; k++) {					\
		unsigned long long int outer_idx =			\
				outer_offset0 +	outer_inc0 * *nzoffs;	\
		int k2 = nzcount_buf[outer_idx]++;			\
		out_nzvals_p[outer_idx][k2] = *nzvals_p;		\
		nzvals_p++;						\
		quick_out_nzoffs_p[outer_idx][k2] = inner_idx;		\
		nzoffs++;						\
	}								\
	return;								\
}

static inline void spray_ans_with_ints
	FUNDEF_spray_ans_with_nzvals(int)
static inline void spray_ans_with_doubles
	FUNDEF_spray_ans_with_nzvals(double)
static inline void spray_ans_with_Rcomplexes
	FUNDEF_spray_ans_with_nzvals(Rcomplex)
static inline void spray_ans_with_Rbytes
	FUNDEF_spray_ans_with_nzvals(Rbyte)

static inline void spray_ans_with_CHARSXPs(
		SEXP nzvals, const int *nzoffs,
		unsigned long long int outer_inc0,
		unsigned long long int outer_offset0,
		int *nzcount_buf,
		void *quick_out_nzvals_p, int **quick_out_nzoffs_p,
		int inner_idx)
{
	int n = LENGTH(nzvals);
	SEXP *out_nzvals_p = quick_out_nzvals_p;
	for (int k = 0; k < n; k++) {
		unsigned long long int outer_idx =
				outer_offset0 + outer_inc0 * *nzoffs;
		int k2 = nzcount_buf[outer_idx]++;
		SET_STRING_ELT(out_nzvals_p[outer_idx], k2,
			       STRING_ELT(nzvals, k));
		quick_out_nzoffs_p[outer_idx][k2] = inner_idx;
		nzoffs++;
	}
	return;
}

static inline void spray_ans_with_list_elts(
		SEXP nzvals, const int *nzoffs,
		unsigned long long int outer_inc0,
		unsigned long long int outer_offset0,
		int *nzcount_buf,
		void *quick_out_nzvals_p, int **quick_out_nzoffs_p,
		int inner_idx)
{
	int n = LENGTH(nzvals);
	SEXP *out_nzvals_p = quick_out_nzvals_p;
	for (int k = 0; k < n; k++) {
		unsigned long long int outer_idx =
				outer_offset0 + outer_inc0 * *nzoffs;
		int k2 = nzcount_buf[outer_idx]++;
		SET_VECTOR_ELT(out_nzvals_p[outer_idx], k2,
			       VECTOR_ELT(nzvals, k));
		quick_out_nzoffs_p[outer_idx][k2] = inner_idx;
		nzoffs++;
	}
	return;
}

static void spray_ans_with_leaf(SEXP leaf, SEXPTYPE Rtype,
		unsigned long long int outer_inc0,
		unsigned long long int outer_offset0,
		int *nzcount_buf,
		void *quick_out_nzvals_p, int **quick_out_nzoffs_p,
		int inner_idx)
{
	void (*spray_FUN)(SEXP nzvals, const int *nzoffs,
			  unsigned long long int outer_inc0,
			  unsigned long long int outer_offset0,
			  int *nzcount_buf,
			  void *quick_out_nzvals_p, int **quick_out_nzoffs_p,
			  int inner_idx);
	SEXP nzvals, nzoffs;
	unzip_leaf(leaf, &nzvals, &nzoffs); /* ignore returned nzcount */
	switch (Rtype) {
	    case LGLSXP: case INTSXP: spray_FUN = spray_ans_with_ints; break;
	    case REALSXP: spray_FUN = spray_ans_with_doubles; break;
	    case CPLXSXP: spray_FUN = spray_ans_with_Rcomplexes; break;
	    case RAWSXP: spray_FUN = spray_ans_with_Rbytes; break;
	    case STRSXP: spray_FUN = spray_ans_with_CHARSXPs; break;
	    case VECSXP: spray_FUN = spray_ans_with_list_elts; break;
	    default:
		error("SparseArray internal error in spray_ans_with_leaf():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));
	}
	spray_FUN(nzvals, INTEGER(nzoffs),
		  outer_inc0, outer_offset0,
		  nzcount_buf, quick_out_nzvals_p, quick_out_nzoffs_p,
		  inner_idx);
	return;
}

/* Recursive. */
static void REC_fill_leaves(SEXP SVT, int ndim, SEXPTYPE Rtype,
		const unsigned long long int *outer_incs,
		unsigned long long int outer_offset0,
		int *nzcount_buf,
		void *quick_out_nzvals_p, int **quick_out_nzoffs_p)
{
	static int inner_idx;

	if (SVT == R_NilValue)
		return;
	unsigned long long int outer_inc = outer_incs[ndim - 1];
	if (ndim == 1) {
		spray_ans_with_leaf(SVT, Rtype,
			outer_inc, outer_offset0,
			nzcount_buf, quick_out_nzvals_p, quick_out_nzoffs_p,
			inner_idx);
		return;
	}
	int SVT_len = LENGTH(SVT);
	for (int i = 0; i < SVT_len; i++, outer_offset0 += outer_inc) {
		/* 'outer_inc == 0' means we're looping along the dimension
		   of SVT that becomes the **leftmost** dimension (a.k.a.
		   innermost dimension) after permutation, that is the leftmost
		   dimension in 'ans'. */
		if (outer_inc == 0)
			inner_idx = i;
		SEXP subSVT = VECTOR_ELT(SVT, i);
		REC_fill_leaves(subSVT, ndim - 1, Rtype,
				outer_incs, outer_offset0,
				nzcount_buf,
				quick_out_nzvals_p, quick_out_nzoffs_p);
	}
	return;
}

static SEXP aperm_with_new_SVT_leaves(
		SEXP SVT, const int *dim, int ndim, SEXPTYPE Rtype,
		const int *ans_dim, Aperm0Bufs *A0Bufs)
{
	SEXP ans;

	/* 1st pass */
	memset(A0Bufs->nzcount_buf, 0, sizeof(int) * A0Bufs->nzcount_buf_len);
	REC_count_SVT_nzvals(SVT, ndim, A0Bufs->outer_incs, 0,
			     A0Bufs->nzcount_buf);
	/* 2nd pass */
	ans = PROTECT(REC_grow_tree_and_alloc_leaves(ans_dim, ndim, Rtype,
						A0Bufs->nzcount_buf_incs,
						A0Bufs->nzcount_buf,
						A0Bufs->quick_out_nzvals_p,
						A0Bufs->quick_out_nzoffs_p));
	/* 3rd pass */
	memset(A0Bufs->nzcount_buf, 0, sizeof(int) * A0Bufs->nzcount_buf_len);
	REC_fill_leaves(SVT, ndim, Rtype,
			A0Bufs->outer_incs, 0,
			A0Bufs->nzcount_buf,
			A0Bufs->quick_out_nzvals_p,
			A0Bufs->quick_out_nzoffs_p);
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_aperm0_SVT()
 */

/* Can handle any array permutation. However, it is optimized for the "no
   outer margin" case, that is, when the permutation vector has no outer
   margin. When the permutation vector has a nonzero outer margin,
   aperm0_SVT() won't be as efficient as C_aperm_SVT().
   In other words, C_aperm_SVT() should always be used instead of
   aperm0_SVT() as it will always handle the general case optimally.
   C_aperm_SVT() is based on aperm0_SVT().
   Note that the bigger the outer margin, the more efficient C_aperm_SVT()
   will be compared to aperm0_SVT(). */
static SEXP aperm0_SVT(SEXP SVT, const int *dim, int ndim, SEXPTYPE Rtype,
		const int *perm, const int *ans_dim,
		int inner_margin, int *coords0_buf, Aperm0Bufs *A0Bufs)
{
	if (perm[0] == 1) {
		/* EASY CASE: Leaves will be preserved i.e. the permuted
		   array will reuse the leaves of the original array.
		   It is very memory efficient because the original leaves
		   can be reused as-is without even the need to copy them! */
		return aperm_with_same_SVT_leaves(SVT, dim, ndim,
						  perm, ans_dim,
						  inner_margin, coords0_buf);
	}

	/* HARD CASE: Leaves won't be preserved i.e. the permuted array will
	   use entirely new leaves. It is much less efficient than the EASY
	   CASE because memory needs to be allocated for the new leaves. */
	return aperm_with_new_SVT_leaves(SVT, dim, ndim, Rtype,
					 ans_dim, A0Bufs);
}

/* --- .Call ENTRY POINT ---
 * Just a wrapper for aperm0_SVT(). Only provided for convenience testing of
 * aperm0_SVT() from the R command line.
 */
SEXP C_aperm0_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP perm)
{
	SEXPTYPE Rtype = _get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("SparseArray internal error in C_aperm0_SVT():\n"
		      "    SVT_SparseArray object has invalid type");

	int ndim = LENGTH(x_dim);
	const int *dim = INTEGER(x_dim);
	check_perm(perm, ndim);
	int *ans_dim = (int *) R_alloc(ndim, sizeof(int));
	int inner_margin, outer_margin;
	compute_ans_dim_and_perm_margins(dim, ndim, INTEGER(perm),
					 ans_dim, &inner_margin, &outer_margin);
	if (outer_margin == ndim || x_SVT == R_NilValue)
		return x_SVT;

	/* Will allocate either the 'coords0_buf' buffer or the 'A0Bufs'
	   buffers. 'ndim' is guaranteed to be >= 2. */
	Aperm0Bufs A0Bufs;
	int *coords0_buf = alloc_aperm0_SVT_bufs(dim, ndim,
					    Rtype, INTEGER(perm),
					    inner_margin, &A0Bufs);
	return aperm0_SVT(x_SVT, dim, ndim, Rtype,
			  INTEGER(perm), ans_dim,
			  inner_margin, coords0_buf, &A0Bufs);
}


/****************************************************************************
 * C_aperm_SVT()
 */

/* Recursive. */
static SEXP REC_aperm_SVT(SEXP SVT, const int *dim, int ndim, SEXPTYPE Rtype,
		const int *perm, const int *ans_dim,
		int inner_margin, int *coords0_buf, Aperm0Bufs *A0Bufs)
{
	if (perm[ndim - 1] != ndim) {
		/* Outer margin is 0. */
		return aperm0_SVT(SVT, dim, ndim, Rtype,
				  perm, ans_dim,
				  inner_margin, coords0_buf, A0Bufs);
	}

	/* Outer margin is > 0. */
	int ans_len = LENGTH(SVT);
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	for (int i = 0; i < ans_len; i++) {
		SEXP subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue)
			continue;
		SEXP ans_elt = PROTECT(
			REC_aperm_SVT(subSVT, dim, ndim - 1, Rtype,
				      perm, ans_dim,
				      inner_margin, coords0_buf, A0Bufs)
		);
		SET_VECTOR_ELT(ans, i, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT ---
 * Handles any array permutation optimally, regardless of the outer margin
 * of the permutation vector.
 */
SEXP C_aperm_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP perm)
{
	SEXPTYPE Rtype = _get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("SparseArray internal error in C_aperm_SVT():\n"
		      "    SVT_SparseArray object has invalid type");

	int ndim = LENGTH(x_dim);
	const int *dim = INTEGER(x_dim);
	check_perm(perm, ndim);
	int *ans_dim = (int *) R_alloc(ndim, sizeof(int));
	int inner_margin, outer_margin;
	compute_ans_dim_and_perm_margins(dim, ndim, INTEGER(perm),
					 ans_dim, &inner_margin, &outer_margin);
	if (outer_margin == ndim || x_SVT == R_NilValue)
		return x_SVT;

	/* Will allocate either the 'coords0_buf' buffer or the 'A0Bufs'
	   buffers. 'ndim - outer_margin' is guaranteed to be >= 2. */
	Aperm0Bufs A0Bufs;
	int *coords0_buf = alloc_aperm0_SVT_bufs(dim, ndim - outer_margin,
					    Rtype, INTEGER(perm),
					    inner_margin, &A0Bufs);
	return REC_aperm_SVT(x_SVT, dim, ndim, Rtype,
			     INTEGER(perm), ans_dim,
			     inner_margin, coords0_buf, &A0Bufs);
}

