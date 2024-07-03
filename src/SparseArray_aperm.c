/****************************************************************************
 ****************************************************************************
 **									   **
 **									   **
 **                     aperm() of a SparseArray object                    **
 **									   **
 **									   **
 **             2D case (transposition) & multidimensional case            **
 **									   **
 **									   **
 ****************************************************************************
 ****************************************************************************/
#include "SparseArray_aperm.h"

#include "Rvector_utils.h"
#include "leaf_utils.h"

#include <string.h>  /* for memset() */


static void *alloc_quick_out_nzvals_p(R_xlen_t n, SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case INTSXP: case LGLSXP: return R_alloc(n, sizeof(int *));
	    case REALSXP:             return R_alloc(n, sizeof(double *));
	    case CPLXSXP:             return R_alloc(n, sizeof(Rcomplex *));
	    case RAWSXP:              return R_alloc(n, sizeof(Rbyte *));
	    case STRSXP: case VECSXP: return R_alloc(n, sizeof(SEXP));
	}
	error("SparseArray internal error in alloc_quick_out_nzvals_p():\n"
	      "    unsupported SparseArray type: \"%s\"", type2char(Rtype));
	return NULL;  /* will never reach this */
}

static inline void *shift_quick_out_nzvals_p(const void **p,
		SEXPTYPE Rtype, R_xlen_t by)
{
	switch (Rtype) {
	    case INTSXP: case LGLSXP: return ((int      **) p) + by;
	    case REALSXP:             return ((double   **) p) + by;
	    case CPLXSXP:             return ((Rcomplex **) p) + by;
	    case RAWSXP:              return ((Rbyte    **) p) + by;
	    case STRSXP: case VECSXP: return ((SEXP      *) p) + by;
	}
	error("SparseArray internal error in shift_quick_out_nzvals_p():\n"
	      "    unsupported SparseArray type: \"%s\"", type2char(Rtype));
	return NULL;  /* will never reach this */
}

static inline void set_quick_out_nzvals_p(const void **quick_out_nzvals_p,
					  SEXPTYPE Rtype, SEXP nzvals)
{
	switch (Rtype) {
	    case INTSXP: case LGLSXP: {
		const int      **p = (const int **) quick_out_nzvals_p;
		*p = nzvals == R_NilValue ? NULL : INTEGER(nzvals);
		return;
	    }
	    case REALSXP: {
		const double   **p = (const double **) quick_out_nzvals_p;
		*p = nzvals == R_NilValue ? NULL : REAL(nzvals);
		return;
	    }
	    case CPLXSXP: {
		const Rcomplex **p = (const Rcomplex **) quick_out_nzvals_p;
		*p = nzvals == R_NilValue ? NULL : COMPLEX(nzvals);
		return;
	    }
	    case RAWSXP: {
		const Rbyte    **p = (const Rbyte **) quick_out_nzvals_p;
		*p = nzvals == R_NilValue ? NULL : RAW(nzvals);
		return;
	    }
	    case STRSXP: case VECSXP: {
		SEXP            *p = (SEXP *) quick_out_nzvals_p;
		*p = nzvals;
		return;
	    }
	}
	error("SparseArray internal error in set_quick_out_nzvals_p():\n"
	      "    unsupported SparseArray type: \"%s\"", type2char(Rtype));
}

static SEXP alloc_output_leaf(SEXPTYPE Rtype, int nzcount,
		const int *onecount_buf,
		const void **quick_out_nzvals_p, const int **quick_out_nzoffs_p)
{
	if (nzcount == 0)
		return R_NilValue;
	SEXP nzoffs = PROTECT(NEW_INTEGER(nzcount));
	*quick_out_nzoffs_p = INTEGER(nzoffs);
	SEXP nzvals;
	if (onecount_buf != NULL && *onecount_buf == nzcount) {
		/* Create a lacunar leaf. Will turn it into a standard leaf
		   later if LACUNAR_MODE_IS_ON is set to 0. */
		nzvals = R_NilValue;
	} else {
		nzvals = PROTECT(allocVector(Rtype, nzcount));
	}
	set_quick_out_nzvals_p(quick_out_nzvals_p, Rtype, nzvals);
	SEXP ans = PROTECT(zip_leaf(nzvals, nzoffs));
	UNPROTECT(nzvals == R_NilValue ? 2 : 3);
	return ans;
}

/* Inplace replacement!
   Maybe move this to SVT_SparseArray_class.c. */
static void REC_replace_lacunar_leaves_with_standard_leaves(
		SEXP SVT, const int *dim, int ndim, SEXPTYPE Rtype)
{
	if (SVT == R_NilValue)
		return;
	if (ndim == 1) {
		/* 'SVT' is a leaf (i.e. a 1D SVT). */
		SEXP nzvals, nzoffs;
		int nzcount = unzip_leaf(SVT, &nzvals, &nzoffs);
		if (nzvals != R_NilValue)  /* standard leaf */
			return;
		/* lacunar leaf */
		nzvals = PROTECT(_new_Rvector1(Rtype, nzcount));
		replace_leaf_nzvals(SVT, nzvals);
		UNPROTECT(1);
		return;
	}
	/* 'SVT' is a tree. */
	int SVT_len = LENGTH(SVT);  /* same as 'dim[ndim - 1]' */
	for (int i = 0; i < SVT_len; i++)
		REC_replace_lacunar_leaves_with_standard_leaves(
			VECTOR_ELT(SVT, i), dim, ndim - 1, Rtype);
	return;
}


/****************************************************************************
 *									    *
 *                         2D case (transposition)                          *
 *									    *
 ****************************************************************************/

/*
 * TODO: Do we still need this? C_aperm_SVT() below reimplements this for
 * the multidimensional case and with no significant/measurable differences
 * in terms of performance for the 2D case compared to C_transpose_2D_SVT().
 */

/* Collect number of nonzero elements (and ones) for each row in the input
   matrix (which will be the columns of the output matrix). */
static void collect_stats_on_input_rows(SEXP SVT, int nrow, int ncol,
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
			    _all_Rsubvec_elts_equal_one(nzvals, k, 1))
				onecount_buf[*nzoffs_p]++;
		}
	}
	return;
}

typedef void (*TransposeCol_FUNType)(int col_idx, SEXP col,
		void **quick_out_nzvals_p, int **quick_out_nzoffs_p,
		int *nzcount_buf);

/* Ignores 'nzcount_buf'. */
static void transpose_integer_col(int col_idx, SEXP col,
		void **quick_out_nzvals_p, int **quick_out_nzoffs_p,
		int *nzcount_buf)
{
	SEXP nzvals, nzoffs;
	int n = unzip_leaf(col, &nzvals, &nzoffs);
	const int *nzvals_p = NULL;  /* -Wmaybe-uninitialized */
	int v;
	if (nzvals != R_NilValue) {  /* standard leaf */
		nzvals_p = INTEGER(nzvals);
	} else {  /* lacunar leaf */
		v = int1;
	}
	const int *nzoffs_p = INTEGER(nzoffs);
	int **out_nzvals_p = (int **) quick_out_nzvals_p;
	for (int k = 0; k < n; k++) {
		int row_idx = *nzoffs_p;
		int **p = out_nzvals_p + row_idx;
		if (*p != NULL) {
			if (nzvals_p != NULL)
				v = nzvals_p[k];
			*((*p)++) = v;
		}
		*(quick_out_nzoffs_p[row_idx]++) = col_idx;
		nzoffs_p++;
	}
	return;
}

/* Ignores 'nzcount_buf'. */
static void transpose_double_col(int col_idx, SEXP col,
		void **quick_out_nzvals_p, int **quick_out_nzoffs_p,
		int *nzcount_buf)
{
	SEXP nzvals, nzoffs;
	int n = unzip_leaf(col, &nzvals, &nzoffs);
	const double *nzvals_p = NULL;  /* -Wmaybe-uninitialized */
	double v;
	if (nzvals != R_NilValue) {  /* standard leaf */
		nzvals_p = REAL(nzvals);
	} else {  /* lacunar leaf */
		v = double1;
	}
	const int *nzoffs_p = INTEGER(nzoffs);
	double **out_nzvals_p = (double **) quick_out_nzvals_p;
	for (int k = 0; k < n; k++) {
		int row_idx = *nzoffs_p;
		double **p = out_nzvals_p + row_idx;
		if (*p != NULL) {
			if (nzvals_p != NULL)
				v = nzvals_p[k];
			*((*p)++) = v;
		}
		*(quick_out_nzoffs_p[row_idx]++) = col_idx;
		nzoffs_p++;
	}
	return;
}

/* Ignores 'nzcount_buf'. */
static void transpose_complex_col(int col_idx, SEXP col,
		void **quick_out_nzvals_p, int **quick_out_nzoffs_p,
		int *nzcount_buf)
{
	SEXP nzvals, nzoffs;
	int n = unzip_leaf(col, &nzvals, &nzoffs);
	const Rcomplex *nzvals_p = NULL;  /* -Wmaybe-uninitialized */
	Rcomplex v;
	if (nzvals != R_NilValue) {  /* standard leaf */
		nzvals_p = COMPLEX(nzvals);
	} else {
		v = Rcomplex1;  /* lacunar leaf */
	}
	const int *nzoffs_p = INTEGER(nzoffs);
	Rcomplex **out_nzvals_p = (Rcomplex **) quick_out_nzvals_p;
	for (int k = 0; k < n; k++) {
		int row_idx = *nzoffs_p;
		Rcomplex **p = out_nzvals_p + row_idx;
		if (*p != NULL) {
			if (nzvals_p != NULL)
				v = nzvals_p[k];
			*((*p)++) = v;
		}
		*(quick_out_nzoffs_p[row_idx]++) = col_idx;
		nzoffs_p++;
	}
	return;
}

/* Ignores 'nzcount_buf'. */
static void transpose_raw_col(int col_idx, SEXP col,
		void **quick_out_nzvals_p, int **quick_out_nzoffs_p,
		int *nzcount_buf)
{
	SEXP nzvals, nzoffs;
	int n = unzip_leaf(col, &nzvals, &nzoffs);
	const Rbyte *nzvals_p = NULL;  /* -Wmaybe-uninitialized */
	Rbyte v;
	if (nzvals != R_NilValue) {  /* standard leaf */
		nzvals_p = RAW(nzvals);
	} else {
		v = Rbyte1;  /* lacunar leaf */
	}
	const int *nzoffs_p = INTEGER(nzoffs);
	Rbyte **out_nzvals_p = (Rbyte **) quick_out_nzvals_p;
	for (int k = 0; k < n; k++) {
		int row_idx = *nzoffs_p;
		Rbyte **p = out_nzvals_p + row_idx;
		if (*p != NULL) {
			if (nzvals_p != NULL)
				v = nzvals_p[k];
			*((*p)++) = v;
		}
		*(quick_out_nzoffs_p[row_idx]++) = col_idx;
		nzoffs_p++;
	}
	return;
}

static void transpose_character_col(int col_idx, SEXP col,
		void **quick_out_nzvals_p, int **quick_out_nzoffs_p,
		int *nzcount_buf)
{
	SEXP nzvals, nzoffs;
	int n = unzip_leaf(col, &nzvals, &nzoffs);
	const int *nzoffs_p = INTEGER(nzoffs);
	SEXP *out_nzvals_p = (SEXP *) quick_out_nzvals_p;
	for (int k = 0; k < n; k++) {
		int row_idx = *nzoffs_p;
		int k2 = nzcount_buf[row_idx]++;
		SET_STRING_ELT(out_nzvals_p[row_idx], k2,
			       STRING_ELT(nzvals, k));
		*(quick_out_nzoffs_p[row_idx]++) = col_idx;
		nzoffs_p++;
	}
	return;
}

static void transpose_list_col(int col_idx, SEXP col,
		void **quick_out_nzvals_p, int **quick_out_nzoffs_p,
		int *nzcount_buf)
{
	SEXP nzvals, nzoffs;
	int n = unzip_leaf(col, &nzvals, &nzoffs);
	const int *nzoffs_p = INTEGER(nzoffs);
	SEXP *out_nzvals_p = (SEXP *) quick_out_nzvals_p;
	for (int k = 0; k < n; k++) {
		int row_idx = *nzoffs_p;
		int k2 = nzcount_buf[row_idx]++;
		SET_VECTOR_ELT(out_nzvals_p[row_idx], k2,
			       VECTOR_ELT(nzvals, k));
		*(quick_out_nzoffs_p[row_idx]++) = col_idx;
		nzoffs_p++;
	}
	return;
}

static TransposeCol_FUNType select_transpose_col_FUN(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case INTSXP: case LGLSXP: return transpose_integer_col;
	    case REALSXP:             return transpose_double_col;
	    case CPLXSXP:             return transpose_complex_col;
	    case RAWSXP:              return transpose_raw_col;
	    case STRSXP:              return transpose_character_col;
	    case VECSXP:              return transpose_list_col;
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

	/* 1st pass */
	collect_stats_on_input_rows(SVT, nrow, ncol, nzcount_buf, onecount_buf);

	/* 2nd pass: Allocate transposed SVT. */
	void **quick_out_nzvals_p =
			alloc_quick_out_nzvals_p((R_xlen_t) nrow, Rtype);
	int **quick_out_nzoffs_p = (int **) R_alloc(nrow, sizeof(int *));
	SEXP ans = PROTECT(NEW_LIST(nrow));
	for (int i = 0; i < nrow; i++) {
		const void **p = shift_quick_out_nzvals_p(
					(const void **) quick_out_nzvals_p,
					Rtype, (R_xlen_t) i);
		SEXP ans_elt = alloc_output_leaf(Rtype, nzcount_buf[i],
					onecount_buf + i,
					p,
					(const int **) quick_out_nzoffs_p + i);
		if (ans_elt != R_NilValue) {
			PROTECT(ans_elt);
			SET_VECTOR_ELT(ans, i, ans_elt);
			UNPROTECT(1);
		}
	}

	/* 3rd pass: Fill leaves of transposed SVT. */
	memset(nzcount_buf, 0, sizeof(int) * nrow);
	for (int j = 0; j < ncol; j++) {
		SEXP leaf = VECTOR_ELT(SVT, j);
		if (leaf == R_NilValue)
			continue;
		transpose_col_FUN(j, leaf,
			quick_out_nzvals_p, quick_out_nzoffs_p,
			nzcount_buf);
	}

	/* 4th pass */
	if (onecount_buf != NULL && LACUNAR_MODE_IS_ON == 0) {
		const int ans_dim[2] = {ncol, nrow};
		REC_replace_lacunar_leaves_with_standard_leaves(
			ans, ans_dim, 2, Rtype);
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
 *									    *
 *                         Multidimensional aperm()                         *
 *									    *
 ****************************************************************************/


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
 * Buffers used by aperm_SVT_shattering_leaves()
 */

typedef struct aperm0_bufs_t {
	/* 'nzcount_buf' must be an array of dimensions 'dim(out_array)[-1]'
	   where 'out_array' is the array produced by permuting the
	   dimensions of the input array. Note that 'dim(out_array)'
	   will be 'dim(in_array)[perm]'. */
	int *nzcount_buf;
	int *onecount_buf;  /* same as 'nzcount_buf' or NULL */
	R_xlen_t nzcount_buf_len;  /* prod(dim(out_array)[-1]) */
	R_xlen_t *nzcount_buf_incs;
	R_xlen_t *outer_incs;
	void **quick_out_nzvals_p;
	int **quick_out_nzoffs_p;
} Aperm0Bufs;

/* Will always be called with 'ndim' >= 2. */
static void init_A0Bufs(Aperm0Bufs *A0Bufs, const int *dim, int ndim,
			const int *perm, SEXPTYPE Rtype)
{
	R_xlen_t *nzcount_buf_incs = (R_xlen_t *)
		R_alloc(ndim - 1, sizeof(R_xlen_t));
	R_xlen_t *outer_incs = (R_xlen_t *)
		R_alloc(ndim, sizeof(R_xlen_t));
	outer_incs[perm[0] - 1] = 0;
	R_xlen_t nzcount_buf_len = 1;
	for (int along = 1; along < ndim; along++) {
		int p = perm[along] - 1;
		nzcount_buf_incs[along - 1] = nzcount_buf_len;
		outer_incs[p] = nzcount_buf_len;
		nzcount_buf_len *= dim[p];
	}
/*
	printf("nzcount_buf_len = %lu\n", nzcount_buf_len);
	for (int along = 0; along < ndim - 1; along++)
		printf("nzcount_buf_incs[%d] = %lu\n",
			along, nzcount_buf_incs[along]);
	for (int along = 0; along < ndim; along++)
		printf("outer_incs[%d] = %lu\n",
			along, outer_incs[along]);
*/
	A0Bufs->nzcount_buf = (int *)
		R_alloc(nzcount_buf_len, sizeof(int));
	A0Bufs->onecount_buf = NULL;
	if (Rtype != STRSXP && Rtype != VECSXP)
		A0Bufs->onecount_buf = (int *)
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
   The former is used by aperm_SVT_preserving_leaves() (the easy case) and
   the latter are used by aperm_SVT_shattering_leaves() (the hard case). */
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
 * EASY CASE: aperm_SVT_preserving_leaves()
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
static void REC_aperm_SVT_preserving_leaves(
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
			REC_aperm_SVT_preserving_leaves(
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

static SEXP aperm_SVT_preserving_leaves(
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
	REC_aperm_SVT_preserving_leaves(SVT, dim, ndim, perm,
					    ans_dim, ndim,
					    inner_margin, coords0_buf, ans);
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * HARD CASE: aperm_SVT_shattering_leaves()
 */

static inline void scan_input_leaf(SEXP leaf,
		R_xlen_t outer_inc0, R_xlen_t outer_offset0,
		int *nzcount_buf, int *onecount_buf)
{
	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
	const int *nzoffs_p = INTEGER(nzoffs);
	for (int k = 0; k < nzcount; k++, nzoffs_p++) {
		R_xlen_t outer_idx = outer_offset0 + outer_inc0 * *nzoffs_p;
		nzcount_buf[outer_idx]++;
		if (onecount_buf == NULL)
			continue;
		if (nzvals == R_NilValue ||  /* lacunar leaf */
		    _all_Rsubvec_elts_equal_one(nzvals, k, 1))
			onecount_buf[outer_idx]++;
	}
	return;
}

/* Recursive. */
static void REC_collect_stats_on_output_leaves(SEXP SVT, int ndim,
		const R_xlen_t *outer_incs, R_xlen_t outer_offset0,
		int *nzcount_buf, int *onecount_buf)
{
	if (SVT == R_NilValue)
		return;
	R_xlen_t outer_inc = outer_incs[ndim - 1];
	if (ndim == 1) {
		scan_input_leaf(SVT, outer_inc, outer_offset0,
				nzcount_buf, onecount_buf);
		return;
	}
	int SVT_len = LENGTH(SVT);
	for (int i = 0; i < SVT_len; i++, outer_offset0 += outer_inc) {
		SEXP subSVT = VECTOR_ELT(SVT, i);
		REC_collect_stats_on_output_leaves(subSVT, ndim - 1,
						   outer_incs, outer_offset0,
						   nzcount_buf, onecount_buf);
	}
	return;
}

/* Recursive. */
static SEXP REC_grow_output_tree(const int *dim, int ndim,
		SEXPTYPE Rtype,
		const R_xlen_t *nzcount_buf_incs,
		const int *nzcount_buf, const int *onecount_buf,
		const void **quick_out_nzvals_p,
		const int **quick_out_nzoffs_p)
{
	if (ndim == 1)
		return alloc_output_leaf(Rtype, *nzcount_buf,
					 onecount_buf,
					 quick_out_nzvals_p,
					 quick_out_nzoffs_p);
	int ans_len = dim[ndim - 1];
	R_xlen_t buf_inc = nzcount_buf_incs[ndim - 2];
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	int is_empty = 1;
	for (int i = 0; i < ans_len; i++) {
		SEXP ans_elt = REC_grow_output_tree(dim, ndim - 1,
					Rtype, nzcount_buf_incs,
					nzcount_buf, onecount_buf,
					quick_out_nzvals_p, quick_out_nzoffs_p);
		if (ans_elt != R_NilValue) {
			PROTECT(ans_elt);
			SET_VECTOR_ELT(ans, i, ans_elt);
			UNPROTECT(1);
			is_empty = 0;
		}
		nzcount_buf += buf_inc;
		if (onecount_buf != NULL)
			onecount_buf += buf_inc;
		quick_out_nzoffs_p += buf_inc;
		quick_out_nzvals_p = shift_quick_out_nzvals_p(
						quick_out_nzvals_p,
						Rtype, buf_inc);
	}
	UNPROTECT(1);
	return is_empty ? R_NilValue : ans;
}

#define	FUNDEF_spray_leaf(type)						\
	(SEXP leaf, int inner_idx,					\
	 R_xlen_t outer_inc0, R_xlen_t outer_offset0,			\
	 int *nzcount_buf,						\
	 void *quick_out_nzvals_p, int **quick_out_nzoffs_p)		\
{									\
	SEXP nzvals, nzoffs;						\
	int n = unzip_leaf(leaf, &nzvals, &nzoffs);			\
	const type *nzvals_p = NULL;  /* -Wmaybe-uninitialized */	\
	type v;								\
	if (nzvals != R_NilValue) {  /* standard leaf */		\
		nzvals_p = (const type *) DATAPTR(nzvals);		\
	} else {  /* lacunar leaf */					\
		v = type ## 1;						\
	}								\
	const int *nzoffs_p = INTEGER(nzoffs);				\
	type **out_nzvals_p = (type **) quick_out_nzvals_p;		\
	for (int k = 0; k < n; k++) {					\
		R_xlen_t outer_idx = outer_offset0 +			\
				     outer_inc0 * *nzoffs_p;		\
		int k2 = nzcount_buf[outer_idx]++;			\
		type *p = out_nzvals_p[outer_idx];			\
		if (p != NULL) {					\
			if (nzvals_p != NULL)				\
				v = nzvals_p[k];			\
			p[k2] = v;					\
		}							\
		quick_out_nzoffs_p[outer_idx][k2] = inner_idx;		\
		nzoffs_p++;						\
	}								\
	return;								\
}

static inline void spray_integer_leaf FUNDEF_spray_leaf(int)
static inline void spray_double_leaf FUNDEF_spray_leaf(double)
static inline void spray_complex_leaf FUNDEF_spray_leaf(Rcomplex)
static inline void spray_raw_leaf FUNDEF_spray_leaf(Rbyte)

static inline void spray_character_leaf(SEXP leaf, int inner_idx,
		R_xlen_t outer_inc0, R_xlen_t outer_offset0,
		int *nzcount_buf,
		void *quick_out_nzvals_p, int **quick_out_nzoffs_p)
{
	SEXP nzvals, nzoffs;
	int n = unzip_leaf(leaf, &nzvals, &nzoffs);
	const int *nzoffs_p = INTEGER(nzoffs);
	SEXP *out_nzvals_p = (SEXP *) quick_out_nzvals_p;
	for (int k = 0; k < n; k++) {
		R_xlen_t outer_idx = outer_offset0 + outer_inc0 * *nzoffs_p;
		int k2 = nzcount_buf[outer_idx]++;
		SET_STRING_ELT(out_nzvals_p[outer_idx], k2,
			       STRING_ELT(nzvals, k));
		quick_out_nzoffs_p[outer_idx][k2] = inner_idx;
		nzoffs_p++;
	}
	return;
}

static inline void spray_list_leaf(SEXP leaf, int inner_idx,
		R_xlen_t outer_inc0, R_xlen_t outer_offset0,
		int *nzcount_buf,
		void *quick_out_nzvals_p, int **quick_out_nzoffs_p)
{
	SEXP nzvals, nzoffs;
	int n = unzip_leaf(leaf, &nzvals, &nzoffs);
	const int *nzoffs_p = INTEGER(nzoffs);
	SEXP *out_nzvals_p = (SEXP *) quick_out_nzvals_p;
	for (int k = 0; k < n; k++) {
		R_xlen_t outer_idx = outer_offset0 + outer_inc0 * *nzoffs_p;
		int k2 = nzcount_buf[outer_idx]++;
		SET_VECTOR_ELT(out_nzvals_p[outer_idx], k2,
			       VECTOR_ELT(nzvals, k));
		quick_out_nzoffs_p[outer_idx][k2] = inner_idx;
		nzoffs_p++;
	}
	return;
}

static void spray_input_leaf_on_output_leaves(SEXP leaf, SEXPTYPE Rtype,
		R_xlen_t outer_inc0, R_xlen_t outer_offset0,
		int *nzcount_buf,
		void *quick_out_nzvals_p, int **quick_out_nzoffs_p,
		int inner_idx)
{
	void (*spray_FUN)(SEXP leaf, int inner_idx,
			  R_xlen_t outer_inc0, R_xlen_t outer_offset0,
			  int *nzcount_buf,
			  void *quick_out_nzvals_p, int **quick_out_nzoffs_p);
	switch (Rtype) {
	    case INTSXP: case LGLSXP: spray_FUN = spray_integer_leaf;   break;
	    case REALSXP:             spray_FUN = spray_double_leaf;    break;
	    case CPLXSXP:             spray_FUN = spray_complex_leaf;   break;
	    case RAWSXP:              spray_FUN = spray_raw_leaf;       break;
	    case STRSXP:              spray_FUN = spray_character_leaf; break;
	    case VECSXP:              spray_FUN = spray_list_leaf;      break;
	    default:
		error("SparseArray internal error in "
		      "spray_input_leaf_on_output_leaves():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));
	}
	spray_FUN(leaf, inner_idx,
		  outer_inc0, outer_offset0,
		  nzcount_buf,
		  quick_out_nzvals_p, quick_out_nzoffs_p);
	return;
}

/* Recursive. */
static void REC_spray_input_leaves_on_output_leaves(
		SEXP SVT, int ndim, SEXPTYPE Rtype,
		const R_xlen_t *outer_incs, R_xlen_t outer_offset0,
		int *nzcount_buf,
		void *quick_out_nzvals_p, int **quick_out_nzoffs_p)
{
	static int inner_idx;

	if (SVT == R_NilValue)
		return;
	R_xlen_t outer_inc = outer_incs[ndim - 1];
	if (ndim == 1) {
		spray_input_leaf_on_output_leaves(SVT, Rtype,
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
		REC_spray_input_leaves_on_output_leaves(
					subSVT, ndim - 1, Rtype,
					outer_incs, outer_offset0,
					nzcount_buf,
					quick_out_nzvals_p, quick_out_nzoffs_p);
	}
	return;
}

static SEXP aperm_SVT_shattering_leaves(
		SEXP SVT, const int *dim, int ndim, SEXPTYPE Rtype,
		const int *ans_dim, Aperm0Bufs *A0Bufs)
{
	SEXP ans;

	/* 1st pass */
	memset(A0Bufs->nzcount_buf, 0, sizeof(int) * A0Bufs->nzcount_buf_len);
	if (A0Bufs->onecount_buf != NULL)
		memset(A0Bufs->onecount_buf, 0,
		       sizeof(int) * A0Bufs->nzcount_buf_len);
	REC_collect_stats_on_output_leaves(SVT, ndim,
				A0Bufs->outer_incs, 0,
				A0Bufs->nzcount_buf,
				A0Bufs->onecount_buf);
	/* 2nd pass */
	ans = PROTECT(REC_grow_output_tree(ans_dim, ndim, Rtype,
				A0Bufs->nzcount_buf_incs,
				A0Bufs->nzcount_buf,
				A0Bufs->onecount_buf,
				(const void **) A0Bufs->quick_out_nzvals_p,
				(const int **) A0Bufs->quick_out_nzoffs_p));
	/* 3rd pass */
	memset(A0Bufs->nzcount_buf, 0, sizeof(int) * A0Bufs->nzcount_buf_len);
	REC_spray_input_leaves_on_output_leaves(SVT, ndim, Rtype,
				A0Bufs->outer_incs, 0,
				A0Bufs->nzcount_buf,
				A0Bufs->quick_out_nzvals_p,
				A0Bufs->quick_out_nzoffs_p);

	/* 4th pass */
	if (A0Bufs->onecount_buf != NULL && LACUNAR_MODE_IS_ON == 0)
		REC_replace_lacunar_leaves_with_standard_leaves(
			ans, ans_dim, ndim, Rtype);

	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_aperm0_SVT()
 */

/* Can handle any array permutation. However, it is optimized for the "no
   outer margin" case, that is, for the case when the permutation vector
   has no outer margin. When the permutation vector has a nonzero outer
   margin, aperm0_SVT() won't be as efficient as C_aperm_SVT().
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
		return aperm_SVT_preserving_leaves(SVT, dim, ndim,
						   perm, ans_dim,
						   inner_margin, coords0_buf);
	}

	/* HARD CASE: Leaves won't be preserved i.e. the permuted array will
	   use entirely new leaves. It is much less efficient than the EASY
	   CASE because memory needs to be allocated for the new leaves. */
	return aperm_SVT_shattering_leaves(SVT, dim, ndim, Rtype,
					   ans_dim, A0Bufs);
}

/* --- .Call ENTRY POINT ---
 * Just a wrapper around aperm0_SVT(). Provided only for convenience testing
 * of aperm0_SVT() from the R command line.
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

