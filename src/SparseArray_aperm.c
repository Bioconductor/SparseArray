/****************************************************************************
 *                  Transposition of a SparseArray object                   *
 ****************************************************************************/
#include "SparseArray_aperm.h"

#include "Rvector_utils.h"
#include "leaf_vector_utils.h"
#include "SBT_utils.h"


/****************************************************************************
 * C_transpose_2D_SVT()
 *
 * TODO: Do we still need this? C_transpose_SVT() below reimplements this for
 * the multidimensional case and with no significant/measurable differences
 * from a performance point of view with C_transpose_2D_SVT() in the 2D case.
 */

static void count_nonzero_vals_per_row(SEXP SVT, int nrow, int ncol,
		int *nzcounts)
{
	int j, lv_len, k;
	SEXP subSVT, lv_offs, lv_vals;
	const int *p;

	memset(nzcounts, 0, sizeof(int) * nrow);
	for (j = 0; j < ncol; j++) {
		subSVT = VECTOR_ELT(SVT, j);
		if (subSVT == R_NilValue)
			continue;
		/* 'subSVT' is a "leaf vector". */
		lv_len = _split_leaf_vector(subSVT, &lv_offs, &lv_vals);
		if (lv_len < 0)
			error("SparseArray internal error in "
			      "count_nonzero_vals_per_row():\n"
			      "    invalid SVT_SparseMatrix object");
		for (k = 0, p = INTEGER(lv_offs); k < lv_len; k++, p++)
			nzcounts[*p]++;
	}
	return;
}

static void **set_quick_out_vals_p(SEXP out_SVT, SEXPTYPE Rtype)
{
	int out_SVT_len, i;
	SEXP lv;

	out_SVT_len = LENGTH(out_SVT);
	switch (Rtype) {
	    case LGLSXP: case INTSXP: {
		int **vals_p, **p;
		vals_p = (int **) R_alloc(out_SVT_len, sizeof(int *));
		for (i = 0, p = vals_p; i < out_SVT_len; i++, p++) {
			lv = VECTOR_ELT(out_SVT, i);
			if (lv != R_NilValue)
				*p = INTEGER(VECTOR_ELT(lv, 1));
		}
		return (void **) vals_p;
	    }
	    case REALSXP: {
		double **vals_p, **p;
		vals_p = (double **) R_alloc(out_SVT_len, sizeof(double *));
		for (i = 0, p = vals_p; i < out_SVT_len; i++, p++) {
			lv = VECTOR_ELT(out_SVT, i);
			if (lv != R_NilValue)
				*p = REAL(VECTOR_ELT(lv, 1));
		}
		return (void **) vals_p;
	    }
	    case CPLXSXP: {
		Rcomplex **vals_p, **p;
		vals_p = (Rcomplex **) R_alloc(out_SVT_len, sizeof(Rcomplex *));
		for (i = 0, p = vals_p; i < out_SVT_len; i++, p++) {
			lv = VECTOR_ELT(out_SVT, i);
			if (lv != R_NilValue)
				*p = COMPLEX(VECTOR_ELT(lv, 1));
		}
		return (void **) vals_p;
	    }
	    case RAWSXP: {
		Rbyte **vals_p, **p;
		vals_p = (Rbyte **) R_alloc(out_SVT_len, sizeof(Rbyte *));
		for (i = 0, p = vals_p; i < out_SVT_len; i++, p++) {
			lv = VECTOR_ELT(out_SVT, i);
			if (lv != R_NilValue)
				*p = RAW(VECTOR_ELT(lv, 1));
		}
		return (void **) vals_p;
	    }
	}
	/* STRSXP and VECSXP cases. */
	return NULL;
}

typedef void (*TransposeCol_FUNType)(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcount_buf);

/* Ignores 'out_SVT' and 'nzcount_buf'. */
static void transpose_INTEGER_col(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcount_buf)
{
	int **vals_p;
	int lv_len, k, row_idx;
	const int *v;

	vals_p = (int **) quick_out_vals_p;
	lv_len = LENGTH(lv_vals);
	for (k = 0, v = INTEGER(lv_vals); k < lv_len; k++, v++) {
		row_idx = *offs;
		*(quick_out_offs_p[row_idx]++) = col_idx;
		*(vals_p[row_idx]++) = *v;
		offs++;
	}
	return;
}

/* Ignores 'out_SVT' and 'nzcount_buf'. */
static void transpose_NUMERIC_col(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcount_buf)
{
	double **vals_p;
	int lv_len, k, row_idx;
	const double *v;

	vals_p = (double **) quick_out_vals_p;
	lv_len = LENGTH(lv_vals);
	for (k = 0, v = REAL(lv_vals); k < lv_len; k++, v++) {
		row_idx = *offs;
		*(quick_out_offs_p[row_idx]++) = col_idx;
		*(vals_p[row_idx]++) = *v;
		offs++;
	}
	return;
}

/* Ignores 'out_SVT' and 'nzcount_buf'. */
static void transpose_COMPLEX_col(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcount_buf)
{
	Rcomplex **vals_p;
	int lv_len, k, row_idx;
	const Rcomplex *v;

	vals_p = (Rcomplex **) quick_out_vals_p;
	lv_len = LENGTH(lv_vals);
	for (k = 0, v = COMPLEX(lv_vals); k < lv_len; k++, v++) {
		row_idx = *offs;
		*(quick_out_offs_p[row_idx]++) = col_idx;
		*(vals_p[row_idx]++) = *v;
		offs++;
	}
	return;
}

/* Ignores 'out_SVT' and 'nzcount_buf'. */
static void transpose_RAW_col(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcount_buf)
{
	Rbyte **vals_p;
	int lv_len, k, row_idx;
	const Rbyte *v;

	vals_p = (Rbyte **) quick_out_vals_p;
	lv_len = LENGTH(lv_vals);
	for (k = 0, v = RAW(lv_vals); k < lv_len; k++, v++) {
		row_idx = *offs;
		*(quick_out_offs_p[row_idx]++) = col_idx;
		*(vals_p[row_idx]++) = *v;
		offs++;
	}
	return;
}

/* Ignores 'quick_out_vals_p'. */
static void transpose_CHARACTER_col(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcount_buf)
{
	int lv_len, k, row_idx;
	SEXP out_lv;

	lv_len = LENGTH(lv_vals);
	for (k = 0; k < lv_len; k++) {
		row_idx = *offs;
		*(quick_out_offs_p[row_idx]++) = col_idx;
		out_lv = VECTOR_ELT(out_SVT, row_idx);
		_copy_CHARACTER_elt(lv_vals, (R_xlen_t) k,
			VECTOR_ELT(out_lv, 1),
			(R_xlen_t) nzcount_buf[row_idx]++);
		offs++;
	}
	return;
}

/* Ignores 'quick_out_vals_p'. */
static void transpose_LIST_col(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcount_buf)
{
	int lv_len, k, row_idx;
	SEXP out_lv;

	lv_len = LENGTH(lv_vals);
	for (k = 0; k < lv_len; k++) {
		row_idx = *offs;
		*(quick_out_offs_p[row_idx]++) = col_idx;
		out_lv = VECTOR_ELT(out_SVT, row_idx);
		_copy_LIST_elt(lv_vals, (R_xlen_t) k,
			VECTOR_ELT(out_lv, 1),
			(R_xlen_t) nzcount_buf[row_idx]++);
		offs++;
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
		int *nzcount_buf)
{
	TransposeCol_FUNType transpose_col_FUN;
	SEXP ans, ans_elt, subSVT, lv_offs, lv_vals;
	int **quick_out_offs_p;
	void **quick_out_vals_p;
	int i, j, lv_len;

	/* 1st pass: Count the number of nonzero values per row in the
	   input object. */
	count_nonzero_vals_per_row(SVT, nrow, ncol, nzcount_buf);

	/* 2nd pass: Build the transposed SVT. */
	transpose_col_FUN = select_transpose_col_FUN(Rtype);
	if (transpose_col_FUN == NULL)
		error("SparseArray internal error in "
		      "transpose_2D_SVT():\n"
		      "    SVT_SparseMatrix object has invalid type");

	ans = PROTECT(NEW_LIST(nrow));
	quick_out_offs_p = (int **) R_alloc(nrow, sizeof(int *));
	for (i = 0; i < nrow; i++) {
		lv_len = nzcount_buf[i];
		if (lv_len != 0) {
			ans_elt = PROTECT(_alloc_leaf_vector(lv_len, Rtype));
			SET_VECTOR_ELT(ans, i, ans_elt);
			UNPROTECT(1);
			quick_out_offs_p[i] = INTEGER(VECTOR_ELT(ans_elt, 0));
		}
	}
	quick_out_vals_p = set_quick_out_vals_p(ans, Rtype);

	memset(nzcount_buf, 0, sizeof(int) * nrow);
	for (j = 0; j < ncol; j++) {
		subSVT = VECTOR_ELT(SVT, j);
		if (subSVT == R_NilValue)
			continue;
		/* 'subSVT' is a "leaf vector". */
		lv_len = _split_leaf_vector(subSVT, &lv_offs, &lv_vals);
		if (lv_len < 0) {
			UNPROTECT(1);
			error("SparseArray internal error in "
			      "transpose_2D_SVT():\n"
			      "    invalid SVT_SparseMatrix object");
		}
		transpose_col_FUN(j,
			INTEGER(lv_offs), lv_vals,
			quick_out_offs_p, quick_out_vals_p,
			ans, nzcount_buf);
	}
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_transpose_2D_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT)
{
	SEXPTYPE x_Rtype;
	int x_nrow, x_ncol, *nzcount_buf;

	x_Rtype = _get_Rtype_from_Rstring(x_type);
	if (x_Rtype == 0)
		error("SparseArray internal error in "
		      "C_transpose_2D_SVT():\n"
		      "    SVT_SparseMatrix object has invalid type");

	if (LENGTH(x_dim) != 2)
		error("object to transpose must have exactly 2 dimensions");

	if (x_SVT == R_NilValue)
		return x_SVT;

	x_nrow = INTEGER(x_dim)[0];
	x_ncol = INTEGER(x_dim)[1];
	nzcount_buf = (int *) R_alloc(x_nrow, sizeof(int));

	return transpose_2D_SVT(x_SVT, x_nrow, x_ncol, x_Rtype, nzcount_buf);
}


/****************************************************************************
 * Basic manipulation of the transpose_SVT() buffers
 */

static void *alloc_quick_out_vals_p(unsigned long long int n, SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return R_alloc(n, sizeof(int *));
	    case REALSXP:             return R_alloc(n, sizeof(double *));
	    case CPLXSXP:             return R_alloc(n, sizeof(Rcomplex *));
	    case RAWSXP:              return R_alloc(n, sizeof(Rbyte *));
	    case STRSXP: case VECSXP: return R_alloc(n, sizeof(SEXP));
	}
	error("SparseArray internal error in alloc_quick_out_vals_p():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return NULL;  /* will never reach this */
}

static inline void init_quick_out_vals_p(void *p, SEXPTYPE Rtype, SEXP vals)
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

static inline void *inc_quick_out_vals_p(void *p, SEXPTYPE Rtype,
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

typedef struct transpose_bufs_t {
	int *nzcount_buf;     /* array of dims 'rev(head(dim(x), n=-1))' */
	unsigned long long int nzcount_buf_len;
	unsigned long long int *nzcount_buf_incs;
	unsigned long long int *outer_incs;
	int **quick_out_offs_p;
	void *quick_out_vals_p;
} TransposeBufs;

static TransposeBufs new_TBufs(const int *dim, int ndim, SEXPTYPE Rtype)
{
	TransposeBufs TBufs;

	unsigned long long int *nzcount_buf_incs;
	unsigned long long int *outer_incs;
	unsigned long long int nzcount_buf_len;
	int along, d;

	nzcount_buf_incs = (unsigned long long int *)
		R_alloc(ndim - 1, sizeof(unsigned long long int));
	outer_incs = (unsigned long long int *)
		R_alloc(ndim, sizeof(unsigned long long int));
	outer_incs[ndim - 1] = 0;
	nzcount_buf_len = 1;
	for (along = ndim - 2; along >= 0; along--) {
		d = dim[along];
		nzcount_buf_incs[ndim - 2 - along] = nzcount_buf_len;
		outer_incs[along] = nzcount_buf_len;
		nzcount_buf_len *= d;
	}

	TBufs.nzcount_buf = (int *)
		R_alloc(nzcount_buf_len, sizeof(int));
	TBufs.nzcount_buf_len = nzcount_buf_len;
	TBufs.nzcount_buf_incs = nzcount_buf_incs;
	TBufs.outer_incs = outer_incs;
	TBufs.quick_out_offs_p = (int **)
		R_alloc(nzcount_buf_len, sizeof(int *));
	TBufs.quick_out_vals_p =
		alloc_quick_out_vals_p(nzcount_buf_len, Rtype);
	return TBufs;
}


/****************************************************************************
 * Multidimensional transposition: C_transpose_SVT()
 */

static inline void count_lv_nzvals(SEXP lv,
		unsigned long long int outer_inc0,
		unsigned long long int outer_offset0,
		int *nzcount_buf)
{
	int lv_len, k;
	SEXP lv_offs, lv_vals;
	const int *lv_offs_p;
	unsigned long long int outer_idx;

	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	if (lv_len < 0)
		error("SparseArray internal error in "
		      "count_lv_nzvals():\n"
		      "    invalid leaf vector");
	lv_offs_p = INTEGER(lv_offs);
	for (k = 0; k < lv_len; k++) {
		outer_idx = outer_offset0 + outer_inc0 * *lv_offs_p;
		nzcount_buf[outer_idx]++;
		lv_offs_p++;
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
		count_lv_nzvals(SVT, outer_inc, outer_offset0, nzcount_buf);
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
static SEXP REC_grow_and_alloc_leaves(const int *dim, int ndim,
		SEXPTYPE Rtype,
		const unsigned long long int *nzcount_buf_incs,
		const int *nzcount_buf,
		int **quick_out_offs_p, void *quick_out_vals_p)
{
	int ans_len, is_empty, i;
	unsigned long long int buf_inc;
	SEXP ans, lv_offs, lv_vals, ans_elt;

	if (ndim == 1) {
		if (*nzcount_buf == 0)
			return R_NilValue;
		ans = PROTECT(
			_alloc_and_split_leaf_vector(*nzcount_buf,
						     Rtype, &lv_offs, &lv_vals)
		);
		*quick_out_offs_p = INTEGER(lv_offs);
		init_quick_out_vals_p(quick_out_vals_p, Rtype, lv_vals);
		UNPROTECT(1);
		return ans;
	}
	ans_len = dim[ndim - 1];
	buf_inc = nzcount_buf_incs[ndim - 2];
	ans = PROTECT(NEW_LIST(ans_len));
	is_empty = 1;
	for (i = 0; i < ans_len; i++) {
		ans_elt = REC_grow_and_alloc_leaves(dim, ndim - 1, Rtype,
					nzcount_buf_incs, nzcount_buf,
					quick_out_offs_p, quick_out_vals_p);
		if (ans_elt != R_NilValue) {
			PROTECT(ans_elt);
			SET_VECTOR_ELT(ans, i, ans_elt);
			UNPROTECT(1);
			is_empty = 0;
		}
		nzcount_buf += buf_inc;
		quick_out_offs_p += buf_inc;
		quick_out_vals_p = inc_quick_out_vals_p(quick_out_vals_p,
							Rtype, buf_inc);
	}
	UNPROTECT(1);
	return is_empty ? R_NilValue : ans;
}

#define	FUNDEF_spray_ans_with_vals(type)				    \
	(const int *offs, SEXP vals,					    \
	 unsigned long long int outer_inc0,				    \
	 unsigned long long int outer_offset0,				    \
	 int *nzcount_buf, int **quick_out_offs_p, void *quick_out_vals_p,  \
	 int inner_idx)							    \
{									    \
	int n, k, k2;							    \
	const type *vals_p;						    \
	type **out_vals_p;						    \
	unsigned long long int outer_idx;				    \
									    \
	n = LENGTH(vals);						    \
	vals_p = DATAPTR(vals);						    \
	out_vals_p = quick_out_vals_p;					    \
	for (k = 0; k < n; k++) {					    \
		outer_idx = outer_offset0 + outer_inc0 * *offs;		    \
		k2 = nzcount_buf[outer_idx]++;				    \
		quick_out_offs_p[outer_idx][k2] = inner_idx;		    \
		out_vals_p[outer_idx][k2] = *vals_p;			    \
		offs++;							    \
		vals_p++;						    \
	}								    \
	return;								    \
}

static inline void spray_ans_with_ints
	FUNDEF_spray_ans_with_vals(int)
static inline void spray_ans_with_doubles
	FUNDEF_spray_ans_with_vals(double)
static inline void spray_ans_with_Rcomplexes
	FUNDEF_spray_ans_with_vals(Rcomplex)
static inline void spray_ans_with_Rbytes
	FUNDEF_spray_ans_with_vals(Rbyte)

static inline void spray_ans_with_CHARSXPs(
	const int *offs, SEXP vals,
	unsigned long long int outer_inc0,
	unsigned long long int outer_offset0,
	int *nzcount_buf, int **quick_out_offs_p, void *quick_out_vals_p,
	int inner_idx)
{
	int n, k, k2;
	SEXP *out_vals_p;
	unsigned long long int outer_idx;

	n = LENGTH(vals);
	out_vals_p = quick_out_vals_p;
	for (k = 0; k < n; k++) {
		outer_idx = outer_offset0 + outer_inc0 * *offs;
		k2 = nzcount_buf[outer_idx]++;
		quick_out_offs_p[outer_idx][k2] = inner_idx;
		SET_STRING_ELT(out_vals_p[outer_idx], k2, STRING_ELT(vals, k));
		offs++;
	}
	return;
}

static inline void spray_ans_with_list_elts(
	const int *offs, SEXP vals,
	unsigned long long int outer_inc0,
	unsigned long long int outer_offset0,
	int *nzcount_buf, int **quick_out_offs_p, void *quick_out_vals_p,
	int inner_idx)
{
	int n, k, k2;
	SEXP *out_vals_p;
	unsigned long long int outer_idx;

	n = LENGTH(vals);
	out_vals_p = quick_out_vals_p;
	for (k = 0; k < n; k++) {
		outer_idx = outer_offset0 + outer_inc0 * *offs;
		k2 = nzcount_buf[outer_idx]++;
		quick_out_offs_p[outer_idx][k2] = inner_idx;
		SET_VECTOR_ELT(out_vals_p[outer_idx], k2, VECTOR_ELT(vals, k));
		offs++;
	}
	return;
}

static void spray_ans_with_lv(SEXP lv, SEXPTYPE Rtype,
	 unsigned long long int outer_inc0,
	 unsigned long long int outer_offset0,
	 int *nzcount_buf, int **quick_out_offs_p, void *quick_out_vals_p,
	 int inner_idx)
{
	int lv_len;
	SEXP lv_offs, lv_vals;
	void (*spray_FUN)(const int *offs, SEXP vals,
			  unsigned long long int outer_inc0,
			  unsigned long long int outer_offset0,
			  int *nzcount_buf,
			  int **quick_out_offs_p, void *quick_out_vals_p,
			  int inner_idx);

	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	if (lv_len < 0)	
		error("SparseArray internal error in spray_ans_with_lv():\n"
		      "    invalid leaf vector");
	switch (Rtype) {
	    case LGLSXP: case INTSXP: spray_FUN = spray_ans_with_ints; break;
	    case REALSXP: spray_FUN = spray_ans_with_doubles; break;
	    case CPLXSXP: spray_FUN = spray_ans_with_Rcomplexes; break;
	    case RAWSXP: spray_FUN = spray_ans_with_Rbytes; break;
	    case STRSXP: spray_FUN = spray_ans_with_CHARSXPs; break;
	    case VECSXP: spray_FUN = spray_ans_with_list_elts; break;
	    default:
		error("SparseArray internal error in spray_ans_with_lv():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));
	}
	spray_FUN(INTEGER(lv_offs), lv_vals,
		  outer_inc0, outer_offset0,
		  nzcount_buf, quick_out_offs_p, quick_out_vals_p,
		  inner_idx);
	return;
}

/* Recursive. */
static void REC_spray_ans_with_SVT(SEXP SVT, int ndim, SEXPTYPE Rtype,
		const unsigned long long int *outer_incs,
		unsigned long long int outer_offset0,
		int *nzcount_buf,
		int **quick_out_offs_p, void *quick_out_vals_p)
{
	unsigned long long int outer_inc;
	int SVT_len, i;
	SEXP subSVT;
	static int inner_idx;

	if (SVT == R_NilValue)
		return;
	outer_inc = outer_incs[ndim - 1];
	if (ndim == 1) {
		spray_ans_with_lv(SVT, Rtype,
			outer_inc, outer_offset0,
			nzcount_buf, quick_out_offs_p, quick_out_vals_p,
			inner_idx);
		return;
	}
	SVT_len = LENGTH(SVT);
	for (i = 0; i < SVT_len; i++, outer_offset0 += outer_inc) {
		/* 'outer_inc == 0' means we're looping on the **rightmost**
		   dimension (a.k.a. outermost dimension) of 'SVT'. */
		if (outer_inc == 0)
			inner_idx = i;
		subSVT = VECTOR_ELT(SVT, i);
		REC_spray_ans_with_SVT(subSVT, ndim - 1, Rtype,
				outer_incs, outer_offset0,
				nzcount_buf,
				quick_out_offs_p, quick_out_vals_p);
	}
	return;
}

static SEXP transpose_SVT(SEXP SVT, const int *dim, int ndim, SEXPTYPE Rtype,
		const int *ans_dim, TransposeBufs *TBufs)
{
	SEXP ans;

	if (ndim <= 1 || SVT == R_NilValue)
		return SVT;

	/* 1st pass */
	memset(TBufs->nzcount_buf, 0, sizeof(int) * TBufs->nzcount_buf_len);
	REC_count_SVT_nzvals(SVT, ndim, TBufs->outer_incs, 0,
			     TBufs->nzcount_buf);
	/* 2nd pass */
	ans = PROTECT(REC_grow_and_alloc_leaves(ans_dim, ndim, Rtype,
						TBufs->nzcount_buf_incs,
						TBufs->nzcount_buf,
						TBufs->quick_out_offs_p,
						TBufs->quick_out_vals_p));
	/* 3rd pass */
	memset(TBufs->nzcount_buf, 0, sizeof(int) * TBufs->nzcount_buf_len);
	REC_spray_ans_with_SVT(SVT, ndim, Rtype,
			       TBufs->outer_incs, 0,
			       TBufs->nzcount_buf,
			       TBufs->quick_out_offs_p,
			       TBufs->quick_out_vals_p);
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_transpose_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT)
{
	SEXPTYPE Rtype;
	int ndim, *ans_dim, along;
	const int *dim;
	TransposeBufs TBufs;

	Rtype = _get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("SparseArray internal error in "
		      "C_transpose_SVT():\n"
		      "    SVT_SparseArray object has invalid type");

	ndim = LENGTH(x_dim);
	dim = INTEGER(x_dim);
	ans_dim = (int *) R_alloc(ndim, sizeof(int));
	for (along = 0; along < ndim; along++)
		ans_dim[along] = dim[ndim - 1 - along];
	TBufs = new_TBufs(dim, ndim, Rtype);
	return transpose_SVT(x_SVT, dim, ndim, Rtype, ans_dim, &TBufs);
}


/****************************************************************************
 * Permutation of the dimensions (most general case)
 */

/* Returns the "perm rank", that is, the number of leftmost dimensions that
   actually get permuted. In other words, the "perm rank" is the smallest
   number 'r' for which 'tail(perm, n=-r) ' is identical to '(r+1):ndim'.
   A "perm rank" of 0 means that 'perm' is '1:ndim', which is the identity
   permutation (no-op). */
static int compute_aperm_ans_dim(const int *dim, int ndim,
				 const int *perm, int *ans_dim)
{
	int *p_is_taken, along, p;

	p_is_taken = (int *) R_alloc(ndim, sizeof(int));
	memset(p_is_taken, 0, sizeof(int) * ndim);
	for (along = 0; along < ndim; along++) {
		p = perm[along];
		if (p == NA_INTEGER || p < 1 || p > ndim)
			error("invalid 'perm' argument");
		p--;
		if (p_is_taken[p])
			error("'perm' cannot contain duplicates");
		ans_dim[along] = dim[p];
		p_is_taken[p] = 1;
	}
	/* Compute the "perm rank". */
	for (along = ndim - 1; along >= 0; along--) {
		p = perm[along];
		if (p != along + 1)
			break;
	}
	//printf("perm_rank = %d\n", along + 1);

	/* TEMPORARY RESTRICTION: We support only dim permutations that
	   transpose the **leftmost dimensions** (a.k.a. innermost or
	   fastest moving dimensions in the memory layout). In other words,
	   we support only a 'perm' that is identical to 'c(r:1, (r+1):ndim)',
	   where 'r' is the "perm rank" ('r' is 'along' + 1).
	   TODO: Support any permutation. */
	for (int along2 = 0; along2 <= along; along2++)
		if (perm[along2] - 1 + along2 != along)
			error("we only support transposing the leftmost "
			      "dimensions at the moment, sorry");

	return along + 1;  /* "perm rank" */
}

/* Recursive. 'ndim' is guaranteed to be >= "perm rank". */
static SEXP REC_aperm_SVT(SEXP SVT, const int *dim, int ndim, SEXPTYPE Rtype,
		const int *perm, const int *ans_dim, TransposeBufs *TBufs)
{
	int ans_len, i;
	SEXP ans, subSVT, ans_elt;

	if (perm[ndim - 1] != ndim) {
		/* 'ndim' is equal to the "perm rank".
		   We know 'head(perm, n=ndim)' is guaranteed to represent a
		   transposition, i.e., that it's identical to 'ndim:1'.
		   See TEMPORARY RESTRICTION above. */
		return transpose_SVT(SVT, dim, ndim, Rtype, ans_dim, TBufs);
	}

	/* 'ndim' > "perm rank". */
	if (SVT == R_NilValue)
		return SVT;

	ans_len = LENGTH(SVT);
	ans = PROTECT(NEW_LIST(ans_len));
	for (i = 0; i < ans_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		ans_elt = PROTECT(
			REC_aperm_SVT(subSVT, dim, ndim - 1, Rtype, perm,
				      ans_dim, TBufs)
		);
		SET_VECTOR_ELT(ans, i, ans_elt);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_aperm_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP perm)
{
	SEXPTYPE Rtype;
	int ndim, *ans_dim, perm_rank;
	const int *dim;
	TransposeBufs TBufs;

	Rtype = _get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("SparseArray internal error in "
		      "C_aperm_SVT():\n"
		      "    SVT_SparseArray object has invalid type");

	ndim = LENGTH(x_dim);
	dim = INTEGER(x_dim);
	if (!IS_INTEGER(perm))
		error("'perm' must be an integer vector");
	if (LENGTH(perm) != ndim)
		error("'length(perm)' not equal to number "
		      "of dimensions of object to transpose");
	ans_dim = (int *) R_alloc(ndim, sizeof(int));
	perm_rank = compute_aperm_ans_dim(dim, ndim,
					  INTEGER(perm), ans_dim);
	if (perm_rank == 0)  /* identity permutation (no-op) */
		return x_SVT;

	TBufs = new_TBufs(dim, perm_rank, Rtype);
	return REC_aperm_SVT(x_SVT, dim, ndim, Rtype, INTEGER(perm),
			     ans_dim, &TBufs);
}

