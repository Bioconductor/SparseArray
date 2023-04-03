/****************************************************************************
 *               Core manipulation of SVT_SparseArray objects               *
 ****************************************************************************/
#include "SVT_SparseArray_class.h"

#include "S4Vectors_interface.h"  /* for sort_ints() */

#include "Rvector_utils.h"
#include "coerceVector2.h"  /* for _CoercionWarning() */
#include "leaf_vector_utils.h"

#include <limits.h>  /* for INT_MAX */
#include <string.h>  /* for memset() */


/* General purpose copy function. */
static inline int copy_Rvector_elts(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset,
		R_xlen_t nelt)
{
	SEXPTYPE Rtype;
	CopyRVectorElts_FUNType copy_Rvector_elts_FUN;

	Rtype = TYPEOF(in);
	copy_Rvector_elts_FUN = _select_copy_Rvector_elts_FUN(Rtype);
	if (copy_Rvector_elts_FUN == NULL)
		return -1;
	if (TYPEOF(out) != Rtype)
		return -1;
	if (in_offset  + nelt > XLENGTH(in))
		return -1;
	if (out_offset + nelt > XLENGTH(out))
		return -1;
	copy_Rvector_elts_FUN(in, in_offset, out, out_offset, nelt);
	return 0;
}


/****************************************************************************
 * C_get_SVT_SparseArray_nzcount()
 */

/* Recursive. */
static R_xlen_t REC_get_SVT_nzcount(SEXP SVT, int ndim)
{
	R_xlen_t nzcount;
	int SVT_len, i;
	SEXP subSVT;

	if (SVT == R_NilValue)
		return 0;

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		return XLENGTH(VECTOR_ELT(SVT, 0));
	}

	/* 'SVT' is a regular node (list). */
	nzcount = 0;
	SVT_len = LENGTH(SVT);
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		nzcount += REC_get_SVT_nzcount(subSVT, ndim - 1);
	}
	return nzcount;
}

/* --- .Call ENTRY POINT --- */
SEXP C_get_SVT_SparseArray_nzcount(SEXP x_dim, SEXP x_SVT)
{
	R_xlen_t nzcount;

	nzcount = REC_get_SVT_nzcount(x_SVT, LENGTH(x_dim));
	if (nzcount > INT_MAX)
		return ScalarReal((double) nzcount);
	return ScalarInteger((int) nzcount);
}


/****************************************************************************
 * type() setter
 */

/* Recursive. */
static int REC_set_SVT_type(SEXP SVT, const int *dim, int ndim,
		SEXPTYPE new_Rtype, int *warn, int *offs_buf)
{
	SEXP new_lv, subSVT;
	int SVT_len, is_empty, i, ret;

	if (SVT == R_NilValue)
		return 1;

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		new_lv = _coerce_leaf_vector(SVT, new_Rtype, warn, offs_buf);
		if (new_lv == R_NilValue)
			return 1;
		PROTECT(new_lv);
		SET_VECTOR_ELT(SVT, 0, VECTOR_ELT(new_lv, 0));
		SET_VECTOR_ELT(SVT, 1, VECTOR_ELT(new_lv, 1));
		UNPROTECT(1);
		return 0;
	}

	/* 'SVT' is a regular node (list). */
	SVT_len = LENGTH(SVT);

	/* Sanity check (should never fail). */
	if (SVT_len != dim[ndim - 1])
		return -1;

	is_empty = 1;
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		ret = REC_set_SVT_type(subSVT, dim, ndim - 1,
				       new_Rtype, warn, offs_buf);
		if (ret < 0)
			return -1;
		if (ret == 1) {
			SET_VECTOR_ELT(SVT, i, R_NilValue);
		} else {
			is_empty = 0;
		}
	}
	return is_empty;
}

/* --- .Call ENTRY POINT --- */
SEXP C_set_SVT_SparseArray_type(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		SEXP new_type)
{
	SEXPTYPE new_Rtype;
	int warn, *offs_buf, ret;
	SEXP ans;

	new_Rtype = _get_Rtype_from_Rstring(new_type);
	if (new_Rtype == 0)
		error("invalid supplied type");

	if (x_SVT == R_NilValue)
		return x_SVT;

	warn = 0;
	offs_buf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	ans = PROTECT(duplicate(x_SVT));
	ret = REC_set_SVT_type(ans, INTEGER(x_dim), LENGTH(x_dim),
			       new_Rtype, &warn, offs_buf);
	if (ret < 0) {
		UNPROTECT(1);
		error("SparseArray internal error in "
		      "C_set_SVT_SparseArray_type():\n"
		      "    REC_set_SVT_type() returned an error");
	}
	if (warn)
		_CoercionWarning(warn);
	UNPROTECT(1);
	return ret == 1 ? R_NilValue : ans;
}


/****************************************************************************
 * Going from SVT_SparseArray to ordinary R array
 */

/* Recursive. */
static int REC_dump_SVT_to_Rsubarray(SEXP SVT,
		const int *dim, int ndim,
		SEXP Rarray, R_xlen_t subarr_offset, R_xlen_t subarr_len)
{
	int SVT_len, i, ret;
	SEXP subSVT;

	if (SVT == R_NilValue)
		return 0;

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		return _expand_leaf_vector(SVT, Rarray, subarr_offset);
	}

	/* 'SVT' is a regular node (list). */
	SVT_len = LENGTH(SVT);
	if (SVT_len != dim[ndim - 1])
		return -1;

	subarr_len /= SVT_len;
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		ret = REC_dump_SVT_to_Rsubarray(subSVT,
				dim, ndim - 1,
				Rarray, subarr_offset, subarr_len);
		if (ret < 0)
			return -1;
		subarr_offset += subarr_len;
	}
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_from_SVT_SparseArray_to_Rarray(SEXP x_dim, SEXP x_dimnames,
		SEXP x_type, SEXP x_SVT)
{
	SEXPTYPE Rtype;
	SEXP ans;
	int ret;

	Rtype = _get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("SparseArray internal error in "
		      "C_from_SVT_SparseArray_to_Rarray():\n"
		      "    SVT_SparseArray object has invalid type");

	ans = PROTECT(_new_Rarray(Rtype, x_dim, x_dimnames));
	ret = REC_dump_SVT_to_Rsubarray(x_SVT,
				INTEGER(x_dim), LENGTH(x_dim),
				ans, 0, XLENGTH(ans));
	UNPROTECT(1);
	if (ret < 0)
		error("SparseArray internal error "
		      "in C_from_SVT_SparseArray_to_Rarray():\n"
		      "    invalid SVT_SparseArray object");
	return ans;
}


/****************************************************************************
 * Going from ordinary R array to SVT_SparseArray
 */

/* Recursive. */
static SEXP REC_build_SVT_from_Rsubarray(
		SEXP Rarray, R_xlen_t subarr_offset, R_xlen_t subarr_len,
		const int *dim, int ndim,
		SEXPTYPE ans_Rtype, int *warn, int *offs_buf)
{
	SEXP ans, ans_elt;
	int SVT_len, is_empty, i;

	if (ndim == 1) {
		/* Sanity check (should never fail). */
		if (dim[0] != subarr_len)
			error("SparseArray internal error "
			      "in REC_build_SVT_from_Rsubarray():\n"
			      "    dim[0] != subarr_len");
		ans = _make_leaf_vector_from_Rsubvec(
					Rarray, subarr_offset, dim[0],
					offs_buf, 1);
		if (ans_Rtype == TYPEOF(Rarray) || ans == R_NilValue)
			return ans;
		PROTECT(ans);
		ans = _coerce_leaf_vector(ans, ans_Rtype, warn, offs_buf);
		UNPROTECT(1);
		return ans;
	}

	SVT_len = dim[ndim - 1];  /* cannot be 0 so safe to divide below */
	subarr_len /= SVT_len;
	ans = PROTECT(NEW_LIST(SVT_len));
	is_empty = 1;
	for (i = 0; i < SVT_len; i++) {
		ans_elt = REC_build_SVT_from_Rsubarray(
					Rarray, subarr_offset, subarr_len,
					dim, ndim - 1,
					ans_Rtype, warn, offs_buf);
		if (ans_elt != R_NilValue) {
			PROTECT(ans_elt);
			SET_VECTOR_ELT(ans, i, ans_elt);
			UNPROTECT(1);
			is_empty = 0;
		}
		subarr_offset += subarr_len;
	}
	UNPROTECT(1);
	return is_empty ? R_NilValue : ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_build_SVT_from_Rarray(SEXP x, SEXP ans_type)
{
	SEXPTYPE ans_Rtype;
	int x_ndim, warn, *offs_buf;
	R_xlen_t x_len;
	SEXP x_dim, ans;

	ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (ans_Rtype == 0)
		error("invalid requested type");

	x_len = XLENGTH(x);
	if (x_len == 0)  /* means that 'any(dim(x) == 0)' is TRUE */
		return R_NilValue;

	x_dim = GET_DIM(x);  /* does not contain zeros */
	x_ndim = LENGTH(x_dim);
	warn = 0;
	offs_buf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	ans = REC_build_SVT_from_Rsubarray(x, 0, x_len,
					   INTEGER(x_dim), x_ndim,
					   ans_Rtype, &warn, offs_buf);
	if (warn) {
		if (ans != R_NilValue)
			PROTECT(ans);
		_CoercionWarning(warn);
		if (ans != R_NilValue)
			UNPROTECT(1);
	}
	return ans;
}


/****************************************************************************
 * Going from SVT_SparseMatrix to [d|l]gCMatrix
 */

/* Return nb of nonzero values in column. */
static int dump_col_to_CsparseMatrix_slots(SEXP SVT, int col_idx,
		SEXP ans_i, SEXP ans_x, int offset)
{
	SEXP subSVT, lv_offs, lv_vals;
	int lv_len, ret;

	subSVT = VECTOR_ELT(SVT, col_idx);
	if (subSVT == R_NilValue)
		return 0;

	/* 'subSVT' is a "leaf vector". */
	lv_len = _split_leaf_vector(subSVT, &lv_offs, &lv_vals);
	if (lv_len < 0)
		return -1;

	/* Copy 0-based row indices from 'lv_offs' to 'ans_i'. */
	_copy_INTEGER_elts(lv_offs, (R_xlen_t) 0,
			ans_i, (R_xlen_t) offset,
			XLENGTH(lv_offs));

	ret = copy_Rvector_elts(lv_vals, (R_xlen_t) 0,
			ans_x, (R_xlen_t) offset,
			XLENGTH(lv_vals));
	if (ret < 0)
		return -1;

	return lv_len;
}

static int dump_SVT_to_CsparseMatrix_slots(SEXP x_SVT, int x_ncol,
		SEXP ans_p, SEXP ans_i, SEXP ans_x)
{
	int offset, j, nzcount;

	INTEGER(ans_p)[0] = 0;
	offset = 0;
	for (j = 0; j < x_ncol; j++) {
		nzcount = dump_col_to_CsparseMatrix_slots(x_SVT, j,
						ans_i, ans_x, offset);
		if (nzcount < 0)
			return -1;
		offset += nzcount;
		INTEGER(ans_p)[j + 1] = offset;
	}
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_from_SVT_SparseMatrix_to_CsparseMatrix(SEXP x_dim,
		SEXP x_type, SEXP x_SVT)
{
	R_xlen_t nzcount;
	SEXPTYPE x_Rtype;
	int x_ncol, ret;
	SEXP ans_p, ans_i, ans_x, ans;

	if (LENGTH(x_dim) != 2)
		error("object to coerce to [d|l]gCMatrix "
		      "must have exactly 2 dimensions");

	nzcount = REC_get_SVT_nzcount(x_SVT, 2);
	if (nzcount > INT_MAX)
		error("SVT_SparseMatrix object contains too many nonzero "
		      "values to be turned into a dgCMatrix or lgCMatrix "
		      "object");

	x_Rtype = _get_Rtype_from_Rstring(x_type);
	if (x_Rtype == 0)
		error("SparseArray internal error in "
		      "C_from_SVT_SparseMatrix_to_CsparseMatrix():\n"
		      "    SVT_SparseMatrix object has invalid type");

	x_ncol = INTEGER(x_dim)[1];

	ans_i = PROTECT(NEW_INTEGER(nzcount));
	ans_x = PROTECT(allocVector(x_Rtype, nzcount));
	if (nzcount == 0) {
		ans_p = PROTECT(_new_Rvector(INTSXP, (R_xlen_t) x_ncol + 1));
	} else {
		ans_p = PROTECT(NEW_INTEGER(x_ncol + 1));
		ret = dump_SVT_to_CsparseMatrix_slots(x_SVT, x_ncol,
						      ans_p, ans_i, ans_x);
		if (ret < 0) {
			UNPROTECT(3);
			error("SparseArray internal error in "
			      "C_from_SVT_SparseMatrix_to_CsparseMatrix():\n"
			      "    invalid SVT_SparseMatrix object");
		}
	}

	ans = PROTECT(NEW_LIST(3));
	SET_VECTOR_ELT(ans, 0, ans_p);
	SET_VECTOR_ELT(ans, 1, ans_i);
	SET_VECTOR_ELT(ans, 2, ans_x);
	UNPROTECT(4);
	return ans;
}


/****************************************************************************
 * Going from [d|l]gCMatrix to SVT_SparseMatrix
 */

/* Returns R_NilValue or a "leaf vector" of length <= nzcount.

   Note that, quite surprisingly, and unfortunately, [d|l]gCMatrix objects
   can sometimes have zeros in their "x" slot.
   For example, with Matrix 1.3-4:

       dgcm <- as(matrix(c(11:13, 0, 0, 21:25), ncol=2), "dgCMatrix")
       dgcm[cbind(3:4, 2)] <- 0
       dgcm
       # 5 x 2 sparse Matrix of class "dgCMatrix"
       # [1,] 11 21
       # [2,] 12 22
       # [3,] 13  0
       # [4,]  .  0
       # [5,]  . 25
       dgcm@x
       # [1] 11 12 13 21 22  0  0 25

   This is still considered a valid dgCMatrix object:

       validObject(dgcm)
       [1] TRUE

   Interestingly this doesn't happen when using a linear index in the
   subassignment:

       dgcm[1:2] <- 0
       dgcm@x
       # [1] 13 21 22 25

   We want to make sure that these zeros don't end up in the "leaf vector"
   returned by build_leaf_vector_from_CsparseMatrix_col(). Unfortunately,
   this introduces an additional cost to coercion from [d|l]gCMatrix to
   SVT_SparseMatrix. This cost is a slowdown that is (approx.) between 1.3x
   and 1.5x. */
static SEXP build_leaf_vector_from_CsparseMatrix_col(SEXP x_i, SEXP x_x,
		int offset, int nzcount,
		SEXPTYPE ans_Rtype, int *warn, int *offs_buf)
{
	SEXP ans, ans_offs;
	int ans_len;

	/* Will skip zeros from 'x_x' if any. */
	ans = _make_leaf_vector_from_Rsubvec(x_x, (R_xlen_t) offset, nzcount,
					     offs_buf, 1);
	if (ans == R_NilValue)
		return ans;
	PROTECT(ans);
	ans_offs = VECTOR_ELT(ans, 0);
	ans_len = LENGTH(ans_offs);  /* can only be <= 'nzcount' */
	_copy_selected_ints(INTEGER(x_i) + offset, INTEGER(ans_offs), ans_len,
			    INTEGER(ans_offs));
	if (ans_Rtype != TYPEOF(x_x))
		ans = _coerce_leaf_vector(ans, ans_Rtype, warn, offs_buf);
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_build_SVT_from_CsparseMatrix(SEXP x, SEXP ans_type)
{
	SEXPTYPE ans_Rtype;
	SEXP x_Dim, x_p, x_i, x_x, ans, ans_elt;
	int x_nrow, x_ncol, j, offset, nzcount, warn, *offs_buf, is_empty;

	ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (ans_Rtype == 0)
		error("invalid requested type");

	x_Dim = GET_SLOT(x, install("Dim"));
	x_nrow = INTEGER(x_Dim)[0];
	x_ncol = INTEGER(x_Dim)[1];
	x_p = GET_SLOT(x, install("p"));

	if (INTEGER(x_p)[x_ncol] == 0)
		return R_NilValue;

	x_i = GET_SLOT(x, install("i"));
	x_x = GET_SLOT(x, install("x"));

	warn = 0;
	offs_buf = (int *) R_alloc(x_nrow, sizeof(int));
	ans = PROTECT(NEW_LIST(x_ncol));
	is_empty = 1;
	for (j = 0; j < x_ncol; j++) {
		offset = INTEGER(x_p)[j];
		nzcount = INTEGER(x_p)[j + 1] - offset;
		if (nzcount != 0) {
			ans_elt = build_leaf_vector_from_CsparseMatrix_col(
						x_i, x_x,
						offset, nzcount,
						ans_Rtype, &warn, offs_buf);
			if (ans_elt != R_NilValue) {
				PROTECT(ans_elt);
				SET_VECTOR_ELT(ans, j, ans_elt);
				UNPROTECT(1);
				is_empty = 0;
			}
		}
	}
	if (warn)
		_CoercionWarning(warn);
	UNPROTECT(1);
	return is_empty ? R_NilValue : ans;
}


/****************************************************************************
 * Going from SVT_SparseArray to COO_SparseArray
 */

static SEXP alloc_nzvals(SEXP type, R_xlen_t n)
{
	SEXPTYPE Rtype;

	Rtype = _get_Rtype_from_Rstring(type);
	if (Rtype == 0)
		error("SparseArray internal error in alloc_nzvals():\n"
		      "    SVT_SparseArray object has invalid type");
	return allocVector(Rtype, n);
}

/* Recursive. */
static int REC_extract_nzcoo_and_nzvals_from_SVT(SEXP SVT,
		SEXP nzvals, int *nzdata_offset,
		int *nzcoo, int nzcoo_nrow, int nzcoo_ncol,
		int *rowbuf, int rowbuf_offset)
{
	int SVT_len, i, ret, lv_len, k, *p, j;
	SEXP subSVT, lv_offs, lv_vals;

	if (SVT == R_NilValue)
		return 0;

	if (rowbuf_offset > 0) {
		if (!isVectorList(SVT))  // IS_LIST() is broken
			return -1;
		SVT_len = LENGTH(SVT);
		for (i = 0; i < SVT_len; i++) {
			subSVT = VECTOR_ELT(SVT, i);
			rowbuf[rowbuf_offset] = i + 1;
			ret = REC_extract_nzcoo_and_nzvals_from_SVT(
					subSVT,
					nzvals, nzdata_offset,
					nzcoo, nzcoo_nrow, nzcoo_ncol,
					rowbuf, rowbuf_offset - 1);
			if (ret < 0)
				return -1;
		}
		return 0;
	}

	/* 'SVT' is a "leaf vector". */
	lv_len = _split_leaf_vector(SVT, &lv_offs, &lv_vals);
	if (lv_len < 0)
		return -1;

	ret = copy_Rvector_elts(lv_vals, (R_xlen_t) 0,
				nzvals, (R_xlen_t) *nzdata_offset,
				XLENGTH(lv_vals));
	if (ret < 0)
		return -1;

	for (k = 0; k < lv_len; k++) {
		rowbuf[0] = INTEGER(lv_offs)[k] + 1;

		/* Copy 'rowbuf' to 'nzcoo'. */
		p = nzcoo + *nzdata_offset;
		for (j = 0; j < nzcoo_ncol; j++) {
			*p = rowbuf[j];
			p += nzcoo_nrow;
		}

		(*nzdata_offset)++;
	}
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_from_SVT_SparseArray_to_COO_SparseArray(SEXP x_dim,
		SEXP x_type, SEXP x_SVT)
{
	R_xlen_t nzcount;
	int nzcoo_nrow, nzcoo_ncol, *rowbuf, nzdata_offset, ret;
	SEXP nzcoo, nzvals, ans;

	nzcount = REC_get_SVT_nzcount(x_SVT, LENGTH(x_dim));
	if (nzcount > INT_MAX)
		error("SVT_SparseArray object contains too many nonzero "
		      "values to be turned into a COO_SparseArray object");

	nzvals = PROTECT(alloc_nzvals(x_type, nzcount));

	nzcoo_nrow = (int) nzcount;
	nzcoo_ncol = LENGTH(x_dim);
	rowbuf = (int *) R_alloc(nzcoo_ncol, sizeof(int));
	nzcoo = PROTECT(allocMatrix(INTSXP, nzcoo_nrow, nzcoo_ncol));

	nzdata_offset = 0;
	ret = REC_extract_nzcoo_and_nzvals_from_SVT(x_SVT,
			nzvals, &nzdata_offset,
			INTEGER(nzcoo), nzcoo_nrow, nzcoo_ncol,
			rowbuf, nzcoo_ncol - 1);
	if (ret < 0) {
		UNPROTECT(2);
		error("SparseArray internal error "
		      "in C_from_SVT_SparseArray_to_COO_SparseArray():\n"
		      "    invalid SVT_SparseArray object");
	}

	/* Sanity check (should never fail). */
	if (nzdata_offset != nzcoo_nrow) {
		UNPROTECT(2);
		error("SparseArray internal error "
		      "in C_from_SVT_SparseArray_to_COO_SparseArray():\n"
		      "    *out_offset != nzcoo_nrow");
	}

	ans = PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(ans, 0, nzcoo);
	SET_VECTOR_ELT(ans, 1, nzvals);
	UNPROTECT(3);
	return ans;
}


/****************************************************************************
 * Transposition
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
		SEXP out_SVT, int *nzcounts);

/* Ignores 'out_SVT' and 'nzcounts'. */
static void transpose_INTEGER_col(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcounts)
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

/* Ignores 'out_SVT' and 'nzcounts'. */
static void transpose_NUMERIC_col(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcounts)
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

/* Ignores 'out_SVT' and 'nzcounts'. */
static void transpose_COMPLEX_col(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcounts)
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

/* Ignores 'out_SVT' and 'nzcounts'. */
static void transpose_RAW_col(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcounts)
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
		SEXP out_SVT, int *nzcounts)
{
	int lv_len, k, row_idx;
	SEXP out_lv;

	lv_len = LENGTH(lv_vals);
	for (k = 0; k < lv_len; k++) {
		row_idx = *offs;
		*(quick_out_offs_p[row_idx]++) = col_idx;
		out_lv = VECTOR_ELT(out_SVT, row_idx);
		_copy_CHARACTER_elt(lv_vals, (R_xlen_t) k,
			VECTOR_ELT(out_lv, 1), (R_xlen_t) nzcounts[row_idx]++);
		offs++;
	}
	return;
}

/* Ignores 'quick_out_vals_p'. */
static void transpose_LIST_col(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcounts)
{
	int lv_len, k, row_idx;
	SEXP out_lv;

	lv_len = LENGTH(lv_vals);
	for (k = 0; k < lv_len; k++) {
		row_idx = *offs;
		*(quick_out_offs_p[row_idx]++) = col_idx;
		out_lv = VECTOR_ELT(out_SVT, row_idx);
		_copy_LIST_elt(lv_vals, (R_xlen_t) k,
			VECTOR_ELT(out_lv, 1), (R_xlen_t) nzcounts[row_idx]++);
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

static SEXP transpose_SVT(SEXP SVT, SEXPTYPE Rtype, int nrow, int ncol,
		int *nzcounts)
{
	TransposeCol_FUNType transpose_col_FUN;
	SEXP ans, ans_elt, subSVT, lv_offs, lv_vals;
	int **quick_out_offs_p;
	void **quick_out_vals_p;
	int i, j, lv_len;

	transpose_col_FUN = select_transpose_col_FUN(Rtype);
	if (transpose_col_FUN == NULL)
		error("SparseArray internal error in "
		      "transpose_SVT():\n"
		      "    SVT_SparseMatrix object has invalid type");

	ans = PROTECT(NEW_LIST(nrow));
	quick_out_offs_p = (int **) R_alloc(nrow, sizeof(int *));
	for (i = 0; i < nrow; i++) {
		lv_len = nzcounts[i];
		if (lv_len != 0) {
			ans_elt = PROTECT(_alloc_leaf_vector(lv_len, Rtype));
			SET_VECTOR_ELT(ans, i, ans_elt);
			UNPROTECT(1);
			quick_out_offs_p[i] = INTEGER(VECTOR_ELT(ans_elt, 0));
		}
	}
	quick_out_vals_p = set_quick_out_vals_p(ans, Rtype);

	memset(nzcounts, 0, sizeof(int) * nrow);
	for (j = 0; j < ncol; j++) {
		subSVT = VECTOR_ELT(SVT, j);
		if (subSVT == R_NilValue)
			continue;
		/* 'subSVT' is a "leaf vector". */
		lv_len = _split_leaf_vector(subSVT, &lv_offs, &lv_vals);
		if (lv_len < 0) {
			UNPROTECT(1);
			error("SparseArray internal error in "
			      "transpose_SVT():\n"
			      "    invalid SVT_SparseMatrix object");
		}
		transpose_col_FUN(j,
			INTEGER(lv_offs), lv_vals,
			quick_out_offs_p, quick_out_vals_p,
			ans, nzcounts);
	}
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_transpose_SVT_SparseMatrix(SEXP x_dim, SEXP x_type, SEXP x_SVT)
{
	SEXPTYPE Rtype;
	int x_nrow, x_ncol, *nzcounts;

	Rtype = _get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("SparseArray internal error in "
		      "C_transpose_SVT_SparseMatrix():\n"
		      "    SVT_SparseMatrix object has invalid type");

	if (LENGTH(x_dim) != 2)
		error("object to transpose must have exactly 2 dimensions");

	if (x_SVT == R_NilValue)
		return x_SVT;

	x_nrow = INTEGER(x_dim)[0];
	x_ncol = INTEGER(x_dim)[1];
	nzcounts = (int *) R_alloc(x_nrow, sizeof(int));

	/* 1st pass: Count the number of nonzero values per row in the
	   input object. */
	count_nonzero_vals_per_row(x_SVT, x_nrow, x_ncol, nzcounts);

	/* 2nd pass: Build the transposed SVT. */
	return transpose_SVT(x_SVT, Rtype, x_nrow, x_ncol, nzcounts);
}

