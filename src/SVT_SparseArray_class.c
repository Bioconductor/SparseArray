/****************************************************************************
 *               Core manipulation of SVT_SparseArray objects               *
 ****************************************************************************/
#include "SVT_SparseArray_class.h"

#include "S4Vectors_interface.h"  /* for sort_ints() */

#include "argcheck_utils.h"
#include "Rvector_utils.h"
#include "coerceVector2.h"  /* for _CoercionWarning() */
#include "leaf_utils.h"

#include <limits.h>  /* for INT_MAX */
#include <string.h>  /* for strcmp() and memcpy() */


static void copy_nzvals_elts_to_Rsubvec(SEXP nzvals,
		SEXP out, int out_offset, int nelt)
{
	if (nzvals == R_NilValue) {
		/* lacunar leaf */
		_set_Rsubvec_elts_to_one(out, (R_xlen_t) out_offset,
					 (R_xlen_t) nelt);
	} else {
		/* regular leaf */
		_copy_Rvector_elts(nzvals, 0,
				   out, (R_xlen_t) out_offset,
				   (R_xlen_t) nelt);
	}
	return;
}


/****************************************************************************
 * type() setter
 */

/* Returns 1 if coercion to the requested type produces a NULL 'leaf',
   and 0 otherwise. */
static int INPLACE_modify_leaf_type(SEXP leaf, int na_background,
		SEXPTYPE new_Rtype, int *warn, int *nzoffs_buf)
{
	SEXP new_leaf;
	if (na_background) {
		new_leaf = _coerce_naleaf(leaf, new_Rtype, warn, nzoffs_buf);
	} else {
		new_leaf = _coerce_leaf(leaf, new_Rtype, warn, nzoffs_buf);
	}
	if (new_leaf == R_NilValue)
		return 1;
	PROTECT(new_leaf);
	replace_leaf_nzvals(leaf, get_leaf_nzvals(new_leaf));
	replace_leaf_nzoffs(leaf, get_leaf_nzoffs(new_leaf));
	UNPROTECT(1);
	return 0;
}

/* Recursive. */
static int REC_INPLACE_modify_SVT_type(SEXP SVT, int na_background,
		const int *dim, int ndim,
		SEXPTYPE new_Rtype, int *warn, int *nzoffs_buf)
{
	if (SVT == R_NilValue)
		return 1;

	if (ndim == 1) {
		/* 'SVT' is a leaf (i.e. a 1D SVT). */
		return INPLACE_modify_leaf_type(SVT, na_background,
						new_Rtype, warn, nzoffs_buf);
	}

	/* 'SVT' is a regular node (list). */
	int SVT_len = LENGTH(SVT);

	/* Sanity check (should never fail). */
	if (SVT_len != dim[ndim - 1])
		return -1;

	int is_empty = 1;
	for (int i = 0; i < SVT_len; i++) {
		SEXP subSVT = VECTOR_ELT(SVT, i);
		int ret = REC_INPLACE_modify_SVT_type(subSVT, na_background,
						      dim, ndim - 1,
						      new_Rtype,
						      warn, nzoffs_buf);
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
SEXP C_set_SVT_type(
		SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP x_na_background,
		SEXP new_type)
{
	SEXPTYPE x_Rtype = _get_and_check_Rtype_from_Rstring(x_type,
					"C_set_SVT_type", "x_type");
	SEXPTYPE new_Rtype = _get_Rtype_from_Rstring(new_type);
	if (new_Rtype == 0)
		error("invalid supplied type");

	if (new_Rtype == x_Rtype || x_SVT == R_NilValue)
		return x_SVT;

	int x_has_NAbg = _get_and_check_na_background(x_na_background,
					"C_set_SVT_type", "x_na_background");

	int *nzoffs_buf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	SEXP ans = PROTECT(duplicate(x_SVT));
	int warn = 0;
	int ret = REC_INPLACE_modify_SVT_type(ans, x_has_NAbg,
				INTEGER(x_dim), LENGTH(x_dim),
				new_Rtype, &warn, nzoffs_buf);
	if (ret < 0) {
		UNPROTECT(1);
		error("SparseArray internal error in "
		      "C_set_SVT_type():\n"
		      "    REC_INPLACE_modify_SVT_type() returned an error");
	}
	if (warn)
		_CoercionWarning(warn);
	UNPROTECT(1);
	return ret == 1 ? R_NilValue : ans;
}

/* Used in src/SparseArray_Arith_methods.c by REC_Arith_SVT1_SVT2().
   Assumes that 'to_Rtype' is equal or a bigger type than 'from_Rtype'
   so coercion won't truncate values or discard imaginary parts etc... */
SEXP _coerce_SVT(SEXP SVT, const int *dim, int ndim,
		 SEXPTYPE from_Rtype, SEXPTYPE to_Rtype, int *offs_buf)
{
        if (from_Rtype == to_Rtype)
		return SVT;
	SEXP ans = PROTECT(duplicate(SVT));
	int warn;
	int ret = REC_INPLACE_modify_SVT_type(ans, 0, dim, ndim,
					      to_Rtype, &warn, offs_buf);
	if (ret < 0) {
		UNPROTECT(1);
		error("SparseArray internal error in _coerce_SVT():\n"
		      "    REC_INPLACE_modify_SVT_type() returned an error");
	}
	UNPROTECT(1);
	return ret == 1 ? R_NilValue : ans;
}


/****************************************************************************
 * C_is_nonzero_SVT()
 */

static SEXP REC_is_nonzero_SVT(SEXP SVT, const int *dim, int ndim)
{
	if (SVT == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* 'SVT' is a leaf (i.e. 1D SVT). */
		return _make_lacunar_leaf(get_leaf_nzoffs(SVT));
	}

	/* 'SVT' is a list. */
	int ans_len = dim[ndim - 1];  /* same as 'LENGTH(SVT)' */
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	int is_empty = 1;
	for (int i = 0; i < ans_len; i++) {
		SEXP subSVT = VECTOR_ELT(SVT, i);
		SEXP ans_elt = REC_is_nonzero_SVT(subSVT, dim, ndim - 1);
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
   Simple strips off the "nzvals" part from the non-empty leaves.
   Returns a fully lacunar SVT. */
SEXP C_is_nonzero_SVT(SEXP x_dim, SEXP x_SVT)
{
	return REC_is_nonzero_SVT(x_SVT, INTEGER(x_dim), LENGTH(x_dim));
}


/****************************************************************************
 * C_nzcount_SVT()
 */

/* Recursive. */
R_xlen_t _REC_nzcount_SVT(SEXP SVT, int ndim)
{
	if (SVT == R_NilValue)
		return 0;

	if (ndim == 1) {
		/* 'SVT' is a leaf (i.e. a 1D SVT). */
		return (R_xlen_t) get_leaf_nzcount(SVT);
	}

	/* 'SVT' is a regular node (list). */
	R_xlen_t nzcount = 0;
	int SVT_len = LENGTH(SVT);
	for (int i = 0; i < SVT_len; i++) {
		SEXP subSVT = VECTOR_ELT(SVT, i);
		nzcount += _REC_nzcount_SVT(subSVT, ndim - 1);
	}
	return nzcount;
}

/* --- .Call ENTRY POINT --- */
SEXP C_nzcount_SVT(SEXP x_dim, SEXP x_SVT)
{
	R_xlen_t nzcount = _REC_nzcount_SVT(x_SVT, LENGTH(x_dim));
	if (nzcount > INT_MAX)
		return ScalarReal((double) nzcount);
	return ScalarInteger((int) nzcount);
}


/****************************************************************************
 * C_nzwhich_SVT()
 */

static inline void from_offs_to_int_Lindex(const int *offs, int n,
		int arr_offset, int *Lindex)
{
	for (int k = 0; k < n; k++)
		Lindex[k] = arr_offset + offs[k] + 1;  /* will never overflow */
	return;
}

static inline void from_offs_to_double_Lindex(const int *offs, int n,
		R_xlen_t arr_offset, double *Lindex)
{
	for (int k = 0; k < n; k++)
		Lindex[k] = (double) (arr_offset + offs[k] + 1);
	return;
}

/* Recursive. */
static int REC_nzwhich_SVT_as_Lindex(SEXP SVT,
		const int *dim, const R_xlen_t *dimcumprod, int ndim,
		R_xlen_t arr_offset, SEXP Lindex, R_xlen_t *Lindex_offset)
{
	if (SVT == R_NilValue)
		return 0;

	if (ndim == 1) {
		/* 'SVT' is a leaf (i.e. 1D SVT). */
		SEXP nzvals, nzoffs;
		int nzcount = unzip_leaf(SVT, &nzvals, &nzoffs);
		if (nzcount < 0)
			return -1;
		if (IS_INTEGER(Lindex)) {
			from_offs_to_int_Lindex(INTEGER(nzoffs), nzcount,
					(int) arr_offset,
					INTEGER(Lindex) + *Lindex_offset);
		} else {
			from_offs_to_double_Lindex(INTEGER(nzoffs), nzcount,
					arr_offset,
					REAL(Lindex) + *Lindex_offset);
		}
		(*Lindex_offset) += nzcount;
		return 0;
	}

	/* 'SVT' is a regular node (list). */
	int SVT_len = LENGTH(SVT);

	/* Sanity check (should never fail). */
	if (SVT_len != dim[ndim - 1])
		return -1;

	R_xlen_t subarr_len = dimcumprod[ndim - 2];
	for (int i = 0; i < SVT_len; i++) {
		SEXP subSVT = VECTOR_ELT(SVT, i);
		int ret = REC_nzwhich_SVT_as_Lindex(subSVT,
					dim, dimcumprod, ndim - 1,
					arr_offset, Lindex, Lindex_offset);
		if (ret < 0)
			return -1;
		arr_offset += subarr_len;
	}
	return 0;
}

static SEXP nzwhich_SVT_as_Lindex(SEXP SVT, const int *dim, int ndim,
		R_xlen_t nzcount)
{
	R_xlen_t *dimcumprod = (R_xlen_t *) R_alloc(ndim, sizeof(R_xlen_t));
	R_xlen_t p = 1;
	for (int along = 0; along < ndim; along++) {
		p *= dim[along];
		dimcumprod[along] = p;
	}

	SEXPTYPE ans_Rtype = p > INT_MAX ? REALSXP : INTSXP;
	SEXP ans = PROTECT(allocVector(ans_Rtype, nzcount));
	R_xlen_t Lindex_offset = 0;
	int ret = REC_nzwhich_SVT_as_Lindex(SVT, dim, dimcumprod, ndim,
					    0, ans, &Lindex_offset);
	UNPROTECT(1);
	if (ret < 0)
		error("SparseArray internal error in "
		      "nzwhich_SVT_as_Lindex():\n"
		      "    invalid SVT_SparseArray object");

	/* Sanity check (should never fail). */
	if (Lindex_offset != nzcount) {
		error("SparseArray internal error in "
		      "nzwhich_SVT_as_Lindex():\n"
		      "    Lindex_offset != nzcount");
	}
	return ans;
}

/* Recursive. */
static int REC_extract_nzcoo_and_nzvals_from_SVT(SEXP SVT,
		int *out_nzcoo, int nzcoo_nrow, int nzcoo_ncol,
		int *rowbuf, int rowbuf_offset,
		SEXP out_nzvals, int *nzvals_offset)
{
	if (SVT == R_NilValue)
		return 0;

	if (rowbuf_offset > 0) {
		if (!isVectorList(SVT))  // IS_LIST() is broken
			return -1;
		int SVT_len = LENGTH(SVT);
		for (int i = 0; i < SVT_len; i++) {
			SEXP subSVT = VECTOR_ELT(SVT, i);
			rowbuf[rowbuf_offset] = i + 1;
			int ret = REC_extract_nzcoo_and_nzvals_from_SVT(subSVT,
					out_nzcoo, nzcoo_nrow, nzcoo_ncol,
					rowbuf, rowbuf_offset - 1,
					out_nzvals, nzvals_offset);
			if (ret < 0)
				return -1;
		}
		return 0;
	}

	/* 'SVT' is a leaf (i.e. a 1D SVT). */
	SEXP leaf_nzvals, leaf_nzoffs;
	int leaf_nzcount = unzip_leaf(SVT, &leaf_nzvals, &leaf_nzoffs);

	if (out_nzvals != R_NilValue)
		copy_nzvals_elts_to_Rsubvec(leaf_nzvals,
					    out_nzvals, *nzvals_offset,
					    leaf_nzcount);

	for (int k = 0; k < leaf_nzcount; k++) {
		rowbuf[0] = INTEGER(leaf_nzoffs)[k] + 1;

		/* Copy 'rowbuf' to 'nzcoo'. */
		int *p = out_nzcoo + *nzvals_offset;
		for (int j = 0; j < nzcoo_ncol; j++) {
			*p = rowbuf[j];
			p += nzcoo_nrow;
		}

		(*nzvals_offset)++;
	}
	return 0;
}

static SEXP extract_nzcoo_and_nzvals_from_SVT(SEXP SVT,
		int nzcoo_nrow, int nzcoo_ncol,
		SEXP out_nzvals)
{
	int *rowbuf = (int *) R_alloc(nzcoo_ncol, sizeof(int));
	SEXP out_nzcoo = PROTECT(allocMatrix(INTSXP, nzcoo_nrow, nzcoo_ncol));

	int nzvals_offset = 0;
	int ret = REC_extract_nzcoo_and_nzvals_from_SVT(SVT,
			INTEGER(out_nzcoo), nzcoo_nrow, nzcoo_ncol,
			rowbuf, nzcoo_ncol - 1,
			out_nzvals, &nzvals_offset);
	if (ret < 0) {
		UNPROTECT(1);
		error("SparseArray internal error in "
		      "extract_nzcoo_and_nzvals_from_SVT():\n"
		      "    invalid SVT_SparseArray object");
	}

	/* Sanity check (should never fail). */
	if (nzvals_offset != nzcoo_nrow) {
		UNPROTECT(1);
		error("SparseArray internal error in "
		      "extract_nzcoo_and_nzvals_from_SVT():\n"
		      "    nzvals_offset != nzcoo_nrow");
	}

	UNPROTECT(1);
	return out_nzcoo;
}

/* --- .Call ENTRY POINT --- */
SEXP C_nzwhich_SVT(SEXP x_dim, SEXP x_SVT, SEXP arr_ind)
{
	int x_ndim = LENGTH(x_dim);
	R_xlen_t nzcount = _REC_nzcount_SVT(x_SVT, x_ndim);

	if (!LOGICAL(arr_ind)[0]) {
		/* Return linear indices of nonzero array elements in an
		   integer or numeric vector representing an L-index. */
		return nzwhich_SVT_as_Lindex(x_SVT, INTEGER(x_dim), x_ndim,
					     nzcount);
	}

	/* Return coordinates of nonzero array elements in an integer matrix
	   representing an M-index. */
	if (nzcount > INT_MAX)
		error("too many nonzero elements in SVT_SparseArray "
		      "object to return their \"array\n  coordinates\" "
		      "(n-tuples) in a matrix");
	return extract_nzcoo_and_nzvals_from_SVT(x_SVT, (int) nzcount, x_ndim,
						 R_NilValue);
}


/****************************************************************************
 * Going from SVT_SparseArray to ordinary R array
 */

/* Recursive. */
static int REC_unroll_SVT_into_Rarray(SEXP SVT,
		const int *dim, int ndim,
		SEXP Rarray, R_xlen_t arr_offset, R_xlen_t subarr_len)
{
	if (SVT == R_NilValue)
		return 0;

	if (ndim == 1) {
		/* 'SVT' is a leaf (i.e. a 1D SVT). */
		_expand_leaf(SVT, Rarray, arr_offset);
		return 0;
	}

	/* 'SVT' is a regular node (list). */
	int SVT_len = LENGTH(SVT);

	/* Sanity check (should never fail). */
	if (SVT_len != dim[ndim - 1])
		return -1;

	subarr_len /= SVT_len;
	for (int i = 0; i < SVT_len; i++) {
		SEXP subSVT = VECTOR_ELT(SVT, i);
		int ret = REC_unroll_SVT_into_Rarray(subSVT,
				dim, ndim - 1,
				Rarray, arr_offset, subarr_len);
		if (ret < 0)
			return -1;
		arr_offset += subarr_len;
	}
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_from_SVT_SparseArray_to_Rarray(SEXP x_dim, SEXP x_dimnames,
		SEXP x_type, SEXP x_SVT, SEXP x_na_background)
{
	SEXPTYPE Rtype = _get_and_check_Rtype_from_Rstring(x_type,
			"C_from_SVT_SparseArray_to_Rarray", "x_type");

	int x_has_NAbg = _get_and_check_na_background(x_na_background,
			"C_from_SVT_SparseArray_to_Rarray", "x_na_background");

	SEXP ans;
	if (x_has_NAbg) {
		ans = _new_RarrayNA(Rtype, x_dim, x_dimnames);
	} else {
		ans = _new_Rarray0(Rtype, x_dim, x_dimnames);
	}
	PROTECT(ans);
	int ret = REC_unroll_SVT_into_Rarray(x_SVT,
				INTEGER(x_dim), LENGTH(x_dim),
				ans, 0, XLENGTH(ans));
	UNPROTECT(1);
	if (ret < 0)
		error("SparseArray internal error in "
		      "C_from_SVT_SparseArray_to_Rarray():\n"
		      "    invalid SVT_SparseArray object");
	return ans;
}


/****************************************************************************
 * Going from ordinary R array to SVT_SparseArray
 */

/* Recursive. */
static SEXP REC_build_SVT_from_Rsubarray(
		SEXP Rarray, R_xlen_t arr_offset, R_xlen_t subarr_len,
		const int *dim, int ndim,
		SEXPTYPE ans_Rtype, int ans_na_background,
		int *warn, int *offs_buf)
{
	if (ndim == 1) {
		/* Sanity check (should never fail). */
		if (dim[0] != subarr_len)
			error("SparseArray internal error in "
			      "REC_build_SVT_from_Rsubarray():\n"
			      "    dim[0] != subarr_len");
		SEXP ans;
		if (ans_na_background) {
			ans = _make_naleaf_from_Rsubvec(
					Rarray, arr_offset, dim[0],
					offs_buf, 1);
		} else {
			ans = _make_leaf_from_Rsubvec(
					Rarray, arr_offset, dim[0],
					offs_buf, 1);
		}
		if (ans_Rtype == TYPEOF(Rarray) || ans == R_NilValue)
			return ans;
		PROTECT(ans);
		if (ans_na_background) {
			ans = _coerce_naleaf(ans, ans_Rtype, warn, offs_buf);
		} else {
			ans = _coerce_leaf(ans, ans_Rtype, warn, offs_buf);
		}
		UNPROTECT(1);
		return ans;
	}

	int SVT_len = dim[ndim - 1];  /* cannot be 0 so safe to divide below */
	subarr_len /= SVT_len;
	SEXP ans = PROTECT(NEW_LIST(SVT_len));
	int is_empty = 1;
	for (int i = 0; i < SVT_len; i++) {
		SEXP ans_elt = REC_build_SVT_from_Rsubarray(
					Rarray, arr_offset, subarr_len,
					dim, ndim - 1,
					ans_Rtype, ans_na_background,
					warn, offs_buf);
		if (ans_elt != R_NilValue) {
			PROTECT(ans_elt);
			SET_VECTOR_ELT(ans, i, ans_elt);
			UNPROTECT(1);
			is_empty = 0;
		}
		arr_offset += subarr_len;
	}
	UNPROTECT(1);
	return is_empty ? R_NilValue : ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_build_SVT_from_Rarray(SEXP x, SEXP ans_type, SEXP ans_na_background)
{
	SEXPTYPE ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (ans_Rtype == 0)
		error("invalid requested type");

	int ans_has_NAbg = _get_and_check_na_background(ans_na_background,
				"C_build_SVT_from_Rarray", "ans_na_background");

	R_xlen_t x_len = XLENGTH(x);
	if (x_len == 0)  /* means that 'any(dim(x) == 0)' is TRUE */
		return R_NilValue;

	SEXP x_dim = GET_DIM(x);  /* does not contain zeros */
	int x_ndim = LENGTH(x_dim);
	int *offs_buf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	int warn = 0;
	SEXP ans = REC_build_SVT_from_Rsubarray(x, 0, x_len,
				INTEGER(x_dim), x_ndim,
				ans_Rtype, ans_has_NAbg,
				&warn, offs_buf);
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
 * Going from SVT_SparseMatrix to [d|l|n]gCMatrix
 */

/* Returns nb of nonzero elements in column. */
static int dump_leaf_to_ix(SEXP leaf,
		SEXP sloti, SEXP slotx, int ix_offset)
{
	if (leaf == R_NilValue)
		return 0;

	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);

	/* Copy 0-based row indices from 'nzoffs' to 'sloti'. */
	memcpy(INTEGER(sloti) + ix_offset, INTEGER(nzoffs),
	       sizeof(int) * nzcount);

	if (slotx == R_NilValue)
		return nzcount;

	/* Copy 'nzvals' to 'slotx'. */
	copy_nzvals_elts_to_Rsubvec(nzvals, slotx, ix_offset, nzcount);
	return nzcount;
}

static int dump_SVT_to_CsparseMatrix_slots(SEXP x_SVT, int x_ncol,
		SEXP slotp, SEXP sloti, SEXP slotx)
{
	INTEGER(slotp)[0] = 0;
	int ix_offset = 0;
	for (int j = 0; j < x_ncol; j++) {
		SEXP leaf = VECTOR_ELT(x_SVT, j);
		int nzcount = dump_leaf_to_ix(leaf, sloti, slotx, ix_offset);
		if (nzcount < 0)
			return -1;
		ix_offset += nzcount;
		INTEGER(slotp)[j + 1] = ix_offset;
	}
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_from_SVT_SparseMatrix_to_CsparseMatrix(SEXP x_dim,
		SEXP x_type, SEXP x_SVT, SEXP as_ngCMatrix)
{
	if (LENGTH(x_dim) != 2)
		error("object to coerce to [d|l]gCMatrix "
		      "must have exactly 2 dimensions");

	R_xlen_t x_nzcount = _REC_nzcount_SVT(x_SVT, LENGTH(x_dim));
	if (x_nzcount > INT_MAX)
		error("SVT_SparseMatrix object contains too many nonzero "
		      "values to be turned into a dgCMatrix or lgCMatrix "
		      "object");

	SEXPTYPE x_Rtype = _get_and_check_Rtype_from_Rstring(x_type,
			"C_from_SVT_SparseMatrix_to_CsparseMatrix", "x_type");

	int x_ncol = INTEGER(x_dim)[1];
	SEXP sloti = PROTECT(NEW_INTEGER(x_nzcount));
	int drop_slotx = LOGICAL(as_ngCMatrix)[0];
	SEXP slotx = R_NilValue;
	if (!drop_slotx)
		slotx = PROTECT(allocVector(x_Rtype, x_nzcount));
	SEXP slotp;
	if (x_SVT == R_NilValue) {
		slotp = PROTECT(_new_Rvector0(INTSXP, (R_xlen_t) x_ncol + 1));
	} else {
		slotp = PROTECT(NEW_INTEGER(x_ncol + 1));
		int ret = dump_SVT_to_CsparseMatrix_slots(x_SVT, x_ncol,
						slotp, sloti, slotx);
		if (ret < 0) {
			UNPROTECT(3);
			error("SparseArray internal error in "
			      "C_from_SVT_SparseMatrix_to_CsparseMatrix():\n"
			      "    invalid SVT_SparseMatrix object");
		}
	}

	SEXP ans = PROTECT(NEW_LIST(3));
	SET_VECTOR_ELT(ans, 0, slotp);
	SET_VECTOR_ELT(ans, 1, sloti);
	SET_VECTOR_ELT(ans, 2, slotx);
	UNPROTECT(drop_slotx ? 3 : 4);
	return ans;
}


/****************************************************************************
 * Going from [d|l|n]gCMatrix to SVT_SparseMatrix
 */

static char get_gCMatrix_subtype(SEXP x)
{
	SEXP class_attrib = GET_CLASS(x);
	const char *class = CHAR(STRING_ELT(class_attrib, 0));
	if (strcmp(class, "dgCMatrix") == 0)
		return 'd';
	if (strcmp(class, "lgCMatrix") == 0)
		return 'l';
	if (strcmp(class, "ngCMatrix") == 0)
		return 'n';
	error("'x' must be a [d|l|n]gCMatrix object");
}

static SEXP build_leaf_from_ngCsparseMatrix_col(const int *sloti,
		R_xlen_t ix_offset, int col_nzcount, SEXPTYPE ans_Rtype)
{
	SEXP ans_nzoffs = PROTECT(NEW_INTEGER(col_nzcount));
	memcpy(INTEGER(ans_nzoffs), sloti + ix_offset,
	       sizeof(int) * col_nzcount);
	SEXP ans_nzvals = LACUNAR_MODE_IS_ON ?
		R_NilValue : PROTECT(_new_Rvector1(ans_Rtype, col_nzcount));
	SEXP ans = zip_leaf(ans_nzvals, ans_nzoffs, 0);
	UNPROTECT(LACUNAR_MODE_IS_ON ? 1 : 2);
	return ans;
}

/* Returns R_NilValue or a leaf with an nzcount <= 'col_nzcount'.

   Note that [d|l]gCMatrix objects are allowed to have zeros in their '@x'
   slot. See IMPORTANT NOTE in R/is_nonzero.R.
   So we want to make sure that these zeros don't end up in the leaf returned
   by build_leaf_from_CsparseMatrix_col(). Unfortunately, this introduces
   an additional cost to coercion from [d|l]gCMatrix to SVT_SparseMatrix.
   E.g. some quick testing seems to indicate a slowdown for this coercion
   that is (approx.) between 1.3x and 1.5x. */
static SEXP build_leaf_from_CsparseMatrix_col(const int *sloti, SEXP slotx,
		R_xlen_t ix_offset, int col_nzcount,
		SEXPTYPE ans_Rtype, int *warn, int *nzoffs_buf)
{
	/* 'slotx' can contain zeros. See above. */
	SEXP ans = _make_leaf_from_Rsubvec(slotx, ix_offset, col_nzcount,
					   nzoffs_buf, 1);
	if (ans == R_NilValue)
		return ans;

	PROTECT(ans);

	/* Replace offsets in 'ans_nzoffs' with offsets from 'sloti'. */
	SEXP ans_nzoffs = get_leaf_nzoffs(ans);
	int ans_nzcount = LENGTH(ans_nzoffs);  /* always <= 'col_nzcount' */
	_copy_selected_int_elts(sloti + ix_offset,
				INTEGER(ans_nzoffs), ans_nzcount,
				INTEGER(ans_nzoffs));
	if (ans_Rtype != TYPEOF(slotx))
		ans = _coerce_leaf(ans, ans_Rtype, warn, nzoffs_buf);

	UNPROTECT(1);
	return ans;
}

#define	GET_SLOTP_ELT(slotp, j) \
	(IS_INTEGER(slotp) ? INTEGER(slotp)[j] : REAL(slotp)[j])

/* 'slotp' must be an integer or double vector of length 'ncol' + 1.
   'slotx' can be R_NilValue. If not, 'slotx' and 'sloti' must be parallel. */
SEXP build_SVT_from_CSC(int nrow, int ncol, SEXP slotp,
		SEXP slotx, const int *sloti, int sloti_is_1based,
		SEXPTYPE ans_Rtype,
		int *order_buf, unsigned short int *rxbuf1, int *rxbuf2)
{
	if (!((IS_INTEGER(slotp) || IS_NUMERIC(slotp)) &&
	      LENGTH(slotp) == ncol + 1 &&
	      GET_SLOTP_ELT(slotp, 0) == 0))
	{
		error("SparseArray internal error in build_SVT_from_CSC():\n"
		      "    invalid 'slotp'");
	}
	/* 'slotp[ncol]' is the common length of 'slotx' and 'sloti'. */
	R_xlen_t ix_len = (R_xlen_t) GET_SLOTP_ELT(slotp, ncol);
	if (ix_len == 0)
		return R_NilValue;

	int *nzoffs_buf = (int *) R_alloc(nrow, sizeof(int));
	SEXP ans = PROTECT(NEW_LIST(ncol));
	int warn = 0;
	int is_empty = 1;
	for (int j = 0; j < ncol; j++) {
		R_xlen_t ix_offset = GET_SLOTP_ELT(slotp, j);
		int col_nzcount = GET_SLOTP_ELT(slotp, j + 1) - ix_offset;
		if (col_nzcount == 0)
			continue;
		SEXP ans_elt = slotx == R_NilValue ?
			build_leaf_from_ngCsparseMatrix_col(sloti,
					ix_offset, col_nzcount, ans_Rtype
			)
			:
			build_leaf_from_CsparseMatrix_col(sloti, slotx,
					ix_offset, col_nzcount, ans_Rtype,
					&warn, nzoffs_buf);
		if (ans_elt != R_NilValue) {
			PROTECT(ans_elt);
			if (order_buf != NULL)
				_INPLACE_order_leaf_by_nzoff(ans_elt,
							     order_buf,
							     rxbuf1, rxbuf2);
			/* We trust that 'sloti' contains no duplicates
			   within columns. If that's the case then the nzoffs
			   in 'ans_elt' should now be in strictly ascending
			   order. If we no longer want to trust 'sloti', here
			   would be a good place to check that the nzoffs have
			   no repeated values, and raise an error if they have.
			*/
			if (sloti_is_1based) {
				SEXP nzoffs = get_leaf_nzoffs(ans_elt);
				int *nzoffs_p = INTEGER(nzoffs);
				int nzcount = LENGTH(nzoffs);
				for (int k = 0; k < nzcount; k++)
					nzoffs_p[k]--;
			}
			SET_VECTOR_ELT(ans, j, ans_elt);
			UNPROTECT(1);
			is_empty = 0;
		}
	}
	if (warn)
		_CoercionWarning(warn);
	UNPROTECT(1);
	return is_empty ? R_NilValue : ans;
}

/* --- .Call ENTRY POINT ---
   'indptr' can be of type "integer" or "double".
   'indices' is expected to contain 1-based row indices. */
SEXP C_build_SVT_from_CSC(SEXP dim, SEXP indptr, SEXP data, SEXP indices,
			  SEXP indices_are_1based)
{
	if (!(IS_INTEGER(dim) && LENGTH(dim) == 2))
		error("SparseArray internal error in C_build_SVT_from_CSC():\n"
		      "    invalid 'dim'");
	int nrow = INTEGER(dim)[0];
	int ncol = INTEGER(dim)[1];
	if (!(IS_INTEGER(indices) && LENGTH(indices) == LENGTH(data)))
		error("SparseArray internal error in C_build_SVT_from_CSC():\n"
		      "    invalid 'indices'");
	int one_based = LOGICAL(indices_are_1based)[0];
	int *order_buf = NULL;
	unsigned short int *rxbuf1 = NULL;
	int *rxbuf2 = NULL;
	if (nrow >= 2) {
		order_buf = (int *) R_alloc(nrow, sizeof(int));
		rxbuf1 = (unsigned short int *)
			 R_alloc(nrow, sizeof(unsigned short int));
		rxbuf2 = (int *) R_alloc(nrow, sizeof(int));
	}
	return build_SVT_from_CSC(nrow, ncol, indptr,
				  data, INTEGER(indices), one_based,
				  TYPEOF(data), order_buf, rxbuf1, rxbuf2);
}

/* --- .Call ENTRY POINT --- */
SEXP C_build_SVT_from_CsparseMatrix(SEXP x, SEXP ans_type)
{
	char x_type = get_gCMatrix_subtype(x);
	SEXPTYPE ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (ans_Rtype == 0)
		error("invalid requested type");
	SEXP x_Dim = GET_SLOT(x, install("Dim"));
	int x_nrow = INTEGER(x_Dim)[0];
	int x_ncol = INTEGER(x_Dim)[1];
	SEXP x_slotp = GET_SLOT(x, install("p"));
	SEXP x_slotx = x_type == 'n' ? R_NilValue : GET_SLOT(x, install("x"));
	SEXP x_sloti = GET_SLOT(x, install("i"));
	return build_SVT_from_CSC(x_nrow, x_ncol, x_slotp,
				  x_slotx, INTEGER(x_sloti), 0,
				  ans_Rtype, NULL, NULL, NULL);
}


/****************************************************************************
 * Going from SVT_SparseArray to COO_SparseArray
 */

static SEXP alloc_nzvals(SEXP type, R_xlen_t n)
{
	SEXPTYPE Rtype = _get_and_check_Rtype_from_Rstring(type,
					"alloc_nzvals", "type");
	return allocVector(Rtype, n);
}

/* --- .Call ENTRY POINT --- */
SEXP C_from_SVT_SparseArray_to_COO_SparseArray(SEXP x_dim,
		SEXP x_type, SEXP x_SVT)
{
	R_xlen_t nzcount = _REC_nzcount_SVT(x_SVT, LENGTH(x_dim));
	if (nzcount > INT_MAX)
		error("SVT_SparseArray object contains too many nonzero "
		      "values to be turned into a COO_SparseArray object");
	SEXP ans_nzvals = PROTECT(alloc_nzvals(x_type, nzcount));
	SEXP ans_nzcoo = PROTECT(
		extract_nzcoo_and_nzvals_from_SVT(x_SVT,
				(int) nzcount, LENGTH(x_dim), ans_nzvals)
	);
	SEXP ans = PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(ans, 0, ans_nzcoo);
	SET_VECTOR_ELT(ans, 1, ans_nzvals);
	UNPROTECT(3);
	return ans;
}

