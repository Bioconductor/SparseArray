/****************************************************************************
 *          Drop/add ineffective dims from/to a SparseArray object          *
 ****************************************************************************/
#include "SparseArray_ineffective_dims.h"

#include "Rvector_utils.h"
#include "leaf_vector_utils.h"  /* for _split_leaf_vector() */

#include <string.h>  /* for memset() */


/****************************************************************************
 * Low-level handling of the "dimension tuner" ('dim_tuner' argument)
 *
 * The 'dim_tuner' argument is a logical vector (possibly with NAs) that
 * indicates what to do with the dimensions in 'dims'. It can contain any
 * number of NAs anywhere, but it **must** have exactly 'ndim' non-NA
 * values, that is, 'ndim' TRUE/FALSE values.
 * Meaning of values in 'dim_tuner':
 *   o TRUE: Keep corresponding dimension in 'dims'.
 *   o FALSE: Drop corresponding dimension in 'dims'. This is allowed only
 *     if the dimension to drop is ineffective (i.e. has an extent of 1).
 *   o NA: Add ineffective dim.
 *
 * Additionally, 'dim_tuner' must contain at least one TRUE value (i.e.
 * at least one of the original dimension must be retained).
 *
 * The following patterns:
 *   - NA followed by a FALSE,
 *   - FALSE followed by an NA,
 * are both equivalent to a single TRUE in 'dim_tuner'. 'dim_tuner'
 * is considered to be "normalized" if it doesn't contain these patterns.
 */

static int dim_tuner_is_normalized(const int *dim_tuner, int dim_tuner_len)
{
	int prev_akd, i, akd;

	prev_akd = dim_tuner[0];  // 'dim_tuner_len' is guaranteed to be >= 1
	for (i = 1; i < dim_tuner_len; i++) {
		akd = dim_tuner[i];
		if (akd == NA_INTEGER) {
			if (prev_akd != NA_INTEGER && !prev_akd)
				return 0;
		} else if (!akd) {
			if (prev_akd == NA_INTEGER)
				return 0;
		}
		prev_akd = akd;
	}
	return 1;
}

/* Return the "new" number of dimensions i.e. the number of dims that we
   will get after tuning the current vector of dimensions. */
static int validate_dim_tuner(const int *dim_tuner, int dim_tuner_len,
		const int *dims, int ndim,
		int *cumallTRUE, int *cumallFALSE)
{
	int along1, along2, nkept, i, akd;

	if (cumallTRUE != NULL)
		memset(cumallTRUE, 0, sizeof(int) * ndim);
	if (cumallFALSE != NULL)
		memset(cumallFALSE, 0, sizeof(int) * ndim);
	along1 = along2 = nkept = 0;
	for (i = 0; i < dim_tuner_len; i++) {
		/* 'akd' can be NA (add ineffective dim), TRUE (keep dim),
		   or FALSE (drop ineffective dim). */
		akd = dim_tuner[i];
		if (akd == NA_INTEGER) {
			/* Add ineffective dimension. */
			along2++;
			continue;
		}
		if (along1 >= ndim)
			error("SparseArray internal error in "
			      "validate_dim_tuner():\n"
			      "    number of non-NA values in 'dim_tuner' "
			      "is > 'length(dim(x))'");
		if (akd) {
			/* Keep dimension. */
			if (cumallTRUE != NULL &&
			    i == along1 && (i == 0 || cumallTRUE[i - 1]))
				cumallTRUE[i] = 1;
			along2++;
			nkept++;
			along1++;
			continue;
		}
		/* Drop ineffective dimension. */
		if (dims[along1] != 1)
			error("SparseArray internal error in "
			      "validate_dim_tuner():\n"
			      "    'dim_tuner[%d]' (= FALSE) is "
			      "mapped to 'dim(x)[%d]' (= %d)\n"
			      "    which cannot be dropped",
			      i + 1, along1 + 1, dims[along1]);
		if (cumallFALSE != NULL &&
		    i == along1 && (i == 0 || cumallFALSE[i - 1]))
			cumallFALSE[i] = 1;
		along1++;
	}
	if (along1 < ndim)
		error("SparseArray internal error in "
		      "validate_dim_tuner():\n"
		      "    number of non-NA values in 'dim_tuner' "
		      "is < 'length(dim(x))'");
	if (nkept == 0)
		error("SparseArray internal error in "
		      "validate_dim_tuner():\n"
		      "    'dim_tuner' must contain at least one TRUE");
	return along2;
}


/****************************************************************************
 * C_tune_dims() and C_tune_dimnames()
 */

/* 'dim_names' must be a character vector (or R_NilValue) that contains the
   names on the vector of dimensions. This is NOT the same as the dimnames! */
static SEXP tune_dims(const int *dims, SEXP dim_names,
		const int *dim_tuner, int dim_tuner_len, int ans_len)
{
	SEXP ans, ans_names;
	int along1, along2, i, akd;

	ans = PROTECT(NEW_INTEGER(ans_len));
	if (dim_names != R_NilValue)
		ans_names = PROTECT(NEW_CHARACTER(ans_len));
	along1 = along2 = 0;
	for (i = 0; i < dim_tuner_len; i++) {
		akd = dim_tuner[i];
		if (akd == NA_INTEGER) {
			/* Add ineffective dimension. */
			INTEGER(ans)[along2] = 1;
			along2++;
			continue;
		}
		if (akd) {
			/* Keep dimension. */
			INTEGER(ans)[along2] = dims[along1];
			if (dim_names != R_NilValue)
				SET_STRING_ELT(ans_names, along2,
					       STRING_ELT(dim_names, along1));
			along2++;
		}
		along1++;
	}
	if (dim_names != R_NilValue) {
		SET_NAMES(ans, ans_names);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

static SEXP tune_dimnames(SEXP dimnames,
		const int *dim_tuner, int dim_tuner_len, int ans_len)
{
	int along1, along2, i, akd;
	SEXP ans;

	ans = PROTECT(NEW_LIST(ans_len));
	along1 = along2 = 0;
	for (i = 0; i < dim_tuner_len; i++) {
		akd = dim_tuner[i];
		if (akd == NA_INTEGER) {
			/* Add ineffective dimension. */
			along2++;
			continue;
		}
		if (akd) {
			/* Keep dimension. */
			SET_VECTOR_ELT(ans, along2,
				       VECTOR_ELT(dimnames, along1));
			along2++;
		}
		along1++;
	}
	UNPROTECT(1);
	return ans;
}

/* Return the length of the tuned 'dimnames' or 0 if the tuning retains no
   dimnames. */
static int compute_tuned_dimnames_length(SEXP dimnames,
		const int *dim_tuner, int dim_tuner_len)
{
	int ndim, along1, along2, any_retained, i, akd;

	if (dimnames == R_NilValue)
		return 0;
	ndim = LENGTH(dimnames);
	along1 = along2 = any_retained = 0;
	for (i = 0; i < dim_tuner_len; i++) {
		akd = dim_tuner[i];
		if (akd == NA_INTEGER) {
			along2++;
			continue;
		}
		if (along1 >= ndim)
			error("SparseArray internal error in "
			      "compute_tuned_dimnames_length():\n"
			      "    number of non-NA values in 'dim_tuner' "
			      "is > 'length(dim(x))'");
		if (akd) {
			if (VECTOR_ELT(dimnames, along1) != R_NilValue)
				any_retained = 1;
			along2++;
		}
		along1++;
	}
	return any_retained ? along2 : 0;
}

/* --- .Call ENTRY POINT ---
   Unlike C_tune_SVT_dims() below in this file, C_tune_dims() accepts
   a 'dim_tuner' that is not normalized. */
SEXP C_tune_dims(SEXP dim, SEXP dim_tuner)
{
	int ndim, dim_tuner_len, ans_len;
	const int *dims, *dim_tuner_p;

	ndim = LENGTH(dim);
	dims = INTEGER(dim);
	dim_tuner_len = LENGTH(dim_tuner);
	dim_tuner_p = LOGICAL(dim_tuner);
	ans_len = validate_dim_tuner(dim_tuner_p, dim_tuner_len,
				      dims, ndim,
				      NULL, NULL);
	return tune_dims(dims, GET_NAMES(dim),
			 dim_tuner_p, dim_tuner_len, ans_len);
}

/* --- .Call ENTRY POINT ---
   Unlike C_tune_SVT_dims() below in this file, C_tune_dimnames() accepts
   a 'dim_tuner' that is not normalized. */
SEXP C_tune_dimnames(SEXP dimnames, SEXP dim_tuner)
{
	int dim_tuner_len, ans_len;
	const int *dim_tuner_p;

	dim_tuner_len = LENGTH(dim_tuner);
	dim_tuner_p = LOGICAL(dim_tuner);
	ans_len = compute_tuned_dimnames_length(dimnames,
						dim_tuner_p, dim_tuner_len);
	if (ans_len == 0)
		return R_NilValue;
	return tune_dimnames(dimnames, dim_tuner_p, dim_tuner_len, ans_len);
}


/****************************************************************************
 * Go back and forth between a "leaf vector" and a 1xN SVT
 */

/* 'lv' is assumed to be a "leaf vector" that represents a sparse vector
   of length N. Turn it into an SVT that represents a 1xN matrix. */
static SEXP make_1xN_SVT_from_lv(SEXP lv)
{
	error("make_1xN_SVT_from_lv() not ready yet!");
	return R_NilValue;
}

/* 'SVT' is assumed to represent a 1xN matrix. Turn it into a "leaf vector"
   that represents a sparse vector of length N. */
static SEXP make_lv_from_1xN_SVT(SEXP SVT,
		SEXPTYPE Rtype, CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	int SVT_len, ans_len, i, lv_len;
	SEXP subSVT, ans_offs, ans_vals, lv_vals, ans;

	SVT_len = LENGTH(SVT);
	ans_len = 0;
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue)
			continue;
		/* 'subSVT' is a "leaf vector" of length 1. */
		lv_len = LENGTH(VECTOR_ELT(subSVT, 0));
		if (lv_len != 1)
			error("SparseArray internal error in "
			      "make_lv_from_1xN_SVT():\n"
			      "    lv_len != 1");
		ans_len++;
	}
	if (ans_len == 0)
		error("SparseArray internal error in "
		      "make_lv_from_1xN_SVT():\n"
		      "    ans_len == 0");
	ans_offs = PROTECT(NEW_INTEGER(ans_len));
	ans_vals = PROTECT(allocVector(Rtype, ans_len));
	ans_len = 0;
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue)
			continue;
		/* 'subSVT' is a "leaf vector" of length 1. */
		INTEGER(ans_offs)[ans_len] = i;
		lv_vals = VECTOR_ELT(subSVT, 1);
		copy_Rvector_elt_FUN(lv_vals, 0, ans_vals, ans_len);
		ans_len++;
	}
	ans = _new_leaf_vector(ans_offs, ans_vals);
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * C_tune_SVT_dims()
 */

/* 'lv' is assumed to represent a sparse vector of length N.
   unroll_lv_as_SVT() turns it into an SVT that represents a
   1x1x..xN array. */
static SEXP unroll_lv_as_SVT(SEXP lv, int N, int ndim_to_insert)
{
	SEXP ans;

	/* Turn 'lv' into an SVT that represents an 1xN matrix. */
	ans = PROTECT(make_1xN_SVT_from_lv(lv));
	/* Insert 'ndim_to_insert' additional inner ineffective dimensions. */
	return R_NilValue;
}

/* 'SVT' is assumed to represent a 1x1x..xN array.
   More precisely: 'ndim' is assumed to be >= 2. Except maybe for its
   outermost dimension, all the dimensions in 'SVT' are assumed to be
   ineffective.
   roll_SVT_into_lv() turns 'SVT' into a "leaf vector" that represents a
   sparse vector of length N. */
static SEXP roll_SVT_into_lv(SEXP SVT, const int *dims, int ndim)
{
	error("roll_SVT_into_lv() not ready yet!");
	return R_NilValue;
}

/* Assumes that 'dim_tuner' is normalized.
   Recursive. */
static SEXP REC_tune_SVT_dims(SEXP SVT, const int *dims, int ndim,
		const int *dim_tuner, int dim_tuner_len,
		const int *cumallTRUE, const int *cumallFALSE)
{
	int akd, ans_len, i;
	SEXP ans_elt, ans, subSVT;

	if (SVT == R_NilValue || dim_tuner_len == ndim && cumallTRUE[ndim - 1])
		return SVT;

	akd = dim_tuner[dim_tuner_len - 1];
	if (akd == NA_INTEGER) {
		/* Add ineffective dimension (as outermost dimension). */
		ans_elt = PROTECT(
			REC_tune_SVT_dims(SVT, dims, ndim,
					  dim_tuner, dim_tuner_len - 1,
					  cumallTRUE, cumallFALSE)
		);
		ans = PROTECT(NEW_LIST(1));
		SET_VECTOR_ELT(ans, 0, ans_elt);
		UNPROTECT(1);
		return ans;
	}
	if (akd) {
		/* Keep dimension. */
		if (ndim == 1) {
			/* 'dim_tuner[dim_tuner_len - 1]' is a TRUE, preceded
			   by NAs only. 'SVT' is a "leaf vector". */
			return unroll_lv_as_SVT(SVT, dims[0],
						dim_tuner_len - 2);
		}
		if (dim_tuner_len == ndim && cumallFALSE[ndim - 2]) {
			/* 'dim_tuner[dim_tuner_len - 1]' is a TRUE, preceded
			   by FALSEs only. Return a "leaf vector". */
			return roll_SVT_into_lv(SVT, dims, ndim);
		}
		ans_len = dims[ndim - 1];
		ans = PROTECT(NEW_LIST(ans_len));
		for (i = 0; i < ans_len; i++) {
			subSVT = VECTOR_ELT(SVT, i);
			ans_elt = PROTECT(
				REC_tune_SVT_dims(subSVT,
						  dims, ndim - 1,
						  dim_tuner, dim_tuner_len - 1,
						  cumallTRUE, cumallFALSE)
			);
			SET_VECTOR_ELT(ans, i, ans_elt);
			UNPROTECT(1);
		}
		UNPROTECT(1);
		return ans;
	}
	/* Drop ineffective dimension.
	   Because 'dim_tuner' is normalized, it's guaranteed to contain
	   at least one TRUE before the FALSE found at the 'dim_tuner_len - 1'
	   position. Furthermore, the closest TRUE (i.e. highest TRUE position
	   that is < 'dim_tuner_len - 1') is guaranteed to be separated from
	   the FALSE at 'dim_tuner_len - 1' by nothing but other FALSEs.
	   In particular, this means that 'ndim' is guaranteed to be >= 2
	   so 'SVT' cannot be a "leaf vector". */
	return REC_tune_SVT_dims(VECTOR_ELT(SVT, 0), dims, ndim - 1,
				 dim_tuner, dim_tuner_len - 1,
				 cumallTRUE, cumallFALSE);
}

/* --- .Call ENTRY POINT ---
   See 'Low-level handling of the "dimension tuner"' at the top of this
   file for a description of the 'dim_tuner' argument. */
SEXP C_tune_SVT_dims(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP dim_tuner)
{
	int ndim, dim_tuner_len, *cumallTRUE, *cumallFALSE;
	const int *dims, *dim_tuner_p;

	/* Make sure that: 1 <= ndim <= dim_tuner_len. */
	ndim = LENGTH(x_dim);
	if (ndim == 0)
		error("SparseArray internal error in "
		      "C_tune_SVT_dims():\n"
		      "    'dim(x)' cannot be empty");
	dim_tuner_len = LENGTH(dim_tuner);
	if (dim_tuner_len < ndim)
		error("SparseArray internal error in "
		      "C_tune_SVT_dims():\n"
		      "    length(dim_tuner) < length(dim(x))");

	/* REC_tune_SVT_dims() assumes that 'dim_tuner' is normalized.
	   Note that we have no use case for a 'dim_tuner' that is not
	   normalized at the moment. */
	dim_tuner_p = LOGICAL(dim_tuner);
	if (!dim_tuner_is_normalized(dim_tuner_p, dim_tuner_len))
		error("SparseArray internal error in "
		      "C_tune_SVT_dims():\n"
		      "    'dim_tuner' is not normalized");

	dims = INTEGER(x_dim);
	cumallTRUE = (int *) R_alloc(ndim, sizeof(int));
	cumallFALSE = (int *) R_alloc(ndim, sizeof(int));
	validate_dim_tuner(dim_tuner_p, dim_tuner_len,
			   dims, ndim,
			   cumallTRUE, cumallFALSE);

	/* Compute 'ans_SVT'. */
	return REC_tune_SVT_dims(x_SVT, dims, ndim,
				 dim_tuner_p, dim_tuner_len,
				 cumallTRUE, cumallFALSE);
}


/****************************************************************************
 * C_drop_SVT_ineffective_dims()
 */

static int count_effective_dims(SEXP x_dim)
{
	int count, x_ndim, along, d;

	count = 0;
	x_ndim = LENGTH(x_dim);
	for (along = 0; along < x_ndim; along++) {
		d = INTEGER(x_dim)[along];
		if (d != 1)
			count++;
	}
	return count;
}

static SEXP drop_ineffective_dims(SEXP x_dim)
{
	int ans_ndim, x_ndim, along1, along2, d;
	SEXP ans_dim;

	ans_ndim = count_effective_dims(x_dim);
	if (ans_ndim == 0)
		return ScalarInteger(1);

	ans_dim = PROTECT(NEW_INTEGER(ans_ndim));
	along1 = 0;
	x_ndim = LENGTH(x_dim);
	for (along2 = 0; along2 < x_ndim; along2++) {
		d = INTEGER(x_dim)[along2];
		if (d != 1) {
			INTEGER(ans_dim)[along1] = d;
			along1++;
		}
	}
	UNPROTECT(1);
	return ans_dim;
}

static SEXP drop_dimnames_along_ineffective_dims(SEXP x_dimnames, SEXP x_dim)
{
	int ans_ndim, x_ndim, along1, along2, d;
	SEXP ans_dimnames;

	ans_ndim = count_effective_dims(x_dim);
	if (ans_ndim == 0)
		return NEW_LIST(1);

	ans_dimnames = PROTECT(NEW_LIST(ans_ndim));
	along1 = 0;
	x_ndim = LENGTH(x_dim);
	for (along2 = 0; along2 < x_ndim; along2++) {
		d = INTEGER(x_dim)[along2];
		if (d != 1) {
			SET_VECTOR_ELT(ans_dimnames, along1,
				       VECTOR_ELT(x_dimnames, along2));
			along1++;
		}
	}
	UNPROTECT(1);
	return ans_dimnames;
}

/* Recursive. */
static SEXP REC_remove_SVT_unary_nodes(SEXP SVT, const int *dim, int ndim,
		const int *all_effective_dim,
		const int *any_effective_dim,
		SEXPTYPE Rtype, CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	int SVT_len, i;
	SEXP subSVT, ans, ans_elt;

	if (SVT == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		return SVT;
	}

	/* 'SVT' is a regular node (list). */

	/* Sanity check (should never fail). */
	if (!isVectorList(SVT))  // IS_LIST() is broken
		 error("SparseArray internal error in "
		       "REC_remove_SVT_unary_nodes():\n"
		       "    SVT node is not a list");

	SVT_len = LENGTH(SVT);

	/* Sanity check (should never fail). */
	if (SVT_len != dim[ndim - 1])
		 error("SparseArray internal error in "
		       "REC_remove_SVT_unary_nodes():\n"
		       "    SVT_len != dim[ndim - 1]");

	if (SVT_len == 1)
		return REC_remove_SVT_unary_nodes(
				VECTOR_ELT(SVT, 0), dim, ndim - 1,
				all_effective_dim, any_effective_dim,
				Rtype, copy_Rvector_elt_FUN);

	if (all_effective_dim[ndim - 2])
		return SVT;  /* nothing else to do */

	ans = PROTECT(NEW_LIST(SVT_len));
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue)
			continue;
		ans_elt = PROTECT(
			REC_remove_SVT_unary_nodes(
				subSVT, dim, ndim - 1,
				all_effective_dim, any_effective_dim,
				Rtype, copy_Rvector_elt_FUN)
		);
		SET_VECTOR_ELT(ans, i, ans_elt);
		UNPROTECT(1);
	}

	if (any_effective_dim[ndim - 2]) {
		UNPROTECT(1);
		return ans;  /* nothing else to do */
	}

	/* 'ans' represents a 1xN matrix. Turn it into a "leaf vector". */
	ans = make_lv_from_1xN_SVT(ans, Rtype, copy_Rvector_elt_FUN);
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_drop_SVT_ineffective_dims(SEXP x_dim, SEXP x_dimnames,
		SEXP x_type, SEXP x_SVT)
{
	SEXPTYPE Rtype;
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;
	int x_ndim, ok1, ok2, along, d;
	int *all_effective_dim, *any_effective_dim;
	SEXP ans_dim, ans_dimnames, ans_SVT, ans;

	Rtype = _get_Rtype_from_Rstring(x_type);
	copy_Rvector_elt_FUN = _select_copy_Rvector_elt_FUN(Rtype);
	if (copy_Rvector_elt_FUN == NULL)
		error("SparseArray internal error in "
		      "C_drop_SVT_ineffective_dims():\n"
		      "    SVT_SparseArray object has invalid type");

	x_ndim = LENGTH(x_dim);
	all_effective_dim = (int *) R_alloc(x_ndim, sizeof(int));
	any_effective_dim = (int *) R_alloc(x_ndim, sizeof(int));
	ok1 = 1;
	ok2 = 0;
	for (along = 0; along < x_ndim; along++) {
		d = INTEGER(x_dim)[along];
		if (d == 1)
			ok1 = 0;
		else
			ok2 = 1;
		all_effective_dim[along] = ok1;
		any_effective_dim[along] = ok2;
	}

	ans_dim = PROTECT(drop_ineffective_dims(x_dim));
	ans_dimnames = PROTECT(
		drop_dimnames_along_ineffective_dims(x_dimnames, x_dim)
	);

	if (x_SVT != R_NilValue) {
		ans_SVT = PROTECT(
			REC_remove_SVT_unary_nodes(x_SVT,
					INTEGER(x_dim), x_ndim,
					all_effective_dim,
					any_effective_dim,
					Rtype, copy_Rvector_elt_FUN)
		);
	}

	ans = PROTECT(NEW_LIST(3));
	SET_VECTOR_ELT(ans, 0, ans_dim);
	SET_VECTOR_ELT(ans, 1, ans_dimnames);
	if (x_SVT != R_NilValue) {
		SET_VECTOR_ELT(ans, 2, ans_SVT);
		UNPROTECT(1);
	}
	UNPROTECT(3);
	return ans;
}

