/****************************************************************************
 *          Drop/add ineffective dims from/to a SparseArray object          *
 ****************************************************************************/
#include "SparseArray_ineffective_dims.h"

#include "Rvector_utils.h"
#include "leaf_vector_utils.h"  /* for _split_leaf_vector() */


/****************************************************************************
 * C_select_SVT_dims()
 */

/* 'dim_selector' is a logical vector (possibly with NAs) that indicates
   what to do with the dimensions in 'dims'. It can contain any number of
   NAs anywhere but it **must** have exactly 'ndim' non-NA values, that is,
   'ndim' TRUE/FALSE values. Meaning of values in 'dim_selector':
   - TRUE: Keep corresponding dimension in 'dims'.
   - FALSE: Drop corresponding dimension in 'dims'. This is allowed only if
     the dimension to drop is ineffective (i.e. has an extend of 1).
   - NA: Add ineffective dim.
   Additionally, 'dim_selector' must contain at least one TRUE value (i.e.
   at least one of the original dimension must be retained). */
static int count_selected_dims(const int *dim_selector, int nselector,
		const int *dims, int ndim)
{
	int along1, along2, nkept, i, akd;

	along1 = along2 = nkept = 0;
	for (i = 0; i < nselector; i++) {
		/* 'akd' can be NA (add ineffective dim), TRUE (keep dim),
		   or FALSE (drop ineffective dim). */
		akd = dim_selector[i];
		if (akd == NA_INTEGER) {
			/* Add ineffective dimension. */
			along2++;
			continue;
		}
		if (along1 >= ndim)
			error("SparseArray internal error in "
			      "count_selected_dims():\n"
			      "    number of non-NA values in 'dim_selector' "
			      "is > 'length(dim(x))'");
		if (akd) {
			/* Keep dimension. */
			along2++;
			nkept++;
			along1++;
			continue;
		}
		/* Drop ineffective dimension. */
		if (dims[along1] != 1)
			error("SparseArray internal error in "
			      "count_selected_dims():\n"
			      "    'dim_selector[%d]' (= FALSE) is "
			      "mapped to 'dim(x)[%d]' (= %d)\n"
			      "    which cannot be dropped",
			      i + 1, along1 + 1, dims[along1]);
		along1++;
	}
	if (along1 < ndim)
		error("SparseArray internal error in "
		      "count_selected_dims():\n"
		      "    number of non-NA values in 'dim_selector' "
		      "is < 'length(dim(x))'");
	if (nkept == 0)
		error("SparseArray internal error in "
		      "count_selected_dims():\n"
		      "    'dim_selector' contains no TRUE values");
	return along2;
}

static SEXP select_dims(const int *dims,
			const int *dim_selector, int nselector, int ans_ndim)
{
	int along1, along2, i, akd;
	SEXP ans_dim;

	ans_dim = PROTECT(NEW_INTEGER(ans_ndim));
	along1 = along2 = 0;
	for (i = 0; i < nselector; i++) {
		akd = dim_selector[i];
		if (akd == NA_INTEGER) {
			/* Add ineffective dimension. */
			INTEGER(ans_dim)[along2++] = 1;
			continue;
		}
		if (akd) {
			/* Keep dimension. */
			INTEGER(ans_dim)[along2++] = dims[along1++];
			continue;
		}
		/* Drop ineffective dimension. */
		along1++;
	}
	UNPROTECT(1);
	return ans_dim;
}

static SEXP select_dimnames(SEXP dimnames,
		const int *dim_selector, int nselector, int ans_ndim)
{
	int along1, along2, any_retained, i, akd;
	SEXP ans_dimnames;

	if (dimnames == R_NilValue)
		return R_NilValue;

	/* 1st pass on 'dimnames': do we need to retain any dimnames at all? */
	along1 = any_retained = 0;
	for (i = 0; i < nselector; i++) {
		akd = dim_selector[i];
		if (akd == NA_INTEGER)
			continue;
		if (akd && VECTOR_ELT(dimnames, along1) != R_NilValue) {
			any_retained = 1;
			break;
		}
		along1++;
	}
	if (!any_retained)
		return R_NilValue;

	/* 2nd pass on 'dimnames': collect retained dimnames. */
	ans_dimnames = PROTECT(NEW_LIST(ans_ndim));
	along1 = along2 = 0;
	for (i = 0; i < nselector; i++) {
		akd = dim_selector[i];
		if (akd == NA_INTEGER) {
			/* Add ineffective dimension. */
			along2++;
			continue;
		}
		if (akd) {
			/* Keep dimension. */
			SET_VECTOR_ELT(ans_dimnames, along2++,
				       VECTOR_ELT(dimnames, along1++));
			continue;
		}
		/* Drop ineffective dimension. */
		along1++;
	}
	UNPROTECT(1);
	return ans_dimnames;
}

/* Recursive. */
static SEXP REC_select_SVT_dims(SEXP SVT, const int *dims, int ndim,
		const int *dim_selector, int nselector)
{
	if (SVT == R_NilValue)
		return R_NilValue;
	error("REC_select_SVT_dims() not ready yet!");
	return R_NilValue;
}

/* --- .Call ENTRY POINT ---
   See count_selected_dims() above in this file for a description of
   the 'dim_selector' argument. */
SEXP C_select_SVT_dims(SEXP x_dim, SEXP x_dimnames,
		SEXP x_type, SEXP x_SVT, SEXP dim_selector)
{
	const int *dims;
	int ndim, nselector, ans_ndim;
	SEXP ans_dim, ans_dimnames, ans_SVT, ans;

	dims = INTEGER(x_dim);
	ndim = LENGTH(x_dim);
	nselector = LENGTH(dim_selector);
	ans_ndim = count_selected_dims(LOGICAL(dim_selector), nselector,
				       dims, ndim);

	/* Compute 'ans_dim'. */
	ans_dim = PROTECT(
		select_dims(dims, LOGICAL(dim_selector), nselector, ans_ndim)
	);

	/* Compute 'ans_dimnames'. */
	ans_dimnames = select_dimnames(x_dimnames,
				       LOGICAL(dim_selector), nselector,
				       ans_ndim);
	if (ans_dimnames != R_NilValue)
		PROTECT(ans_dimnames);

	/* Compute 'ans_SVT'. */
	ans_SVT = REC_select_SVT_dims(x_SVT, dims, ndim,
				      LOGICAL(dim_selector), nselector);
	if (ans_SVT != R_NilValue)
		PROTECT(ans_SVT);

	/* Assemble and return 'ans'. */
	ans = PROTECT(NEW_LIST(3));
	SET_VECTOR_ELT(ans, 0, ans_dim);
	SET_VECTOR_ELT(ans, 1, ans_dimnames);
	SET_VECTOR_ELT(ans, 2, ans_SVT);
	UNPROTECT(2 + (ans_dimnames != R_NilValue) + (x_SVT != R_NilValue));
	return ans;
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

/* 'SVT' must be a depth-1 tree representing a 1xN matrix.
   Turn it into a "leaf vector". */
static SEXP make_leaf_vector_from_depth_one_SVT(SEXP SVT,
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
			      "make_leaf_vector_from_depth_one_SVT():\n"
			      "    lv_len != 1");
		ans_len++;
	}
	if (ans_len == 0)
		error("SparseArray internal error in "
		      "make_leaf_vector_from_depth_one_SVT():\n"
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
	ans = make_leaf_vector_from_depth_one_SVT(ans, Rtype,
						  copy_Rvector_elt_FUN);
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

