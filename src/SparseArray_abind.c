/****************************************************************************
 *              Combining multidimensional SparseArray objects              *
 ****************************************************************************/
#include "SparseArray_abind.h"

#include "Rvector_utils.h"
#include "leaf_vector_utils.h"  /* for _split_leaf_vector() */


static SEXP check_and_combine_object_dims(SEXP objects, int along0,
		int *dims_along)
{
	SEXP object, dim, ans_dim;
	int nb_objects, n;

	object = VECTOR_ELT(objects, 0);
	dim = GET_SLOT(object, install("dim"));
	if (along0 < 0 || along0 >= LENGTH(dim))
		error("'along' must be >= 1 and <= the number "
		      "of dimensions of the objects to bind");
	dims_along[0] = INTEGER(dim)[along0];

	ans_dim = PROTECT(duplicate(dim));
	nb_objects = LENGTH(objects);
	for (n = 1; n < nb_objects; n++) {
		object = VECTOR_ELT(objects, n);
		dim = GET_SLOT(object, install("dim"));
		if (XLENGTH(dim) != XLENGTH(ans_dim)) {
			UNPROTECT(1);
			error("all the objects to bind must have "
			      "the same number of dimensions");
		}
		INTEGER(ans_dim)[along0] +=
			dims_along[n] = INTEGER(dim)[along0];
	}
	UNPROTECT(1);
	return ans_dim;
}

static SEXP *prepare_SVTs_buf(SEXP objects, int ndim, int along0)
{
	int nb_objects, n;
	SEXP *SVTs_buf, object;

	nb_objects = LENGTH(objects);
	SVTs_buf = (SEXP *) R_alloc(nb_objects * (ndim - along0), sizeof(SEXP));
	for (n = 0; n < nb_objects; n++) {
		object = VECTOR_ELT(objects, n);
		SVTs_buf[n] = GET_SLOT(object, install("SVT"));
	}
	return SVTs_buf;
}

static int all_NULLs(SEXP *SVTs, int nb_objects)
{
	int n;

	for (n = 0; n < nb_objects; n++)
		if (SVTs[n] != R_NilValue)
			return 0;
	return 1;
}

/* All the SVTs are expected to have their last dimension (a.k.a. rightmost
   dimension, a.k.a. outermost dimension) equal to 'd'. */
static int collect_SVTs_ith_elt(SEXP *SVTs, int nb_objects, int i, int d,
		SEXP *subSVTs_buf)
{
	int n;
	SEXP SVT, subSVT;

	for (n = 0; n < nb_objects; n++) {
		SVT = SVTs[n];
		/* 'SVT' must be NULL or a list of length 'd'. */
		if (SVT == R_NilValue) {
			subSVT = R_NilValue;
		} else {
			if (!isVectorList(SVT))  // IS_LIST() is broken
				return -1;
			if (LENGTH(SVT) != d)
				return -1;
			subSVT = VECTOR_ELT(SVT, i);
		}
		subSVTs_buf[n] = subSVT;
	}
	return 0;
}

/* Each SVT is expected to be either NULL or a regular SVT node (i.e. list).
   They should never be all NULLs.
   Returns a list of length 'sum_dims_along'. */
static SEXP concatenate_SVTs(SEXP *SVTs, int nb_objects,
		const int *dims_along, int sum_dims_along)
{
	SEXP ans, SVT;
	int i1, n, SVT_len, i2;

	ans = PROTECT(NEW_LIST(sum_dims_along));

	i1 = 0;
	for (n = 0; n < nb_objects; n++) {
		SVT = SVTs[n];
		if (SVT == R_NilValue) {
			i1 += dims_along[n];
			continue;
		}
		/* Sanity check (should never fail). */
		if (!isVectorList(SVT))
			error("input object %d is an invalid SVT_SparseArray",
			      n + 1);
		SVT_len = LENGTH(SVT);
		/* Sanity check (should never fail). */
		if (SVT_len != dims_along[n])
			error("input object %d is an invalid SVT_SparseArray",
			      n + 1);
		for (i2 = 0; i2 < SVT_len; i2++, i1++)
			SET_VECTOR_ELT(ans, i1, VECTOR_ELT(SVT, i2));
	}

	UNPROTECT(1);
	/* Sanity check (should never fail). */
	if (i1 != sum_dims_along)
		error("SparseArray internal error in concatenate_SVTs():\n"
		      "    i1 != sum_dims_along");
	return ans;
}

/* Each SVT is expected to be either NULL or a "leaf vector".
   They should never be all NULLs.
   'sum_dims_along' is used for a sanity check only so is not strictly needed */
static SEXP concatenate_leaf_vectors(SEXP *SVTs, int nb_objects,
		const int *dims_along, int sum_dims_along,
		SEXPTYPE ans_Rtype,
		CopyRVectorElts_FUNType copy_Rvector_elts_FUN)
{
	SEXP SVT, ans, ans_offs, ans_vals, lv_offs, lv_vals;
	int ans_len, n, k1, offset, lv_len, k2;

	ans_len = 0;
	for (n = 0; n < nb_objects; n++) {
		SVT = SVTs[n];
		if (SVT == R_NilValue)
			continue;
		/* Sanity check (should never fail). */
		if (!isVectorList(SVT) || LENGTH(SVT) != 2)
			error("input object %d is an invalid SVT_SparseArray",
			      n + 1);
		ans_len += LENGTH(VECTOR_ELT(SVT, 0));
	}

	ans_offs = PROTECT(NEW_INTEGER(ans_len));
	ans_vals = PROTECT(allocVector(ans_Rtype, ans_len));
	k1 = offset = 0;
	for (n = 0; n < nb_objects; n++) {
		if ((SVT = SVTs[n]) != R_NilValue) {
			lv_len = _split_leaf_vector(SVT, &lv_offs, &lv_vals);
			copy_Rvector_elts_FUN(lv_vals, 0, ans_vals, k1, lv_len);
			for (k2 = 0; k2 < lv_len; k2++, k1++)
				INTEGER(ans_offs)[k1] = INTEGER(lv_offs)[k2] +
							offset;
		}
		offset += dims_along[n];
	}
	ans = _new_leaf_vector(ans_offs, ans_vals);
	UNPROTECT(2);
	/* Sanity checks (should never fail). */
	if (k1 != ans_len)
		error("SparseArray internal error in "
		      "concatenate_leaf_vectors():\n"
		      "    k1 != ans_len");
	if (offset != sum_dims_along)
		error("SparseArray internal error in "
		      "concatenate_leaf_vectors():\n"
		      "    offset != sum_dims_along");
	return ans;
}

/* Recursive.
   'along0' will always be >= 0 and < 'ndim'.
   Returns R_NilValue or a list of length 'ans_dim[ndim - 1]'. */
static SEXP REC_abind_SVTs(SEXP *SVTs, int nb_objects,
		const int *ans_dim, int ndim, int along0, const int *dims_along,
		SEXPTYPE ans_Rtype,
		CopyRVectorElts_FUNType copy_Rvector_elts_FUN)
{
	SEXP *subSVTs_buf, ans, ans_elt;
	int ans_len, is_empty, i, ret;

	if (all_NULLs(SVTs, nb_objects))
		return R_NilValue;

	if (ndim == 1) {
		/* This is the case where we are binding SVT_SparseArray
		   objects along their first dimension. */
		return concatenate_leaf_vectors(SVTs, nb_objects,
					dims_along, ans_dim[along0],
					ans_Rtype,
					copy_Rvector_elts_FUN);
	}

	if (along0 == ndim - 1)
		return concatenate_SVTs(SVTs, nb_objects,
					dims_along, ans_dim[along0]);

	subSVTs_buf = SVTs + nb_objects;
	ans_len = ans_dim[ndim - 1];
	ans = PROTECT(NEW_LIST(ans_len));
	is_empty = 1;
	for (i = 0; i < ans_len; i++) {
		ret = collect_SVTs_ith_elt(SVTs, nb_objects, i, ans_len,
					   subSVTs_buf);
		if (ret < 0) {
			UNPROTECT(1);
			error("SparseArray internal error in "
			      "REC_abind_SVTs():\n"
			      "    collect_SVTs_ith_elt() returned an error");
		}
		ans_elt = REC_abind_SVTs(subSVTs_buf, nb_objects,
				ans_dim, ndim - 1, along0, dims_along,
				ans_Rtype, copy_Rvector_elts_FUN);
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

/* --- .Call ENTRY POINT --- */
SEXP C_abind_SVT_SparseArray_objects(SEXP objects, SEXP along, SEXP ans_type)
{
	SEXPTYPE ans_Rtype;
	CopyRVectorElts_FUNType copy_Rvector_elts_FUN;
	int along0, nb_objects, *dims_along, ndim;
	SEXP ans_dim, *SVTs_buf, ans_SVT, ans;

	if (!isVectorList(objects))  // IS_LIST() is broken
		error("'objects' must be a list of SVT_SparseArray objects");

	ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	copy_Rvector_elts_FUN = _select_copy_Rvector_elts_FUN(ans_Rtype);
	if (copy_Rvector_elts_FUN == NULL)
		error("invalid requested type");

	if (!IS_INTEGER(along) || XLENGTH(along) != 1)
		error("'along' must be a single positive integer");
	along0 = INTEGER(along)[0] - 1;

	nb_objects = LENGTH(objects);
	if (nb_objects == 0)
		error("'objects' cannot be an empty list");

	dims_along = (int *) R_alloc(nb_objects, sizeof(int));
	ans_dim = PROTECT(
		check_and_combine_object_dims(objects, along0, dims_along)
	);
	ndim = LENGTH(ans_dim);

	SVTs_buf = prepare_SVTs_buf(objects, ndim, along0);
	ans_SVT = REC_abind_SVTs(SVTs_buf, nb_objects,
			INTEGER(ans_dim), ndim, along0, dims_along,
			ans_Rtype, copy_Rvector_elts_FUN);
	if (ans_SVT != R_NilValue)
		PROTECT(ans_SVT);

	ans = PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(ans, 0, ans_dim);
	if (ans_SVT != R_NilValue) {
		SET_VECTOR_ELT(ans, 1, ans_SVT);
		UNPROTECT(1);
	}
	UNPROTECT(2);
	return ans;
}

