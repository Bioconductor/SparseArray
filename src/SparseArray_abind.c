/****************************************************************************
 ****************************************************************************
 **									   **
 **             Combining multidimensional SparseArray objects             **
 **									   **
 ****************************************************************************
 ****************************************************************************/
#include "SparseArray_abind.h"

#include "Rvector_utils.h"
#include "leaf_utils.h"


static SEXP check_and_combine_object_dims(SEXP objects, int along0,
		int *dims_along)
{
	SEXP object = VECTOR_ELT(objects, 0);
	SEXP dim = GET_SLOT(object, install("dim"));
	if (along0 < 0 || along0 >= LENGTH(dim))
		error("'along' must be >= 1 and <= the number "
		      "of dimensions of the objects to bind");
	dims_along[0] = INTEGER(dim)[along0];

	SEXP ans_dim = PROTECT(duplicate(dim));
	int nb_objects = LENGTH(objects);
	for (int n = 1; n < nb_objects; n++) {
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

static SEXP *prepare_SVTs_buf(SEXP objects, SEXP SVTslotname,
			      int ndim, int along0)
{
	if (!(IS_CHARACTER(SVTslotname) && LENGTH(SVTslotname) == 1))
		error("'SVTslotname' must be a single string");
	SVTslotname = STRING_ELT(SVTslotname, 0);
	if (SVTslotname == NA_STRING)
		error("'SVTslotname' cannot be NA");
	const char *slotname = CHAR(SVTslotname);

	int nb_objects = LENGTH(objects);
	SEXP *SVTs_buf = (SEXP *)
		R_alloc(nb_objects * (ndim - along0), sizeof(SEXP));
	for (int n = 0; n < nb_objects; n++) {
		SEXP object = VECTOR_ELT(objects, n);
		SVTs_buf[n] = GET_SLOT(object, install(slotname));
	}
	return SVTs_buf;
}

static int all_NULLs(SEXP *SVTs, int nb_objects)
{
	for (int n = 0; n < nb_objects; n++)
		if (SVTs[n] != R_NilValue)
			return 0;
	return 1;
}

/* All the SVTs are expected to have their last dimension (a.k.a. rightmost
   dimension, a.k.a. outermost dimension) equal to 'd'. */
static int collect_SVTs_ith_elt(SEXP *SVTs, int nb_objects, int i, int d,
		SEXP *subSVTs_buf)
{
	for (int n = 0; n < nb_objects; n++) {
		SEXP SVT = SVTs[n], subSVT;
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
	SEXP ans = PROTECT(NEW_LIST(sum_dims_along));

	int i1 = 0;
	for (int n = 0; n < nb_objects; n++) {
		SEXP SVT = SVTs[n];
		if (SVT == R_NilValue) {
			i1 += dims_along[n];
			continue;
		}
		/* Sanity check (should never fail). */
		if (!isVectorList(SVT))  // IS_LIST() is broken
			error("input object %d is an invalid SVT_SparseArray",
			      n + 1);
		int SVT_len = LENGTH(SVT);
		/* Sanity check (should never fail). */
		if (SVT_len != dims_along[n])
			error("input object %d is an invalid SVT_SparseArray",
			      n + 1);
		for (int i2 = 0; i2 < SVT_len; i2++, i1++)
			SET_VECTOR_ELT(ans, i1, VECTOR_ELT(SVT, i2));
	}

	UNPROTECT(1);
	/* Sanity check (should never fail). */
	if (i1 != sum_dims_along)
		error("SparseArray internal error in concatenate_SVTs():\n"
		      "    i1 != sum_dims_along");
	return ans;
}

/* The leaves in 'leaves' should never be all NULLs.
   'sum_dims_along' is used for a sanity check only so is not strictly needed */
static SEXP concatenate_leaves(SEXP *leaves, int nb_objects,
		const int *dims_along, int sum_dims_along,
		SEXPTYPE ans_Rtype)
{
	int ans_nzcount = 0, ans_is_lacunar = 1;
	for (int n = 0; n < nb_objects; n++) {
		SEXP leaf = leaves[n];
		if (leaf == R_NilValue)
			continue;
		SEXP nzvals, nzoffs;
		int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
		ans_nzcount += nzcount;
		if (nzvals != R_NilValue &&
		    !_all_Rvector_elts_equal_one(nzvals))
			ans_is_lacunar = 0;
	}

	SEXP ans_nzvals;
	if (ans_is_lacunar && LACUNAR_MODE_IS_ON) {
		ans_nzvals = R_NilValue;
	} else {
		ans_nzvals = PROTECT(_new_Rvector1(ans_Rtype, ans_nzcount));
	}
	SEXP ans_nzoffs = PROTECT(NEW_INTEGER(ans_nzcount));

	int k1 = 0, offset = 0;
	for (int n = 0; n < nb_objects; n++) {
		SEXP leaf = leaves[n];
		if (leaf != R_NilValue) {
			SEXP nzvals, nzoffs;
			int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
			if (ans_nzvals != R_NilValue && nzvals != R_NilValue)
				_copy_Rvector_elts(nzvals, 0,
						   ans_nzvals, k1, nzcount);
			for (int k2 = 0; k2 < nzcount; k2++, k1++)
				INTEGER(ans_nzoffs)[k1] = INTEGER(nzoffs)[k2] +
							  offset;
		}
		offset += dims_along[n];
	}
	SEXP ans = zip_leaf(ans_nzvals, ans_nzoffs);
	UNPROTECT(ans_nzvals == R_NilValue ? 1 : 2);
	/* Sanity checks (should never fail). */
	if (k1 != ans_nzcount)
		error("SparseArray internal error in "
		      "concatenate_leaves():\n"
		      "    k1 != ans_nzcount");
	if (offset != sum_dims_along)
		error("SparseArray internal error in "
		      "concatenate_leaves():\n"
		      "    offset != sum_dims_along");
	return ans;
}

/* Recursive.
   'along0' will always be >= 0 and < 'ndim'.
   Returns R_NilValue or a list of length 'ans_dim[ndim - 1]'. */
static SEXP REC_abind_SVTs(SEXP *SVTs, int nb_objects,
		const int *ans_dim, int ndim, int along0, const int *dims_along,
		SEXPTYPE ans_Rtype)
{
	if (all_NULLs(SVTs, nb_objects))
		return R_NilValue;

	if (ndim == 1) {
		/* This is the case where we are binding SVT_SparseArray
		   objects along their first dimension. */
		return concatenate_leaves(SVTs, nb_objects,
					dims_along, ans_dim[along0],
					ans_Rtype);
	}

	if (along0 == ndim - 1)
		return concatenate_SVTs(SVTs, nb_objects,
					dims_along, ans_dim[along0]);

	SEXP *subSVTs_buf = SVTs + nb_objects;
	int ans_len = ans_dim[ndim - 1];
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	int is_empty = 1;
	for (int i = 0; i < ans_len; i++) {
		int ret = collect_SVTs_ith_elt(SVTs, nb_objects, i, ans_len,
					       subSVTs_buf);
		if (ret < 0) {
			UNPROTECT(1);
			error("SparseArray internal error in "
			      "REC_abind_SVTs():\n"
			      "    collect_SVTs_ith_elt() returned an error");
		}
		SEXP ans_elt = REC_abind_SVTs(subSVTs_buf, nb_objects,
					ans_dim, ndim - 1, along0, dims_along,
					ans_Rtype);
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
SEXP C_abind_SVT_SparseArray_objects(SEXP objects, SEXP SVTslotname,
		SEXP along, SEXP ans_type)
{
	if (!isVectorList(objects))  // IS_LIST() is broken
		error("'objects' must be a list of SVT_SparseArray objects");

	SEXPTYPE ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (ans_Rtype == 0)
		error("invalid requested type");

	if (!IS_INTEGER(along) || XLENGTH(along) != 1)
		error("'along' must be a single positive integer");
	int along0 = INTEGER(along)[0] - 1;

	int nb_objects = LENGTH(objects);
	if (nb_objects == 0)
		error("'objects' cannot be an empty list");

	int *dims_along = (int *) R_alloc(nb_objects, sizeof(int));
	SEXP ans_dim = PROTECT(
		check_and_combine_object_dims(objects, along0, dims_along)
	);
	int ndim = LENGTH(ans_dim);

	SEXP *SVTs_buf = prepare_SVTs_buf(objects, SVTslotname, ndim, along0);
	SEXP ans_SVT = REC_abind_SVTs(SVTs_buf, nb_objects,
				INTEGER(ans_dim), ndim, along0, dims_along,
				ans_Rtype);
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

