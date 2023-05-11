/****************************************************************************
 *                    Dim tuning of a SparseArray object                    *
 ****************************************************************************/
#include "SparseArray_dim_tuning.h"

#include "Rvector_utils.h"
#include "leaf_vector_utils.h"  /* for _split_leaf_vector() */

#include <string.h>  /* for memset() */



/****************************************************************************
 * Dim tuning and the 'dim_tuner' argument
 *
 * Dim tuning
 * ----------
 * Dim tuning is the act of adding and/or dropping ineffective dimensions,
 * (i.e. dimensions that have an extent of 1) to/from an array-like object.
 * Note that dim tuning doesn't change the length (which is prod(dim(.)))
 * or alter the content of the object, and is always reversible (except when
 * it drops ineffective dimensions with names on them).
 *
 * The 'dim_tuner' argument
 * ------------------------
 * The exact action to perform on the dimensions of the object is encoded
 * in the 'dim_tuner' argument. This is an integer vector where each value
 * represents one of three possible operations:
 *   o  0: Keep the dimension.
 *   o -1: Drop the dimension. This operation is allowed only if the
 *         dimension to drop is ineffective (i.e. has an extent of 1).
 *   o  1: Add ineffective dimension.
 * Note that the 'dim_tuner' vector can contain any number of 1's, but the
 * number of non-positive values (i.e. number of 0 and -1 values together)
 * must match the number of dimensions of the array-like object to tune.
 *
 * Additionally, 'dim_tuner' must contain at least one 0. In other words,
 * the tuning must retain at least one of the original dimensions of the
 * object.
 *
 * Normalized 'dim_tuner' vector
 * -----------------------------
 * Note that REC_tune_SVT() does not support a 'dim_tuner' vector where 1
 * and -1 values are neighbors. In other words, if a 'dim_tuner' vector
 * contains both 1 and -1 values, then there must be at least one 0 between
 * them. Such a 'dim_tuner' vector is considered to be "normalized".
 *
 * Reverse tuning
 * --------------
 * To revert a dim tuning, simply tune again with '- dim_tuner' (i.e. minus
 * 'dim_tuner'). More precisely, if 'dim' is a vector of dimensions and
 * if 'dim_tuner' represents dim tuning compatible with 'dim':
 *
 *     tuned_dim <- SparseArray:::.tune_dims(dim, dim_tuner)
 *     dim2 <- SparseArray:::.tune_dims(tuned_dim, - dim_tuner)
 */

#define	KEEP_DIM        0
#define	DROP_DIM       -1
#define	ADD_DIM         1

static int dim_tuner_is_normalized(const int *ops, int nops)
{
	int prev_op, r, op;

	prev_op = ops[0];  /* 'nops' is guaranteed to be >= 1 */
	for (r = 1; r < nops; r++) {
		op = ops[r];  /* -1 <= op <= 1 */
		if (prev_op * op < 0)
			return 0;
		prev_op = op;
	}
	return 1;
}

/* Return the "new" number of dimensions i.e. the number of dims that we
   will get after tuning the current vector of dimensions. Note that this
   is simply the number of non-negative values in 'ops' (i.e. number of
   0 and 1 values together). */
static int validate_dim_tuner(const int *ops, int nops,
		const int *dims, int ndim,
		int *cumallKEEP, int *cumallDROP)
{
	int along1, along2, nkept, r, op;

	if (cumallKEEP != NULL)
		memset(cumallKEEP, 0, sizeof(int) * ndim);
	if (cumallDROP != NULL)
		memset(cumallDROP, 0, sizeof(int) * ndim);
	along1 = along2 = nkept = 0;
	for (r = 0; r < nops; r++) {
		op = ops[r];  /* ADD_DIM, KEEP_DIM, or DROP_DIM */
		if (op == ADD_DIM) {
			along2++;
			continue;
		}
		if (along1 >= ndim)
			error("SparseArray internal error in "
			      "validate_dim_tuner():\n"
			      "    number of 0 (KEEP) or -1 (DROP) values "
			      "in 'dim_tuner' is > 'length(dim(x))'");
		if (op == KEEP_DIM) {
			if (cumallKEEP != NULL &&
			    r == along1 && (r == 0 || cumallKEEP[r - 1]))
				cumallKEEP[r] = 1;
			along2++;
			nkept++;
			along1++;
			continue;
		}
		if (op != DROP_DIM)
			error("SparseArray internal error in "
			      "validate_dim_tuner():\n"
			      "    'dim_tuner' can only contain 0 (KEEP), "
			      "-1 (DROP), or 1 (ADD) values");
		if (dims[along1] != 1)
			error("SparseArray internal error in "
			      "validate_dim_tuner():\n"
			      "    'dim_tuner[%d]' (= -1) is "
			      "mapped to 'dim(x)[%d]' (= %d)\n"
			      "    which cannot be dropped",
			      r + 1, along1 + 1, dims[along1]);
		if (cumallDROP != NULL &&
		    r == along1 && (r == 0 || cumallDROP[r - 1]))
			cumallDROP[r] = 1;
		along1++;
	}
	if (along1 < ndim)
		error("SparseArray internal error in "
		      "validate_dim_tuner():\n"
		      "    number of 0 (KEEP) or -1 (DROP) values "
		      "in 'dim_tuner' is < 'length(dim(x))'");
	if (nkept == 0)
		error("SparseArray internal error in "
		      "validate_dim_tuner():\n"
		      "    'dim_tuner' must contain at least one 0");
	return along2;
}


/****************************************************************************
 * C_tune_dims() and C_tune_dimnames()
 */

/* 'dim_names' must be a character vector (or R_NilValue) that contains the
   names on the vector of dimensions. This is NOT the same as the dimnames! */
static SEXP tune_dims(const int *dims, SEXP dim_names,
		const int *ops, int nops, int ans_len)
{
	SEXP ans, ans_names;
	int along1, along2, r, op;

	ans = PROTECT(NEW_INTEGER(ans_len));
	if (dim_names != R_NilValue)
		ans_names = PROTECT(NEW_CHARACTER(ans_len));
	along1 = along2 = 0;
	for (r = 0; r < nops; r++) {
		op = ops[r];  /* ADD_DIM, KEEP_DIM, or DROP_DIM */
		if (op == ADD_DIM) {
			INTEGER(ans)[along2] = 1;
			along2++;
			continue;
		}
		if (op == KEEP_DIM) {
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
		const int *ops, int nops, int ans_len)
{
	int along1, along2, r, op;
	SEXP ans;

	ans = PROTECT(NEW_LIST(ans_len));
	along1 = along2 = 0;
	for (r = 0; r < nops; r++) {
		op = ops[r];  /* ADD_DIM, KEEP_DIM, or DROP_DIM */
		if (op == ADD_DIM) {
			along2++;
			continue;
		}
		if (op == KEEP_DIM) {
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
		const int *ops, int nops)
{
	int ndim, along1, along2, any_retained, r, op;

	if (dimnames == R_NilValue)
		return 0;
	ndim = LENGTH(dimnames);
	along1 = along2 = any_retained = 0;
	for (r = 0; r < nops; r++) {
		op = ops[r];  /* ADD_DIM, KEEP_DIM, or DROP_DIM */
		if (op == ADD_DIM) {
			along2++;
			continue;
		}
		if (along1 >= ndim)
			error("SparseArray internal error in "
			      "compute_tuned_dimnames_length():\n"
			      "    number of 0 (KEEP) or -1 (DROP) values "
			      "in 'dim_tuner' is > 'length(dim(x))'");
		if (op == KEEP_DIM) {
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
   a 'dim_tuner' vector that is not normalized. */
SEXP C_tune_dims(SEXP dim, SEXP dim_tuner)
{
	int ndim, nops, ans_len;
	const int *dims, *ops;

	ndim = LENGTH(dim);
	dims = INTEGER(dim);
	nops = LENGTH(dim_tuner);
	ops = INTEGER(dim_tuner);
	ans_len = validate_dim_tuner(ops, nops, dims, ndim, NULL, NULL);
	return tune_dims(dims, GET_NAMES(dim), ops, nops, ans_len);
}

/* --- .Call ENTRY POINT ---
   Unlike C_tune_SVT_dims() below in this file, C_tune_dimnames() accepts
   a 'dim_tuner' vector that is not normalized. */
SEXP C_tune_dimnames(SEXP dimnames, SEXP dim_tuner)
{
	int nops, ans_len;
	const int *ops;

	nops = LENGTH(dim_tuner);
	ops = INTEGER(dim_tuner);
	ans_len = compute_tuned_dimnames_length(dimnames, ops, nops);
	if (ans_len == 0)
		return R_NilValue;
	return tune_dimnames(dimnames, ops, nops, ans_len);
}


/****************************************************************************
 * Add/drop ineffective dimensions as outermost dimensions.
 */

/* Add ineffective dimensions as outermost dimensions.
   Caller must PROTECT() the result. */
static SEXP add_outermost_dims(SEXP SVT, int ndim_to_add)
{
	SEXP ans, tmp;
	int along;

	if (ndim_to_add <= 0)
		return SVT;
	ans = PROTECT(NEW_LIST(1));
	SET_VECTOR_ELT(ans, 0, SVT);
	for (along = 1; along < ndim_to_add; along++) {
		tmp = PROTECT(NEW_LIST(1));
		SET_VECTOR_ELT(tmp, 0, VECTOR_ELT(ans, 0));
		SET_VECTOR_ELT(ans, 0, tmp);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

/* Drop outermost ineffective dimensions.
   Caller does NOT need to PROTECT() the result. */
static SEXP drop_outermost_dims(SEXP SVT, int ndim_to_drop)
{
	int along;

	for (along = 0; along < ndim_to_drop; along++) {
		/* Sanity check. */
		if (SVT == R_NilValue || LENGTH(SVT) != 1)
			error("SparseArray internal error in "
			      "drop_outermost_dims():\n"
			      "    'SVT' not as expected");
		SVT = VECTOR_ELT(SVT, 0);
	}
	return SVT;
}


/****************************************************************************
 * Go back and forth between a "leaf vector" and a 1x1x..xN SVT
 */

/* Return a "leaf vector" of length 1. */
static SEXP wrap_Rvector_elt_in_lv1(SEXP in_Rvector, int k,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	SEXP ans_offs, ans_vals, ans;

	ans_offs = PROTECT(NEW_INTEGER(1));
	ans_vals = PROTECT(allocVector(TYPEOF(in_Rvector), 1));
	INTEGER(ans_offs)[0] = 0;
	copy_Rvector_elt_FUN(in_Rvector, k, ans_vals, 0);
	ans = _new_leaf_vector(ans_offs, ans_vals);
	UNPROTECT(2);
	return ans;
}

/* 'lv' must be a "leaf vector" of length 1. */
static void copy_lv1_val_to_Rvector(SEXP lv, SEXP out_Rvector, int k,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	int lv_len;
	SEXP lv_offs, lv_vals;

	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	/* Sanity checks. */
	if (lv_len != 1 || INTEGER(lv_offs)[0] != 0)
		error("SparseArray internal error in "
		      "copy_lv1_val_to_Rvector():\n"
		      "    leaf vector not as expected");
	copy_Rvector_elt_FUN(lv_vals, 0, out_Rvector, k);
	return;
}


/* From "leaf vector" to 1x1x..xN SVT.
   'lv' is assumed to be a "leaf vector" that represents a sparse vector
   of length N. unroll_lv_as_SVT() turns it into an SVT that represents
   a 1x1x..xN array. 'ans_ndim' is the number of dimensions of the result.
   It must be >= 2. */
static SEXP unroll_lv_as_SVT(SEXP lv, int N, int ans_ndim,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	int lv_len, k, i;
	SEXP lv_offs, lv_vals, ans, ans_elt;

	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	ans = PROTECT(NEW_LIST(N));
	for (k = 0; k < lv_len; k++) {
		i = INTEGER(lv_offs)[k];
		ans_elt = PROTECT(
			wrap_Rvector_elt_in_lv1(lv_vals, k,
						copy_Rvector_elt_FUN)
		);
		ans_elt = PROTECT(add_outermost_dims(ans_elt, ans_ndim - 2));
		SET_VECTOR_ELT(ans, i, ans_elt);
		UNPROTECT(2);
	}
	UNPROTECT(1);
	return ans;
}

/* From 1x1x..xN SVT to "leaf vector".
   'SVT' is assumed to represent a 1x1x..xN array.
   More precisely: 'ndim' is assumed to be >= 2. Except maybe for its
   outermost dimension, all the dimensions in 'SVT' are assumed to be
   ineffective.
   roll_SVT_into_lv() turns 'SVT' into a "leaf vector" that represents a
   sparse vector of length N. */
static SEXP roll_SVT_into_lv(SEXP SVT, int ndim, SEXPTYPE Rtype,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	int N, ans_len, i;
	SEXP subSVT, ans_offs, ans_vals, ans;

	N = LENGTH(SVT);
	ans_len = 0;
	for (i = 0; i < N; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue)
			continue;
		ans_len++;
	}
	if (ans_len == 0)
		error("SparseArray internal error in "
		      "roll_SVT_into_lv():\n"
		      "    ans_len == 0");
	ans_offs = PROTECT(NEW_INTEGER(ans_len));
	ans_vals = PROTECT(allocVector(Rtype, ans_len));
	ans_len = 0;
	for (i = 0; i < N; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue)
			continue;
		subSVT = drop_outermost_dims(subSVT, ndim - 2);
		/* 'subSVT' is a "leaf vector" of length 1. */
		copy_lv1_val_to_Rvector(subSVT, ans_vals, ans_len,
					copy_Rvector_elt_FUN);
		INTEGER(ans_offs)[ans_len] = i;
		ans_len++;
	}
	ans = _new_leaf_vector(ans_offs, ans_vals);
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * C_tune_SVT_dims()
 */

/* Assumes that 'dim_tuner' is normalized.
   Recursive. */
static SEXP REC_tune_SVT(SEXP SVT, const int *dims, int ndim,
		const int *ops, int nops,
		const int *cumallKEEP, const int *cumallDROP,
		SEXPTYPE Rtype, CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	int op, ans_len, i;
	SEXP ans_elt, ans, subSVT;

	if (SVT == R_NilValue || nops == ndim && cumallKEEP[ndim - 1])
		return SVT;

	op = ops[nops - 1];
	if (op == ADD_DIM) {
		/* Add ineffective dimension (as outermost dimension). */
		ans_elt = PROTECT(
			REC_tune_SVT(SVT, dims, ndim,
				     ops, nops - 1,
				     cumallKEEP, cumallDROP,
				     Rtype, copy_Rvector_elt_FUN)
		);
		ans = PROTECT(add_outermost_dims(ans_elt, 1));
		UNPROTECT(2);
		return ans;
	}
	if (op == KEEP_DIM) {
		if (ndim == 1) {
			/* 'ops[nops - 1]' is KEEP_DIM, with only ADD_DIM ops
			   on its left. 'SVT' is a "leaf vector". */
			return unroll_lv_as_SVT(SVT, dims[0], nops,
						copy_Rvector_elt_FUN);
		}
		if (nops == ndim && cumallDROP[ndim - 2]) {
			/* 'ops[nops - 1]' is KEEP_DIM, with only DROP_DIM ops
			   on its left. Return a "leaf vector". */
			return roll_SVT_into_lv(SVT, ndim - 2, Rtype,
						copy_Rvector_elt_FUN);
		}
		ans_len = dims[ndim - 1];
		ans = PROTECT(NEW_LIST(ans_len));
		for (i = 0; i < ans_len; i++) {
			subSVT = VECTOR_ELT(SVT, i);
			ans_elt = PROTECT(
				REC_tune_SVT(subSVT, dims, ndim - 1,
					     ops, nops - 1,
					     cumallKEEP, cumallDROP,
					     Rtype, copy_Rvector_elt_FUN)
			);
			SET_VECTOR_ELT(ans, i, ans_elt);
			UNPROTECT(1);
		}
		UNPROTECT(1);
		return ans;
	}
	/* Drop ineffective dimension.
	   Because the 'ops' vector is normalized, it's guaranteed to contain
	   at least one KEEP_DIM op on the left of the DROP_DIM op found at
	   position 'nops - 1'.
	   Furthermore, the closest KEEP_DIM op (i.e. highest KEEP_DIM's
	   position that is < 'nops - 1') is guaranteed to be separated from
	   the DROP_DIM op at position 'nops - 1' by nothing but other
	   DROP_DIM ops.
	   In particular, this means that 'ndim' is guaranteed to be >= 2
	   so 'SVT' cannot be a "leaf vector". */
	return REC_tune_SVT(VECTOR_ELT(SVT, 0), dims, ndim - 1,
			    ops, nops - 1,
			    cumallKEEP, cumallDROP,
			    Rtype, copy_Rvector_elt_FUN);
}

/* --- .Call ENTRY POINT ---
   See "Dim tuning and the 'dim_tuner' argument" at the top of this
   file for a description of the 'dim_tuner' argument. */
SEXP C_tune_SVT_dims(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP dim_tuner)
{
	SEXPTYPE Rtype;
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;
	int ndim, nops, *cumallKEEP, *cumallDROP;
	const int *dims, *ops;

	Rtype = _get_Rtype_from_Rstring(x_type);
	copy_Rvector_elt_FUN = _select_copy_Rvector_elt_FUN(Rtype);
	if (copy_Rvector_elt_FUN == NULL)
		error("SparseArray internal error in "
		      "C_tune_SVT_dims():\n"
		      "    SVT_SparseArray object has invalid type");

	/* Make sure that: 1 <= ndim <= nops. */
	ndim = LENGTH(x_dim);
	if (ndim == 0)
		error("SparseArray internal error in "
		      "C_tune_SVT_dims():\n"
		      "    'dim(x)' cannot be empty");
	nops = LENGTH(dim_tuner);
	if (nops < ndim)
		error("SparseArray internal error in "
		      "C_tune_SVT_dims():\n"
		      "    length(dim_tuner) < length(dim(x))");

	ops = INTEGER(dim_tuner);
	/* REC_tune_SVT() assumes that the 'ops' vector is normalized.
	   Note that we have no use case for an 'ops' vector that is not
	   normalized at the moment. */
	if (!dim_tuner_is_normalized(ops, nops))
		error("SparseArray internal error in "
		      "C_tune_SVT_dims():\n"
		      "    'dim_tuner' is not normalized");

	dims = INTEGER(x_dim);
	cumallKEEP = (int *) R_alloc(ndim, sizeof(int));
	cumallDROP = (int *) R_alloc(ndim, sizeof(int));
	validate_dim_tuner(ops, nops, dims, ndim, cumallKEEP, cumallDROP);

	/* Compute tuned 'SVT'. */
	return REC_tune_SVT(x_SVT, dims, ndim, ops, nops,
			    cumallKEEP, cumallDROP,
			    Rtype, copy_Rvector_elt_FUN);
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
	ans = PROTECT(roll_SVT_into_lv(ans, 2, Rtype, copy_Rvector_elt_FUN));
	UNPROTECT(2);
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

