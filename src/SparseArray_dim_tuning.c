/****************************************************************************
 *                    Dim tuning of a SparseArray object                    *
 ****************************************************************************/
#include "SparseArray_dim_tuning.h"

#include "Rvector_utils.h"
#include "leaf_utils.h"  /* for unzip_leaf() */

#include <string.h>  /* for memset() */


#define	KEEP_DIM        0
#define	DROP_DIM       -1
#define	ADD_DIM         1

/* REC_tune_SVT() below in this file does not support a 'dim_tuner' vector
   where 1 and -1 values are neighbors. In other words, if a 'dim_tuner'
   vector contains both 1 and -1 values, then there must be at least one 0
   between them. Such a 'dim_tuner' vector is considered to be "normalized". */
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

/* set_cumallKEEP_cumallDROP() also validates the content of 'ops' ('dim_tuner'
   argument at the R level). See "Dim tuning and the 'dim_tuner' argument" in
   src/dim_tuning_utils.c in the S4Arrays package for a description of
   this argument. */
static void set_cumallKEEP_cumallDROP(int *cumallKEEP, int *cumallDROP,
		const int *ops, int nops, const int *dim, int ndim)
{
	int along1, along2, nkept, r, op;

	memset(cumallKEEP, 0, sizeof(int) * ndim);
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
			      "set_cumallKEEP_cumallDROP():\n"
			      "    number of 0 (KEEP) or -1 (DROP) values "
			      "in 'dim_tuner' is > 'length(dim(x))'");
		if (op == KEEP_DIM) {
			if (r == along1 && (r == 0 || cumallKEEP[r - 1]))
				cumallKEEP[r] = 1;
			along2++;
			nkept++;
			along1++;
			continue;
		}
		if (op != DROP_DIM)
			error("SparseArray internal error in "
			      "set_cumallKEEP_cumallDROP():\n"
			      "    'dim_tuner' can only contain 0 (KEEP), "
			      "-1 (DROP), or 1 (ADD) values");
		if (dim[along1] != 1)
			error("SparseArray internal error in "
			      "set_cumallKEEP_cumallDROP():\n"
			      "    'dim_tuner[%d]' (= -1) is "
			      "mapped to 'dim(x)[%d]' (= %d)\n"
			      "    which cannot be dropped",
			      r + 1, along1 + 1, dim[along1]);
		if (r == along1 && (r == 0 || cumallDROP[r - 1]))
			cumallDROP[r] = 1;
		along1++;
	}
	if (along1 < ndim)
		error("SparseArray internal error in "
		      "set_cumallKEEP_cumallDROP():\n"
		      "    number of 0 (KEEP) or -1 (DROP) values "
		      "in 'dim_tuner' is < 'length(dim(x))'");
	if (nkept == 0)
		error("SparseArray internal error in "
		      "set_cumallKEEP_cumallDROP():\n"
		      "    'dim_tuner' must contain at least one 0");
	return;
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
 * Go back and forth between a leaf and a 1x1x..xN SVT
 */

/* Returns a leaf representing a single nonzero value ("scalar leaf"). */
static SEXP wrap_Rvector_elt_in_scalar_leaf(SEXP in_Rvector, int k,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	SEXP ans_nzvals = PROTECT(allocVector(TYPEOF(in_Rvector), 1));
	SEXP ans_nzoffs = PROTECT(NEW_INTEGER(1));
	INTEGER(ans_nzoffs)[0] = 0;
	copy_Rvector_elt_FUN(in_Rvector, k, ans_nzvals, 0);
	SEXP ans = zip_leaf(ans_nzvals, ans_nzoffs);
	UNPROTECT(2);
	return ans;
}

static void copy_scalar_leaf_val_to_Rvector(SEXP scalar_leaf,
		SEXP out_Rvector, int k,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(scalar_leaf, &nzvals, &nzoffs);
	/* Sanity checks. */
	if (nzcount != 1 || INTEGER(nzoffs)[0] != 0)
		error("SparseArray internal error in "
		      "copy_scalar_leaf_val_to_Rvector():\n"
		      "    not a scalar leaf");
	copy_Rvector_elt_FUN(nzvals, 0, out_Rvector, k);
	return;
}


/* Transform leaf into 1x1x..xN SVT.
   'leaf' is assumed to represent a sparse vector of length N (i.e. its dense
   form would be of length N). unroll_leaf_as_SVT() turns it into an SVT that
   represents a 1x1x..xN array.
   'ans_ndim' is the number of dimensions of the result. It must be >= 2. */
static SEXP unroll_leaf_as_SVT(SEXP leaf, int N, int ans_ndim,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &nzvals, &nzoffs);
	SEXP ans = PROTECT(NEW_LIST(N));
	for (int k = 0; k < nzcount; k++) {
		int i = INTEGER(nzoffs)[k];
		SEXP ans_elt = PROTECT(
			wrap_Rvector_elt_in_scalar_leaf(nzvals, k,
							copy_Rvector_elt_FUN)
		);
		ans_elt = PROTECT(add_outermost_dims(ans_elt, ans_ndim - 2));
		SET_VECTOR_ELT(ans, i, ans_elt);
		UNPROTECT(2);
	}
	UNPROTECT(1);
	return ans;
}

/* Transform 1x1x..xN SVT into leaf.
   'SVT' is assumed to represent a 1x1x..xN array.
   More precisely: 'ndim' is assumed to be >= 2. Except maybe for its
   outermost dimension, all the dimensions in 'SVT' are assumed to be
   ineffective.
   roll_SVT_into_leaf() turns 'SVT' into a leaf that represents a
   sparse vector of length N (i.e. its dense form would be of length N). */
static SEXP roll_SVT_into_leaf(SEXP SVT, int ndim, SEXPTYPE Rtype,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	int N = LENGTH(SVT);
	int ans_nzcount = 0;
	for (int i = 0; i < N; i++) {
		SEXP subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue)
			continue;
		ans_nzcount++;
	}
	if (ans_nzcount == 0)
		error("SparseArray internal error in roll_SVT_into_leaf():\n"
		      "    ans_nzcount == 0");
	SEXP ans_nzvals, ans_nzoffs;
	SEXP ans = PROTECT(_alloc_and_unzip_leaf(Rtype, ans_nzcount,
						 &ans_nzvals, &ans_nzoffs));
	ans_nzcount = 0;
	for (int i = 0; i < N; i++) {
		SEXP subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue)
			continue;
		subSVT = drop_outermost_dims(subSVT, ndim - 2);
		/* 'subSVT' must be a scalar leaf or this will fail. */
		copy_scalar_leaf_val_to_Rvector(subSVT, ans_nzvals, ans_nzcount,
						copy_Rvector_elt_FUN);
		INTEGER(ans_nzoffs)[ans_nzcount] = i;
		ans_nzcount++;
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_tune_SVT_dims()
 */

/* Assumes that 'dim_tuner' is normalized.
   Recursive. */
static SEXP REC_tune_SVT(SEXP SVT, const int *dim, int ndim,
		const int *ops, int nops,
		const int *cumallKEEP, const int *cumallDROP,
		SEXPTYPE Rtype, CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	int op, ans_len, i;
	SEXP ans_elt, ans, subSVT;

	if (SVT == R_NilValue || (nops == ndim && cumallKEEP[ndim - 1]))
		return SVT;

	op = ops[nops - 1];
	if (op == ADD_DIM) {
		/* Add ineffective dimension (as outermost dimension). */
		ans_elt = PROTECT(
			REC_tune_SVT(SVT, dim, ndim,
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
			   on its left. 'SVT' is a leaf (i.e. 1D SVT). */
			return unroll_leaf_as_SVT(SVT, dim[0], nops,
						  copy_Rvector_elt_FUN);
		}
		if (nops == ndim && cumallDROP[ndim - 2]) {
			/* 'ops[nops - 1]' is KEEP_DIM, with only DROP_DIM ops
			   on its left. Returns a leaf (i.e. 1D SVT). */
			return roll_SVT_into_leaf(SVT, ndim, Rtype,
						  copy_Rvector_elt_FUN);
		}
		ans_len = dim[ndim - 1];
		ans = PROTECT(NEW_LIST(ans_len));
		for (i = 0; i < ans_len; i++) {
			subSVT = VECTOR_ELT(SVT, i);
			ans_elt = PROTECT(
				REC_tune_SVT(subSVT, dim, ndim - 1,
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
	   so 'SVT' cannot be a leaf. */
	return REC_tune_SVT(VECTOR_ELT(SVT, 0), dim, ndim - 1,
			    ops, nops - 1,
			    cumallKEEP, cumallDROP,
			    Rtype, copy_Rvector_elt_FUN);
}

/* --- .Call ENTRY POINT ---
   See "Dim tuning and the 'dim_tuner' argument" in src/dim_tuning_utils.c
   in the S4Arrays package for a description of the 'dim_tuner' argument. */
SEXP C_tune_SVT_dims(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP dim_tuner)
{
	SEXPTYPE Rtype;
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;
	int ndim, nops, *cumallKEEP, *cumallDROP;
	const int *dim, *ops;

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

	dim = INTEGER(x_dim);
	cumallKEEP = (int *) R_alloc(ndim, sizeof(int));
	cumallDROP = (int *) R_alloc(ndim, sizeof(int));
	set_cumallKEEP_cumallDROP(cumallKEEP, cumallDROP,
				  ops, nops, dim, ndim);

	/* Compute tuned 'SVT'. */
	return REC_tune_SVT(x_SVT, dim, ndim, ops, nops,
			    cumallKEEP, cumallDROP,
			    Rtype, copy_Rvector_elt_FUN);
}

