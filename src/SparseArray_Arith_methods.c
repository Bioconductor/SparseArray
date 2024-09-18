/****************************************************************************
 *                'Arith' operations on SparseArray objects                 *
 ****************************************************************************/
#include "SparseArray_Arith_methods.h"

#include "argcheck_utils.h"
#include "SparseVec.h"
#include "SparseVec_Arith.h"
#include "leaf_utils.h"
#include "SVT_SparseArray_class.h"  /* for _coerce_SVT() */


/****************************************************************************
 * 'Arith' operations on the tree leaves
 */

/* Assumes that 'ans_Rtype' is equal or bigger than the type of the nonzero
   values in 'leaf'. Performs **in-place** replacement if 'ans_Rtype' is 0!
   Note that this could also have been achieved by calling:

     Arith_leaf1_scalar(MULT_OPCODE, leaf, Rtype, 0, -1, ...)

   but unary_minus_leaf() takes a lot of shortcuts so is A LOT more
   efficient. */
static SEXP unary_minus_leaf(SEXP leaf, SEXPTYPE Rtype, SEXPTYPE ans_Rtype)
{
	SEXP leaf_nzvals, ans_nzvals, nzoffs;
	int nzcount = unzip_leaf(leaf, &leaf_nzvals, &nzoffs);
	if (leaf_nzvals == R_NilValue) {  /* input leaf is lacunar */
		ans_nzvals = PROTECT(
			allocVector(ans_Rtype == 0 ? Rtype : ans_Rtype, nzcount)
		);
		_set_Rvector_elts_to_minus_one(ans_nzvals);
		if (ans_Rtype == 0) {  /* in-place replacement! */
			replace_leaf_nzvals(leaf, ans_nzvals);
			UNPROTECT(1);
			return leaf;
		}
		SEXP ans = zip_leaf(ans_nzvals, nzoffs, 0);
		UNPROTECT(1);
		return ans;
	}
	/* input leaf is standard */
	if (ans_Rtype == 0) {  /* in-place replacement! */
		ans_nzvals = leaf_nzvals;
	} else {
		ans_nzvals = PROTECT(allocVector(ans_Rtype, nzcount));
	}
	_unary_minus_Rvector(leaf_nzvals, ans_nzvals);
	int go_lacunar = LACUNAR_MODE_IS_ON &&
			 _all_Rvector_elts_equal_one(ans_nzvals);
	if (ans_Rtype == 0) {
		if (go_lacunar)
			replace_leaf_nzvals(leaf, R_NilValue);
		return leaf;
	}
	if (go_lacunar)
		ans_nzvals = R_NilValue;
	SEXP ans = zip_leaf(ans_nzvals, nzoffs, 0);
	UNPROTECT(1);
	return ans;
}

static SEXP Arith_leaf1_scalar(int opcode,
		SEXP leaf1, SEXPTYPE Rtype1, SEXP scalar,
		SparseVec *buf_sv, int *ovflow)
{
	const SparseVec sv1 = leaf2SV(leaf1, Rtype1,
				      buf_sv->len, buf_sv->na_background);
	_Arith_sv1_scalar(opcode, &sv1, scalar, buf_sv, ovflow);
	if (buf_sv->nzcount == PROPAGATE_NZOFFS)
		return _make_leaf_with_single_shared_nzval(
					      buf_sv->Rtype, buf_sv->nzvals,
					      get_leaf_nzoffs(leaf1));
	return _make_leaf_from_two_arrays(buf_sv->Rtype, buf_sv->nzvals,
					  buf_sv->nzoffs, buf_sv->nzcount);
}

static SEXP Arith_scalar_leaf2(int opcode,
		SEXP scalar, SEXP leaf2, SEXPTYPE Rtype2,
		SparseVec *buf_sv, int *ovflow)
{
	const SparseVec sv2 = leaf2SV(leaf2, Rtype2,
				      buf_sv->len, buf_sv->na_background);
	_Arith_scalar_sv2(opcode, scalar, &sv2, buf_sv, ovflow);
	if (buf_sv->nzcount == PROPAGATE_NZOFFS)
		return _make_leaf_with_single_shared_nzval(
					      buf_sv->Rtype, buf_sv->nzvals,
					      get_leaf_nzoffs(leaf2));
	return _make_leaf_from_two_arrays(buf_sv->Rtype, buf_sv->nzvals,
					  buf_sv->nzoffs, buf_sv->nzcount);
}

static SEXP Arith_NULL_leaf2(int opcode,
		SEXPTYPE Rtype1, int na_background1,
		SEXP leaf2, SEXPTYPE Rtype2, int na_background2,
		SparseVec *buf_sv)
{
	if (leaf2 == R_NilValue)
		error("SparseArray internal error in "
		      "Arith_NULL_leaf2():\n"
		      "    'leaf2' cannot be NULL");
	const SparseVec sv2 = leaf2SV(leaf2, Rtype2,
				      buf_sv->len, na_background2);
	if (na_background1) {
		_Arith_na_sv2(opcode, Rtype1, &sv2, buf_sv);
	} else {
		if (opcode == SUB_OPCODE) {
			/* Should be A LOT more efficient than the call
			   to _Arith_0_sv(opcode, ...) below. */
			return unary_minus_leaf(leaf2, Rtype2, buf_sv->Rtype);
		}
		_Arith_zero_sv2(opcode, Rtype1, &sv2, buf_sv);
	}
	if (buf_sv->nzcount == PROPAGATE_NZOFFS)
		return _make_leaf_with_single_shared_nzval(
					      buf_sv->Rtype, buf_sv->nzvals,
					      get_leaf_nzoffs(leaf2));
	return _make_leaf_from_two_arrays(buf_sv->Rtype, buf_sv->nzvals,
					  buf_sv->nzoffs, buf_sv->nzcount);
}

static SEXP Arith_leaf1_NULL(int opcode,
		SEXP leaf1, SEXPTYPE Rtype1, int na_background1,
		SEXPTYPE Rtype2, int na_background2,
		SparseVec *buf_sv)
{
	if (leaf1 == R_NilValue)
		error("SparseArray internal error in "
		      "Arith_leaf1_NULL():\n"
		      "    'leaf1' cannot be NULL");
	const SparseVec sv1 = leaf2SV(leaf1, Rtype1,
				      buf_sv->len, na_background1);
	if (na_background2) {
		_Arith_sv1_na(opcode, &sv1, Rtype2, buf_sv);
	} else {
		_Arith_sv1_zero(opcode, &sv1, Rtype2, buf_sv);
	}
	if (buf_sv->nzcount == PROPAGATE_NZOFFS)
		return _make_leaf_with_single_shared_nzval(
					      buf_sv->Rtype, buf_sv->nzvals,
					      get_leaf_nzoffs(leaf1));
	return _make_leaf_from_two_arrays(buf_sv->Rtype, buf_sv->nzvals,
					  buf_sv->nzoffs, buf_sv->nzcount);
}

static SEXP Arith_leaf1_leaf2(int opcode,
		SEXP leaf1, SEXPTYPE Rtype1, int na_background1,
		SEXP leaf2, SEXPTYPE Rtype2, int na_background2,
		SparseVec *buf_sv, int *ovflow)
{
	if (leaf1 == R_NilValue)
		return Arith_NULL_leaf2(opcode,
					Rtype1, na_background1,
					leaf2, Rtype2, na_background2,
					buf_sv);
	if (leaf2 == R_NilValue)
		return Arith_leaf1_NULL(opcode,
					leaf1, Rtype1, na_background1,
					Rtype2, na_background2,
					buf_sv);
	const SparseVec sv1 = leaf2SV(leaf1, Rtype1,
				      buf_sv->len, na_background1);
	const SparseVec sv2 = leaf2SV(leaf2, Rtype2,
				      buf_sv->len, na_background2);
	_Arith_sv1_sv2(opcode, &sv1, &sv2, buf_sv, ovflow);
	return _make_leaf_from_two_arrays(buf_sv->Rtype, buf_sv->nzvals,
					  buf_sv->nzoffs, buf_sv->nzcount);
}


/****************************************************************************
 * Recursive tree traversals
 */

static void REC_unary_minus_SVT(SEXP SVT, SEXPTYPE Rtype,
				const int *dim, int ndim)
{
	if (SVT == R_NilValue)
		return;
	if (ndim == 1) {
		/* 'SVT' is a leaf (i.e. 1D SVT). */
		unary_minus_leaf(SVT, Rtype, 0);  /* inplace unary minus ! */
		return;
	}
	int SVT_len = dim[ndim - 1];
	for (int i = 0; i < SVT_len; i++)
		REC_unary_minus_SVT(VECTOR_ELT(SVT, i), Rtype, dim, ndim - 1);
	return;
}

static SEXP REC_Arith_SVT1_scalar(int opcode,
		SEXP SVT1, SEXPTYPE Rtype1, SEXP scalar,
		const int *dim, int ndim,
		SparseVec *buf_sv, int *ovflow)
{
	if (SVT1 == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* 'SVT1' is a leaf (i.e. 1D SVT). */
		return Arith_leaf1_scalar(opcode,
					  SVT1, Rtype1, scalar,
					  buf_sv, ovflow);
	}

	/* 'SVT1' is a list. */
	int ans_len = dim[ndim - 1];
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	int is_empty = 1;
	for (int i = 0; i < ans_len; i++) {
		SEXP subSVT1 = VECTOR_ELT(SVT1, i);
		SEXP ans_elt = REC_Arith_SVT1_scalar(opcode,
					subSVT1, Rtype1, scalar,
					dim, ndim - 1,
					buf_sv, ovflow);
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

static SEXP REC_Arith_scalar_SVT2(int opcode,
		SEXP scalar, SEXP SVT2, SEXPTYPE Rtype2,
		const int *dim, int ndim,
		SparseVec *buf_sv, int *ovflow)
{
	if (SVT2 == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* 'SVT2' is a leaf (i.e. 1D SVT). */
		return Arith_scalar_leaf2(opcode,
					  scalar, SVT2, Rtype2,
					  buf_sv, ovflow);
	}

	/* 'SVT2' is a list. */
	int ans_len = dim[ndim - 1];
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	int is_empty = 1;
	for (int i = 0; i < ans_len; i++) {
		SEXP subSVT2 = VECTOR_ELT(SVT2, i);
		SEXP ans_elt = REC_Arith_scalar_SVT2(opcode,
					scalar, subSVT2, Rtype2,
					dim, ndim - 1,
					buf_sv, ovflow);
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

static SEXP REC_Arith_SVT1_SVT2(int opcode,
		SEXP SVT1, SEXPTYPE Rtype1, int na_background1,
		SEXP SVT2, SEXPTYPE Rtype2, int na_background2,
		const int *dim, int ndim,
		SparseVec *buf_sv, int *ovflow)
{
	if (SVT1 == R_NilValue && SVT2 == R_NilValue)
		return R_NilValue;

	if (!(na_background1 || na_background2)) {
		if (SVT1 == R_NilValue) {
			if (opcode == ADD_OPCODE)
				return _coerce_SVT(SVT2, dim, ndim, Rtype2,
						   buf_sv->Rtype,
						   buf_sv->nzoffs);
		} else if (SVT2 == R_NilValue) {
			if (opcode == ADD_OPCODE || opcode == SUB_OPCODE)
				return _coerce_SVT(SVT1, dim, ndim, Rtype1,
						   buf_sv->Rtype,
						   buf_sv->nzoffs);
		}
	}

	if (ndim == 1) {
		/* 'SVT1' and 'SVT2' are leaves (i.e. 1D SVTs).
		   Additionally, if !(na_background1 || na_background2):
		   - 'SVT1' can be NULL if 'opcode' is SUB_OPCODE;
		   - either 'SVT1' or 'SVT2' (but not both) can be NULL
		     if 'opcode' is MULT_OPCODE. */
		return Arith_leaf1_leaf2(opcode,
					 SVT1, Rtype1, na_background1,
					 SVT2, Rtype2, na_background2,
					 buf_sv, ovflow);
	}

	/* Each of 'SVT1' and 'SVT2' is either a list or NULL, but they
	   cannot both be NULL. */
	int ans_len = dim[ndim - 1];
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	SEXP subSVT1 = R_NilValue;
	SEXP subSVT2 = R_NilValue;
	int is_empty = 1;
	for (int i = 0; i < ans_len; i++) {
		if (SVT1 != R_NilValue)
			subSVT1 = VECTOR_ELT(SVT1, i);
		if (SVT2 != R_NilValue)
			subSVT2 = VECTOR_ELT(SVT2, i);
		SEXP ans_elt = REC_Arith_SVT1_SVT2(opcode,
					subSVT1, Rtype1, na_background1,
					subSVT2, Rtype2, na_background2,
					dim, ndim - 1,
					buf_sv, ovflow);
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


/****************************************************************************
 * .Call ENTRY POINTS
 */

/* --- .Call ENTRY POINT --- */
SEXP C_unary_minus_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT)
{
	SEXPTYPE x_Rtype = _get_and_check_Rtype_from_Rstring(x_type,
					"C_unary_minus_SVT", "x_type");
	SEXP ans = PROTECT(duplicate(x_SVT));
	REC_unary_minus_SVT(ans, x_Rtype, INTEGER(x_dim), LENGTH(x_dim));
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_Arith_SVT1_v2(
		SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP x_na_background,
		SEXP v2, SEXP op, SEXP ans_type)
{
	SEXPTYPE x_Rtype = _get_and_check_Rtype_from_Rstring(x_type,
					"C_Arith_SVT1_v2", "x_type");

	int x_has_NAbg = _get_and_check_na_background(x_na_background,
					"C_Arith_SVT1_v2", "x_na_background");

	SEXPTYPE ans_Rtype = _get_and_check_Rtype_from_Rstring(ans_type,
					"C_Arith_SVT1_v2", "ans_type");

	int opcode = _get_Arith_opcode(op);
	if (!x_has_NAbg && opcode != MULT_OPCODE &&
			   opcode != DIV_OPCODE &&
			   opcode != POW_OPCODE &&
			   opcode != MOD_OPCODE &&
			   opcode != IDIV_OPCODE)
	{
		error("\"%s\" is not supported between an SVT_SparseArray "
		      "object and a numeric vector", CHAR(STRING_ELT(op, 0)));
	}

	int dim0 = INTEGER(x_dim)[0];
	SparseVec buf_sv = alloc_SparseVec(ans_Rtype, dim0, x_has_NAbg);

	int ovflow = 0;
	SEXP ans = REC_Arith_SVT1_scalar(opcode,
				x_SVT, x_Rtype, v2,
				INTEGER(x_dim), LENGTH(x_dim),
				&buf_sv, &ovflow);
	if (ovflow) {
		PROTECT(ans);
		warning("NAs produced by integer overflow");
		UNPROTECT(1);
	}
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_Arith_v1_SVT2(SEXP v1,
		SEXP y_dim, SEXP y_type, SEXP y_SVT, SEXP y_na_background,
		SEXP op, SEXP ans_type)
{
	SEXPTYPE y_Rtype = _get_and_check_Rtype_from_Rstring(y_type,
					"C_Arith_v1_SVT2", "y_type");

	int y_has_NAbg = _get_and_check_na_background(y_na_background,
					"C_Arith_v1_SVT2", "y_na_background");

	SEXPTYPE ans_Rtype = _get_and_check_Rtype_from_Rstring(ans_type,
					"C_Arith_v1_SVT2", "ans_type");

	int opcode = _get_Arith_opcode(op);
	if (!y_has_NAbg && opcode != MULT_OPCODE &&
			   opcode != DIV_OPCODE &&
			   opcode != POW_OPCODE &&
			   opcode != MOD_OPCODE &&
			   opcode != IDIV_OPCODE)
	{
		error("\"%s\" is not supported between a numeric vector and "
		      "an SVT_SparseArray object", CHAR(STRING_ELT(op, 0)));
	}

	int dim0 = INTEGER(y_dim)[0];
	SparseVec buf_sv = alloc_SparseVec(ans_Rtype, dim0, y_has_NAbg);

	int ovflow = 0;
	SEXP ans = REC_Arith_scalar_SVT2(opcode,
				v1, y_SVT, y_Rtype,
				INTEGER(y_dim), LENGTH(y_dim),
				&buf_sv, &ovflow);
	if (ovflow) {
		PROTECT(ans);
		warning("NAs produced by integer overflow");
		UNPROTECT(1);
	}
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_Arith_SVT1_SVT2(
		SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP x_na_background,
		SEXP y_dim, SEXP y_type, SEXP y_SVT, SEXP y_na_background,
		SEXP op, SEXP ans_type)
{
	_check_array_conformability(x_dim, y_dim);
	SEXPTYPE x_Rtype = _get_and_check_Rtype_from_Rstring(x_type,
					"C_Arith_SVT1_SVT2", "x_type");
	int x_has_NAbg = _get_and_check_na_background(x_na_background,
					"C_Arith_SVT1_SVT2", "x_na_background");
	SEXPTYPE y_Rtype = _get_and_check_Rtype_from_Rstring(y_type,
					"C_Arith_SVT1_SVT2", "y_type");
        int y_has_NAbg = _get_and_check_na_background(y_na_background,
					"C_Arith_SVT1_SVT2", "y_na_background");
	SEXPTYPE ans_Rtype = _get_and_check_Rtype_from_Rstring(ans_type,
					"C_Arith_SVT1_SVT2", "ans_type");

	int opcode = _get_Arith_opcode(op);

	if (!x_has_NAbg && !y_has_NAbg && opcode != ADD_OPCODE &&
					  opcode != SUB_OPCODE &&
					  opcode != MULT_OPCODE)
	{
		error("\"%s\" is not supported between SVT_SparseArray "
		      "objects", CHAR(STRING_ELT(op, 0)));
	}
	int dim0 = INTEGER(x_dim)[0];
	SparseVec buf_sv = alloc_SparseVec(ans_Rtype, dim0,
					   x_has_NAbg || y_has_NAbg);

	int ovflow = 0;
	SEXP ans = REC_Arith_SVT1_SVT2(opcode,
				x_SVT, x_Rtype, x_has_NAbg,
				y_SVT, y_Rtype, y_has_NAbg,
				INTEGER(x_dim), LENGTH(x_dim),
				&buf_sv, &ovflow);
	if (ans != R_NilValue)
		PROTECT(ans);
	if (ovflow)
		warning("NAs produced by integer overflow");
	if (ans != R_NilValue)
		UNPROTECT(1);
	return ans;
}

