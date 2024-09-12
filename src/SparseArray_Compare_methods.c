/****************************************************************************
 *               'Compare' operations on SparseArray objects                *
 ****************************************************************************/
#include "SparseArray_Compare_methods.h"

#include "argcheck_utils.h"
#include "SparseVec.h"
#include "SparseVec_Compare.h"
#include "leaf_utils.h"


static SEXP make_noNA_logical_leaf(SEXP nzoffs)
{
	if (LACUNAR_MODE_IS_ON)
		return _make_lacunar_leaf(nzoffs);
	SEXP nzvals = PROTECT(_new_Rvector1(LGLSXP, LENGTH(nzoffs)));
	SEXP ans = zip_leaf(nzvals, nzoffs, 0);
	UNPROTECT(1);
	return ans;
}

static SEXP make_logical_leaf_with_single_shared_int(int na_background,
		int shared_int, SEXP nzoffs)
{
	if (na_background) {
		/* Sanity check. */
		if (shared_int == NA_INTEGER)
			error("SparseArray internal error in "
			      "make_logical_leaf_with_single_shared_int():\n"
			      "    shared_int == NA_INTEGER");
		return _make_leaf_with_single_shared_nzval(LGLSXP,
						&shared_int, nzoffs);
	}
	/* Sanity check. */
	if (shared_int != int1)
		error("SparseArray internal error in "
		      "make_logical_leaf_with_single_shared_int():\n"
		      "    shared_int != int1");
	return make_noNA_logical_leaf(nzoffs);
}


/****************************************************************************
 * 'Compare' operations on the tree leaves
 */

static SEXP Compare_leaf1_zero(int opcode,
		SEXP leaf1, SEXPTYPE Rtype1,
		int dim0,
		int *nzvals_buf, int *nzoffs_buf)
{
	const SparseVec sv1 = leaf2SV(leaf1, Rtype1, dim0, 0);
	int buf_len = _Compare_sv1_zero(opcode, &sv1,
					nzvals_buf, nzoffs_buf);
	if (buf_len == PROPAGATE_NZOFFS)
		return make_logical_leaf_with_single_shared_int(
				sv1.na_background, nzvals_buf[0],
				get_leaf_nzoffs(leaf1));
	return _make_leaf_from_two_arrays(LGLSXP,
					  nzvals_buf, nzoffs_buf, buf_len);
}

static SEXP Compare_leaf1_scalar(int opcode,
		SEXP leaf1, SEXPTYPE Rtype1, int na_background1, SEXP scalar,
		int dim0,
		int *nzvals_buf, int *nzoffs_buf)
{
	const SparseVec sv1 = leaf2SV(leaf1, Rtype1, dim0, na_background1);
	int buf_len = _Compare_sv1_scalar(opcode, &sv1, scalar,
					  nzvals_buf, nzoffs_buf);
	if (buf_len == PROPAGATE_NZOFFS)
		return make_logical_leaf_with_single_shared_int(
				sv1.na_background, nzvals_buf[0],
				get_leaf_nzoffs(leaf1));
	return _make_leaf_from_two_arrays(LGLSXP,
					  nzvals_buf, nzoffs_buf, buf_len);
}

static SEXP Compare_leaf1_leaf2(int opcode,
		SEXP leaf1, SEXPTYPE Rtype1, SEXP leaf2, SEXPTYPE Rtype2,
		int dim0,
		int *nzvals_buf, int *nzoffs_buf)
{
	if (leaf1 == R_NilValue) {
		if (leaf2 == R_NilValue)
			return R_NilValue;
		return Compare_leaf1_zero(flip_opcode(opcode), leaf2, Rtype2,
					  dim0,
					  nzvals_buf, nzoffs_buf);
	}
	if (leaf2 == R_NilValue)
		return Compare_leaf1_zero(opcode, leaf1, Rtype1, dim0,
					  nzvals_buf, nzoffs_buf);
	const SparseVec sv1 = leaf2SV(leaf1, Rtype1, dim0, 0);
	const SparseVec sv2 = leaf2SV(leaf2, Rtype2, dim0, 0);
	int buf_len = _Compare_sv1_sv2(opcode, &sv1, &sv2,
				       nzvals_buf, nzoffs_buf);
	return _make_leaf_from_two_arrays(LGLSXP,
					  nzvals_buf, nzoffs_buf, buf_len);
}


/****************************************************************************
 * Recursive tree traversals
 */

static SEXP REC_Compare_SVT1_scalar(int opcode,
		SEXP SVT1, SEXPTYPE Rtype1, int na_background1, SEXP scalar,
		const int *dim, int ndim,
		void *nzvals_buf, int *nzoffs_buf)
{
	if (SVT1 == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* 'SVT1' is a leaf (i.e. 1D SVT). */
		return Compare_leaf1_scalar(opcode,
					SVT1, Rtype1, na_background1, scalar,
					dim[0],
					nzvals_buf, nzoffs_buf);
	}

	/* 'SVT1' is a list. */
	int ans_len = dim[ndim - 1];
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	int is_empty = 1;
	for (int i = 0; i < ans_len; i++) {
		SEXP subSVT1 = VECTOR_ELT(SVT1, i);
		SEXP ans_elt = REC_Compare_SVT1_scalar(opcode,
					subSVT1, Rtype1, na_background1, scalar,
					dim, ndim - 1,
					nzvals_buf, nzoffs_buf);
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

static SEXP REC_Compare_SVT1_SVT2(int opcode,
		SEXP SVT1, SEXPTYPE Rtype1, SEXP SVT2, SEXPTYPE Rtype2,
		const int *dim, int ndim,
		int *nzvals_buf, int *nzoffs_buf)
{
	if (SVT1 == R_NilValue && SVT2 == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* 'SVT1' and 'SVT2' are leaves (i.e. 1D SVTs).
		   They cannot both be NULL. */
		return Compare_leaf1_leaf2(opcode, SVT1, Rtype1, SVT2, Rtype2,
					   dim[0],
					   nzvals_buf, nzoffs_buf);
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
		SEXP ans_elt = REC_Compare_SVT1_SVT2(opcode,
					subSVT1, Rtype1, subSVT2, Rtype2,
					dim, ndim - 1,
					nzvals_buf, nzoffs_buf);
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
SEXP C_Compare_SVT1_v2(
		SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP x_na_background,
		SEXP v2, SEXP op)
{
	SEXPTYPE x_Rtype = _get_and_check_Rtype_from_Rstring(x_type,
					"C_Compare_SVT1_v2", "x_type");

	int x_has_NAbg = _get_and_check_na_background(x_na_background,
					"C_Compare_SVT1_v2", "x_na_background");

	int opcode = _get_Compare_opcode(op);
	int dim0 = INTEGER(x_dim)[0];
	int *nzvals_buf = (int *) R_alloc(dim0, sizeof(int));
	int *nzoffs_buf = (int *) R_alloc(dim0, sizeof(int));
	return REC_Compare_SVT1_scalar(opcode,
			x_SVT, x_Rtype, x_has_NAbg, v2,
			INTEGER(x_dim), LENGTH(x_dim),
			nzvals_buf, nzoffs_buf);
}

/* --- .Call ENTRY POINT --- */
SEXP C_Compare_SVT1_SVT2(
		SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP x_na_background,
		SEXP y_dim, SEXP y_type, SEXP y_SVT, SEXP y_na_background,
		SEXP op)
{
	_check_array_conformability(x_dim, y_dim);
	SEXPTYPE x_Rtype = _get_and_check_Rtype_from_Rstring(x_type,
					"C_Compare_SVT1_SVT2", "x_type");
	SEXPTYPE y_Rtype = _get_and_check_Rtype_from_Rstring(y_type,
					"C_Compare_SVT1_SVT2", "y_type");

	int opcode = _get_Compare_opcode(op);
	if (opcode != NE_OPCODE &&
	    opcode != LT_OPCODE &&
	    opcode != GT_OPCODE)
	{
		error("\"%s\" is not supported between SVT_SparseArray "
		      "objects", CHAR(STRING_ELT(op, 0)));
	}
	int dim0 = INTEGER(x_dim)[0];
	int *nzvals_buf = (int *) R_alloc(dim0, sizeof(int));
	int *nzoffs_buf = (int *) R_alloc(dim0, sizeof(int));
	return REC_Compare_SVT1_SVT2(opcode, x_SVT, x_Rtype, y_SVT, y_Rtype,
				     INTEGER(x_dim), LENGTH(x_dim),
				     nzvals_buf, nzoffs_buf);
}

