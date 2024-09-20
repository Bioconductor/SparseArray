/****************************************************************************
 *                'Logic' operations on SparseArray objects                 *
 ****************************************************************************/
#include "SparseArray_Logic_methods.h"

#include "argcheck_utils.h"
#include "SparseVec.h"
#include "SparseVec_Logic.h"
#include "leaf_utils.h"


/****************************************************************************
 * 'Logic' operations on the tree leaves
 */

static void INPLACE_logical_neg_naleaf(SEXP naleaf, SEXPTYPE Rtype)
{
	SEXP nzvals, nzoffs;
	int nzcount = unzip_leaf(naleaf, &nzvals, &nzoffs);
	if (nzvals == R_NilValue) {  /* input leaf is lacunar */
		SEXP new_nzvals = PROTECT(
			_new_Rvector0(Rtype, (R_xlen_t) nzcount)
		);
		replace_leaf_nzvals(naleaf, new_nzvals);
		UNPROTECT(1);
		return;
	}
	/* input leaf is regular */
	if (Rtype != LGLSXP && Rtype != INTSXP)
		error("SparseArray internal error in "
		      "INPLACE_logical_neg_naleaf():\n"
		      "    logical negation (\"!\") of an NaArray object "
		      "of type \"%s\" is not supported", type2char(Rtype));
	int *nzvals_p = INTEGER(nzvals);
	int all_ones = 1;
	for (int k = 0; k < nzcount; k++) {
		/* The input leaf is assumed to be an NA-leaf so 'nzvals_p[k]'
		   should never be NA. So no need to check for NA before
		   negating 'nzvals_p[k]' with ! */
		if ((nzvals_p[k] = !nzvals_p[k]) != int1)
			all_ones = 0;
	}
	if (LACUNAR_MODE_IS_ON && all_ones)
		replace_leaf_nzvals(naleaf, R_NilValue);
	return ;
}

static SEXP Logic_leaf1_na(int opcode,
		SEXP leaf1, SEXPTYPE Rtype1, int na_background1,
		SEXPTYPE Rtype2,
		SparseVec *buf_sv)
{
	if (leaf1 == R_NilValue)
		error("SparseArray internal error in "
		      "Logic_leaf1_na():\n"
		      "    'leaf1' cannot be NULL");
	const SparseVec sv1 = leaf2SV(leaf1, Rtype1,
				      buf_sv->len, na_background1);
	_Logic_intSV_na(opcode, &sv1, Rtype2, buf_sv);
	if (buf_sv->nzcount == PROPAGATE_NZOFFS)
		return _make_leaf_with_single_shared_nzval(
					      buf_sv->Rtype, buf_sv->nzvals,
					      get_leaf_nzoffs(leaf1));
	return SV2leaf(buf_sv);
}

static SEXP Logic_leaf1_leaf2(int opcode,
		SEXP leaf1, SEXPTYPE Rtype1, int na_background1,
		SEXP leaf2, SEXPTYPE Rtype2, int na_background2,
		SparseVec *buf_sv)
{
	if (leaf1 == R_NilValue) {
		if (!na_background1)
			error("SparseArray internal error in "
			      "Logic_leaf1_leaf2():\n"
			      "    'na_background1' is expected to be TRUE");
		return Logic_leaf1_na(opcode,
				      leaf2, Rtype2, na_background2,
				      Rtype1,
				      buf_sv);
	}
	if (leaf2 == R_NilValue) {
		if (!na_background2)
			error("SparseArray internal error in "
			      "Logic_leaf1_leaf2():\n"
			      "    'na_background2' is expected to be TRUE");
		return Logic_leaf1_na(opcode,
				      leaf1, Rtype1, na_background1,
				      Rtype2,
				      buf_sv);
	}
	const SparseVec sv1 = leaf2SV(leaf1, Rtype1,
				      buf_sv->len, na_background1);
	const SparseVec sv2 = leaf2SV(leaf2, Rtype2,
				      buf_sv->len, na_background2);
	_Logic_intSV_intSV(opcode, &sv1, &sv2, buf_sv);
	return SV2leaf(buf_sv);
}


/****************************************************************************
 * Recursive tree traversals
 */

static void REC_logical_neg_NaSVT(SEXP NaSVT, SEXPTYPE Rtype,
				  const int *dim, int ndim)
{
	if (NaSVT == R_NilValue)
		return;
	if (ndim == 1) {
		/* 'NaSVT' is a leaf (i.e. 1D SVT). */
		INPLACE_logical_neg_naleaf(NaSVT, Rtype);
		return;
	}
	int NaSVT_len = dim[ndim - 1];
	for (int i = 0; i < NaSVT_len; i++)
		REC_logical_neg_NaSVT(VECTOR_ELT(NaSVT, i), Rtype,
				      dim, ndim - 1);
	return;
}

static SEXP REC_Logic_SVT1_na(int opcode,
		SEXP SVT1, SEXPTYPE Rtype1, int na_background1,
		const int *dim, int ndim,
		SparseVec *buf_sv)
{
	if (SVT1 == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* 'SVT1' is a leaf (i.e. 1D SVT). */
		return Logic_leaf1_na(opcode,
				      SVT1, Rtype1, na_background1, LGLSXP,
				      buf_sv);
	}

	/* 'SVT1' is a list. */
	int ans_len = dim[ndim - 1];
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	int is_empty = 1;
	for (int i = 0; i < ans_len; i++) {
		SEXP subSVT1 = VECTOR_ELT(SVT1, i);
		SEXP ans_elt = REC_Logic_SVT1_na(opcode,
					subSVT1, Rtype1, na_background1,
					dim, ndim - 1,
					buf_sv);
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

static SEXP REC_Logic_SVT1_SVT2(int opcode,
		SEXP SVT1, SEXPTYPE Rtype1, int na_background1,
		SEXP SVT2, SEXPTYPE Rtype2, int na_background2,
		const int *dim, int ndim,
		SparseVec *buf_sv)
{
	if (SVT1 == R_NilValue && SVT2 == R_NilValue)
		return R_NilValue;

	if (!na_background1 && SVT1 == R_NilValue)
		return opcode == OR_OPCODE ? SVT2 : R_NilValue;
	if (!na_background2 && SVT2 == R_NilValue)
		return opcode == OR_OPCODE ? SVT1 : R_NilValue;

	if (ndim == 1)
		/* 'SVT1' and 'SVT2' are leaves (i.e. 1D SVTs). */
		return Logic_leaf1_leaf2(opcode,
					 SVT1, Rtype1, na_background1,
					 SVT2, Rtype2, na_background2,
					 buf_sv);

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
		SEXP ans_elt = REC_Logic_SVT1_SVT2(opcode,
					subSVT1, Rtype1, na_background1,
					subSVT2, Rtype2, na_background2,
					dim, ndim - 1,
					buf_sv);
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
SEXP C_logical_neg_NaSVT(SEXP x_dim, SEXP x_type, SEXP x_NaSVT)
{
	SEXPTYPE x_Rtype = _get_and_check_Rtype_from_Rstring(x_type,
				"C_logical_neg_NaSVT", "x_type");
	SEXP ans = PROTECT(duplicate(x_NaSVT));
	REC_logical_neg_NaSVT(ans, x_Rtype, INTEGER(x_dim), LENGTH(x_dim));
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_Logic_NaSVT1_na(SEXP x_dim, SEXP x_type, SEXP x_NaSVT, SEXP op)
{
	SEXPTYPE x_Rtype = _get_and_check_Rtype_from_Rstring(x_type,
				"C_logical_neg_NaSVT", "x_type");

	int opcode = _get_Logic_opcode(op);

	int dim0 = INTEGER(x_dim)[0];
	SparseVec buf_sv = alloc_SparseVec(LGLSXP, dim0, 1);
	return REC_Logic_SVT1_na(opcode, x_NaSVT, x_Rtype, 1,
				 INTEGER(x_dim), LENGTH(x_dim),
				 &buf_sv);
}

/* --- .Call ENTRY POINT --- */
SEXP C_Logic_SVT1_SVT2(
		SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP x_na_background,
		SEXP y_dim, SEXP y_type, SEXP y_SVT, SEXP y_na_background,
		SEXP op)
{
	_check_array_conformability(x_dim, y_dim);
	SEXPTYPE x_Rtype = _get_and_check_Rtype_from_Rstring(x_type,
				"C_Logic_SVT1_SVT2", "x_type");
	int x_has_NAbg = _get_and_check_na_background(x_na_background,
				"C_Logic_SVT1_SVT2", "x_na_background");
	SEXPTYPE y_Rtype = _get_and_check_Rtype_from_Rstring(y_type,
				"C_Logic_SVT1_SVT2", "y_type");
	int y_has_NAbg = _get_and_check_na_background(y_na_background,
				"C_Logic_SVT1_SVT2", "y_na_background");

	int opcode = _get_Logic_opcode(op);

	int dim0 = INTEGER(x_dim)[0];
	int out_na_background = 0;
	if (x_has_NAbg && y_has_NAbg) {
		out_na_background = 1;
	} else if (x_has_NAbg || y_has_NAbg) {
		out_na_background = opcode == OR_OPCODE;
	}
	SparseVec buf_sv = alloc_SparseVec(LGLSXP, dim0, out_na_background);
	return REC_Logic_SVT1_SVT2(opcode,
				x_SVT, x_Rtype, x_has_NAbg,
				y_SVT, y_Rtype, y_has_NAbg,
				INTEGER(x_dim), LENGTH(x_dim),
				&buf_sv);
}

