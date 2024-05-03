/****************************************************************************
 *                  'Ops' methods for SparseArray objects                   *
 ****************************************************************************/
#include "SparseArray_Ops_methods.h"

#include "Rvector_utils.h"          /* _get_Rtype_from_Rstring() */
#include "SparseVec.h"
#include "SparseVec_Arith.h"
#include "SparseVec_Compare.h"
#include "SparseVec_Logic.h"
#include "leaf_utils.h"
#include "SVT_SparseArray_class.h"  /* for _coerce_SVT() */

#include <string.h>  /* for memcmp() */


/****************************************************************************
 * 'Arith' operations on the tree leaves
 */

/* Assumes that 'ans_Rtype' is equal or bigger than the type of the nonzero
   values in 'leaf'. Performs **in-place** replacement if 'ans_Rtype' is 0!
   Note that this could also have been achieved by calling:

     Arith_leaf1_scalar(MULT_OPCODE, leaf, -1, ...)

   but unary_minus_leaf() takes a lot of shortcuts so is A LOT more
   efficient. */
static SEXP unary_minus_leaf(SEXP leaf, SEXPTYPE ans_Rtype)
{
	SEXP nzoffs, leaf_nzvals, ans_nzvals;

	int nzcount = unzip_leaf(leaf, &nzoffs, &leaf_nzvals);
	if (ans_Rtype == 0) {
		/* In-place replacement! */
		ans_nzvals = leaf_nzvals;
	} else {
		ans_nzvals = PROTECT(allocVector(ans_Rtype, nzcount));
        }
	const char *errmsg = _unary_minus_Rvector(leaf_nzvals, ans_nzvals);
	if (errmsg != NULL) {
		if (ans_Rtype != 0)
			UNPROTECT(1);
		error("%s", errmsg);
	}
	if (ans_Rtype == 0)
		return leaf;
	SEXP ans = zip_leaf(nzoffs, ans_nzvals);
	UNPROTECT(1);
	return ans;
}

static SEXP Arith_leaf1_scalar(int opcode, SEXP leaf1, SEXP scalar, int dim0,
		SEXPTYPE ans_Rtype,
		int *nzoffs_buf, void *nzvals_buf, int *ovflow)
{
	const SparseVec sv1 = leaf2SV(leaf1, dim0);
	int nzcount = _Arith_sv1_scalar(opcode, &sv1, scalar, ans_Rtype,
					nzoffs_buf, nzvals_buf, ovflow);
	return _make_leaf_from_bufs(ans_Rtype, nzoffs_buf, nzvals_buf, nzcount);
}

static SEXP Arith_leaf1_leaf2(int opcode, SEXP leaf1, SEXP leaf2, int dim0,
		SEXPTYPE ans_Rtype,
		int *nzoffs_buf, void *nzvals_buf, int *ovflow)
{
	int nzcount;

	if (leaf1 == R_NilValue) {
		if (leaf2 == R_NilValue)
			error("SparseArray internal error in "
			      "Arith_leaf1_leaf2():\n"
			      "    'leaf1' and 'leaf2' cannot both be NULL");
		if (opcode == SUB_OPCODE)
			return unary_minus_leaf(leaf2, ans_Rtype);
		if (opcode == MULT_OPCODE) {
			const SparseVec sv2 = leaf2SV(leaf2, dim0);
			nzcount = _mult_SV_zero(&sv2, ans_Rtype,
						nzoffs_buf, nzvals_buf);
		} else {
			error("SparseArray internal error in "
			      "Arith_leaf1_leaf2():\n"
			      "    'op' must be \"-\" or \"*\" "
			      "when 'leaf1' is NULL");
		}
	} else if (leaf2 == R_NilValue) {
		if (opcode == MULT_OPCODE) {
			const SparseVec sv1 = leaf2SV(leaf1, dim0);
			nzcount = _mult_SV_zero(&sv1, ans_Rtype,
						nzoffs_buf, nzvals_buf);
		} else {
			error("SparseArray internal error in "
			      "_Arith_leaf1_leaf2():\n"
			      "    'op' must be \"*\" when 'leaf2' is NULL");
		}
	} else {
		const SparseVec sv1 = leaf2SV(leaf1, dim0);
		const SparseVec sv2 = leaf2SV(leaf2, dim0);
		nzcount = _Arith_sv1_sv2(opcode, &sv1, &sv2, ans_Rtype,
					 nzoffs_buf, nzvals_buf, ovflow);
	}
	return _make_leaf_from_bufs(ans_Rtype, nzoffs_buf, nzvals_buf, nzcount);
}


/****************************************************************************
 * 'Compare' operations on the tree leaves
 */

static SEXP Compare_leaf1_zero(int opcode, SEXP leaf1, int dim0,
		int *nzoffs_buf, int *nzvals_buf)
{
	const SparseVec sv1 = leaf2SV(leaf1, dim0);
	int nzcount = _Compare_sv1_zero(opcode, &sv1,
					nzoffs_buf, nzvals_buf);
	return _make_leaf_from_bufs(LGLSXP, nzoffs_buf, nzvals_buf, nzcount);
}

static SEXP Compare_leaf1_scalar(int opcode, SEXP leaf1, SEXP scalar, int dim0,
		int *nzoffs_buf, int *nzvals_buf)
{
	const SparseVec sv1 = leaf2SV(leaf1, dim0);
	int nzcount = _Compare_sv1_scalar(opcode, &sv1, scalar,
					  nzoffs_buf, nzvals_buf);
	return _make_leaf_from_bufs(LGLSXP, nzoffs_buf, nzvals_buf, nzcount);
}

static SEXP Compare_leaf1_leaf2(int opcode, SEXP leaf1, SEXP leaf2, int dim0,
		int *nzoffs_buf, int *nzvals_buf)
{
	if (leaf1 == R_NilValue) {
		if (leaf2 == R_NilValue)
			return R_NilValue;
		return Compare_leaf1_zero(flip_opcode(opcode), leaf2, dim0,
					  nzoffs_buf, nzvals_buf);
	}
	if (leaf2 == R_NilValue)
		return Compare_leaf1_zero(opcode, leaf1, dim0,
					  nzoffs_buf, nzvals_buf);
	const SparseVec sv1 = leaf2SV(leaf1, dim0);
	const SparseVec sv2 = leaf2SV(leaf2, dim0);
	int nzcount = _Compare_sv1_sv2(opcode, &sv1, &sv2,
				       nzoffs_buf, nzvals_buf);
	return _make_leaf_from_bufs(LGLSXP, nzoffs_buf, nzvals_buf, nzcount);
}


/****************************************************************************
 * 'Logic' operations on the tree leaves
 */

static SEXP Logic_leaf1_leaf2(int opcode, SEXP leaf1, SEXP leaf2, int dim0,
		int *nzoffs_buf, int *nzvals_buf)
{
	if (leaf1 == R_NilValue || leaf2 == R_NilValue) {
		if (opcode == AND_OPCODE)
			return R_NilValue;
		return leaf1 == R_NilValue ? leaf2 : leaf1;
	}
	const SparseVec sv1 = leaf2SV(leaf1, dim0);
	const SparseVec sv2 = leaf2SV(leaf2, dim0);
	int nzcount = _Logic_intSV_intSV(opcode, &sv1, &sv2,
					 nzoffs_buf, nzvals_buf);
	return _make_leaf_from_bufs(LGLSXP, nzoffs_buf, nzvals_buf, nzcount);
}


/****************************************************************************
 * Recursive tree traversals
 */

static void REC_unary_minus_SVT(SEXP SVT, const int *dim, int ndim)
{
	if (SVT == R_NilValue)
		return;
	if (ndim == 1) {
		/* 'SVT' is a leaf (i.e. 1D SVT). */
		unary_minus_leaf(SVT, 0);
		return;
	}
	int SVT_len = dim[ndim - 1];
	for (int i = 0; i < SVT_len; i++)
		REC_unary_minus_SVT(VECTOR_ELT(SVT, i), dim, ndim - 1);
	return;
}

static SEXP REC_Arith_SVT1_v2(int opcode, SEXP SVT1, SEXP v2,
			      const int *dim, int ndim,
			      SEXPTYPE ans_Rtype,
			      int *nzoffs_buf, void *nzvals_buf, int *ovflow)
{
	if (SVT1 == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* 'SVT1' is a leaf (i.e. 1D SVT). */
		return Arith_leaf1_scalar(opcode, SVT1, v2, dim[0],
					  ans_Rtype,
					  nzoffs_buf, nzvals_buf, ovflow);
	}

	/* 'SVT1' is a list. */
	int ans_len = dim[ndim - 1];
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	int is_empty = 1;
	for (int i = 0; i < ans_len; i++) {
		SEXP subSVT1 = VECTOR_ELT(SVT1, i);
		SEXP ans_elt = REC_Arith_SVT1_v2(opcode, subSVT1, v2,
					dim, ndim - 1,
					ans_Rtype,
					nzoffs_buf, nzvals_buf, ovflow);
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

static SEXP REC_Compare_SVT1_v2(int opcode, SEXP SVT1, SEXP v2,
				const int *dim, int ndim,
				int *nzoffs_buf, void *nzvals_buf)
{
	if (SVT1 == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* 'SVT1' is a leaf (i.e. 1D SVT). */
		return Compare_leaf1_scalar(opcode, SVT1, v2, dim[0],
				nzoffs_buf, nzvals_buf);
	}

	/* 'SVT1' is a list. */
	int ans_len = dim[ndim - 1];
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	int is_empty = 1;
	for (int i = 0; i < ans_len; i++) {
		SEXP subSVT1 = VECTOR_ELT(SVT1, i);
		SEXP ans_elt = REC_Compare_SVT1_v2(opcode, subSVT1, v2,
						   dim, ndim - 1,
						   nzoffs_buf, nzvals_buf);
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
				SEXP SVT1, SEXPTYPE Rtype1,
				SEXP SVT2, SEXPTYPE Rtype2,
				const int *dim, int ndim,
				SEXPTYPE ans_Rtype,
				int *nzoffs_buf, void *nzvals_buf, int *ovflow)
{
	if (SVT1 == R_NilValue) {
		if (SVT2 == R_NilValue)
			return R_NilValue;
		if (opcode == ADD_OPCODE)
			return _coerce_SVT(SVT2, dim, ndim,
					   Rtype2, ans_Rtype, nzoffs_buf);
	} else if (SVT2 == R_NilValue) {
		if (opcode == ADD_OPCODE || opcode == SUB_OPCODE)
			return _coerce_SVT(SVT1, dim, ndim,
					   Rtype1, ans_Rtype, nzoffs_buf);
	}

	if (ndim == 1) {
		/* 'SVT1' and 'SVT2' are "leaf vectors", but:
		   - 'SVT1' can be NULL if 'opcode' is SUB_OPCODE;
		   - either 'SVT1' or 'SVT2' (but not both) can be NULL
		     if 'opcode' is MULT_OPCODE. */
		return Arith_leaf1_leaf2(opcode, SVT1, SVT2, dim[0],
					 ans_Rtype,
					 nzoffs_buf, nzvals_buf, ovflow);
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
					subSVT1, Rtype1,
					subSVT2, Rtype2,
					dim, ndim - 1,
					ans_Rtype,
					nzoffs_buf, nzvals_buf, ovflow);
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

static SEXP REC_Compare_SVT1_SVT2(int opcode, SEXP SVT1, SEXP SVT2,
				  const int *dim, int ndim,
				  int *nzoffs_buf, int *nzvals_buf)
{
	if (SVT1 == R_NilValue && SVT2 == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* 'SVT1' and 'SVT2' are leaves (i.e. 1D SVTs).
		   They cannot both be NULL. */
		return Compare_leaf1_leaf2(opcode, SVT1, SVT2, dim[0],
					   nzoffs_buf, nzvals_buf);
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
		SEXP ans_elt = REC_Compare_SVT1_SVT2(opcode, subSVT1, subSVT2,
						     dim, ndim - 1,
						     nzoffs_buf, nzvals_buf);
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

static SEXP REC_Logic_SVT1_SVT2(int opcode, SEXP SVT1, SEXP SVT2,
				const int *dim, int ndim,
				int *nzoffs_buf, int *nzvals_buf)
{
	if (SVT1 == R_NilValue || SVT2 == R_NilValue) {
		if (opcode == AND_OPCODE)
			return R_NilValue;
		return SVT1 == R_NilValue ? SVT2 : SVT1;
	}

	if (ndim == 1)
		return Logic_leaf1_leaf2(opcode, SVT1, SVT2, dim[0],
					 nzoffs_buf, nzvals_buf);

	/* Each of 'SVT1' and 'SVT2' is a list. */
	int ans_len = dim[ndim - 1];
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	int is_empty = 1;
	for (int i = 0; i < ans_len; i++) {
		SEXP subSVT1 = VECTOR_ELT(SVT1, i);
		SEXP subSVT2 = VECTOR_ELT(SVT2, i);
		SEXP ans_elt = REC_Logic_SVT1_SVT2(opcode, subSVT1, subSVT2,
						   dim, ndim - 1,
						   nzoffs_buf, nzvals_buf);
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

/* --- .Call ENTRY POINT ---
  'x_type' is ignored at the moment. */
SEXP C_unary_minus_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT)
{
	SEXP ans = PROTECT(duplicate(x_SVT));
	REC_unary_minus_SVT(ans, INTEGER(x_dim), LENGTH(x_dim));
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_Arith_SVT1_v2(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP v2,
		     SEXP op, SEXP ans_type)
{
	SEXPTYPE x_Rtype = _get_Rtype_from_Rstring(x_type);
	SEXPTYPE ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (x_Rtype == 0 || ans_Rtype == 0)
		error("SparseArray internal error in "
		      "C_Arith_SVT1_v2():\n"
		      "    invalid 'x_type' or 'ans_type' value");
	int opcode = _get_Arith_opcode(op);
	if (opcode != MULT_OPCODE &&
	    opcode != DIV_OPCODE &&
	    opcode != POW_OPCODE &&
	    opcode != MOD_OPCODE &&
	    opcode != IDIV_OPCODE)
	{
		error("\"%s\" is not supported between an SVT_SparseArray "
		      "object and a numeric vector", CHAR(STRING_ELT(op, 0)));
	}
	int dim0 = INTEGER(x_dim)[0];
	int *nzoffs_buf = (int *) R_alloc(dim0, sizeof(int));
	/* Must be big enough to contain ints or doubles. */
	double *nzvals_buf = (double *) R_alloc(dim0, sizeof(double));
	int ovflow = 0;
	SEXP ans = REC_Arith_SVT1_v2(opcode, x_SVT, v2,
				     INTEGER(x_dim), LENGTH(x_dim),
				     ans_Rtype,
				     nzoffs_buf, nzvals_buf, &ovflow);
	if (ovflow) {
		PROTECT(ans);
		warning("NAs produced by integer overflow");
		UNPROTECT(1);
	}
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_Compare_SVT1_v2(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP v2, SEXP op)
{
	SEXPTYPE x_Rtype = _get_Rtype_from_Rstring(x_type);
	if (x_Rtype == 0)
		error("SparseArray internal error in "
		      "C_Compare_SVT1_v2():\n"
		      "    invalid 'x_type'");
	int opcode = _get_Compare_opcode(op);
	int dim0 = INTEGER(x_dim)[0];
	int *nzoffs_buf = (int *) R_alloc(dim0, sizeof(int));
	int *nzvals_buf = (int *) R_alloc(dim0, sizeof(int));
	return REC_Compare_SVT1_v2(opcode, x_SVT, v2,
				   INTEGER(x_dim), LENGTH(x_dim),
				   nzoffs_buf, nzvals_buf);
}

static void check_array_conformability(SEXP x_dim, SEXP y_dim)
{
	int ndim = LENGTH(x_dim);
	if (ndim != LENGTH(y_dim) ||
	    memcmp(INTEGER(x_dim), INTEGER(y_dim), sizeof(int) * ndim) != 0)
		error("non-conformable arrays");
	return;
}

/* --- .Call ENTRY POINT --- */
SEXP C_Arith_SVT1_SVT2(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		       SEXP y_dim, SEXP y_type, SEXP y_SVT,
		       SEXP op, SEXP ans_type)
{
	check_array_conformability(x_dim, y_dim);
	SEXPTYPE x_Rtype = _get_Rtype_from_Rstring(x_type);
	SEXPTYPE y_Rtype = _get_Rtype_from_Rstring(y_type);
	SEXPTYPE ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (x_Rtype == 0 || y_Rtype == 0 || ans_Rtype == 0)
		error("SparseArray internal error in "
		      "C_Arith_SVT1_SVT2():\n"
		      "    invalid 'x_type', 'y_type', or 'ans_type' value");
	int opcode = _get_Arith_opcode(op);
	if (opcode != ADD_OPCODE &&
	    opcode != SUB_OPCODE &&
	    opcode != MULT_OPCODE)
	{
		error("\"%s\" is not supported between SVT_SparseArray "
		      "objects", CHAR(STRING_ELT(op, 0)));
	}
	int dim0 = INTEGER(x_dim)[0];
	int *nzoffs_buf = (int *) R_alloc(dim0, sizeof(int));
	/* Must be big enough to contain ints or doubles. */
	double *nzvals_buf = (double *) R_alloc(dim0, sizeof(double));
	int ovflow = 0;
	SEXP ans = REC_Arith_SVT1_SVT2(opcode, x_SVT, x_Rtype, y_SVT, y_Rtype,
				       INTEGER(x_dim), LENGTH(x_dim),
				       ans_Rtype,
				       nzoffs_buf, nzvals_buf, &ovflow);
	if (ans != R_NilValue)
		PROTECT(ans);
	if (ovflow)
		warning("NAs produced by integer overflow");
	if (ans != R_NilValue)
		UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT ---
   'x_type' and 'y_type' are ignored. */
SEXP C_Compare_SVT1_SVT2(SEXP x_dim, SEXP x_type, SEXP x_SVT,
			 SEXP y_dim, SEXP y_type, SEXP y_SVT,
			 SEXP op)
{
	check_array_conformability(x_dim, y_dim);
	int opcode = _get_Compare_opcode(op);
	if (opcode != NE_OPCODE &&
	    opcode != LT_OPCODE &&
	    opcode != GT_OPCODE)
	{
		error("\"%s\" is not supported between SVT_SparseArray "
		      "objects", CHAR(STRING_ELT(op, 0)));
	}
	int dim0 = INTEGER(x_dim)[0];
	int *nzoffs_buf = (int *) R_alloc(dim0, sizeof(int));
	int *nzvals_buf = (int *) R_alloc(dim0, sizeof(int));
	return REC_Compare_SVT1_SVT2(opcode, x_SVT, y_SVT,
				     INTEGER(x_dim), LENGTH(x_dim),
				     nzoffs_buf, nzvals_buf);
}

/* --- .Call ENTRY POINT ---
   'x_type' and 'y_type' are ignored. */
SEXP C_Logic_SVT1_SVT2(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		       SEXP y_dim, SEXP y_type, SEXP y_SVT,
		       SEXP op)
{
	check_array_conformability(x_dim, y_dim);
	int opcode = _get_Logic_opcode(op);
	int dim0 = INTEGER(x_dim)[0];
	int *nzoffs_buf = (int *) R_alloc(dim0, sizeof(int));
	int *nzvals_buf = (int *) R_alloc(dim0, sizeof(int));
	return REC_Logic_SVT1_SVT2(opcode, x_SVT, y_SVT,
				   INTEGER(x_dim), LENGTH(x_dim),
				   nzoffs_buf, nzvals_buf);
}

