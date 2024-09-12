/****************************************************************************
 *                'Logic' operations on SparseArray objects                 *
 ****************************************************************************/
#include "SparseArray_Logic_methods.h"

#include "Rvector_utils.h"          /* _get_Rtype_from_Rstring() */
#include "SparseVec.h"
#include "SparseVec_Logic.h"
#include "leaf_utils.h"

#include <string.h>  /* for memcmp() */


/****************************************************************************
 * 'Logic' operations on the tree leaves
 */

static SEXP Logic_leaf1_leaf2(int opcode,
		SEXP leaf1, SEXPTYPE Rtype1, SEXP leaf2, SEXPTYPE Rtype2,
		int dim0,
		int *nzvals_buf, int *nzoffs_buf)
{
	if (leaf1 == R_NilValue || leaf2 == R_NilValue) {
		if (opcode == AND_OPCODE)
			return R_NilValue;
		return leaf1 == R_NilValue ? leaf2 : leaf1;
	}
	const SparseVec sv1 = leaf2SV(leaf1, Rtype1, dim0, 0);
	const SparseVec sv2 = leaf2SV(leaf2, Rtype2, dim0, 0);
	int buf_len = _Logic_intSV_intSV(opcode, &sv1, &sv2,
					 nzvals_buf, nzoffs_buf);
	return _make_leaf_from_two_arrays(LGLSXP,
					  nzvals_buf, nzoffs_buf, buf_len);
}


/****************************************************************************
 * Recursive tree traversals
 */

static SEXP REC_Logic_SVT1_SVT2(int opcode,
		SEXP SVT1, SEXPTYPE Rtype1, SEXP SVT2, SEXPTYPE Rtype2,
		const int *dim, int ndim,
		int *nzvals_buf, int *nzoffs_buf)
{
	if (SVT1 == R_NilValue || SVT2 == R_NilValue) {
		if (opcode == AND_OPCODE)
			return R_NilValue;
		return SVT1 == R_NilValue ? SVT2 : SVT1;
	}

	if (ndim == 1)
		/* 'SVT1' and 'SVT2' are leaves (i.e. 1D SVTs).
		    They should not be NULL. */
		return Logic_leaf1_leaf2(opcode, SVT1, Rtype1, SVT2, Rtype2,
					 dim[0],
					 nzvals_buf, nzoffs_buf);

	/* Each of 'SVT1' and 'SVT2' is a list. */
	int ans_len = dim[ndim - 1];
	SEXP ans = PROTECT(NEW_LIST(ans_len));
	int is_empty = 1;
	for (int i = 0; i < ans_len; i++) {
		SEXP subSVT1 = VECTOR_ELT(SVT1, i);
		SEXP subSVT2 = VECTOR_ELT(SVT2, i);
		SEXP ans_elt = REC_Logic_SVT1_SVT2(opcode,
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

static void check_array_conformability(SEXP x_dim, SEXP y_dim)
{
	int ndim = LENGTH(x_dim);
	if (ndim != LENGTH(y_dim) ||
	    memcmp(INTEGER(x_dim), INTEGER(y_dim), sizeof(int) * ndim) != 0)
		error("non-conformable arrays");
	return;
}

/* --- .Call ENTRY POINT --- */
SEXP C_Logic_SVT1_SVT2(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		       SEXP y_dim, SEXP y_type, SEXP y_SVT,
		       SEXP op)
{
	check_array_conformability(x_dim, y_dim);
	SEXPTYPE x_Rtype = _get_Rtype_from_Rstring(x_type);
	SEXPTYPE y_Rtype = _get_Rtype_from_Rstring(y_type);
	if (x_Rtype == 0 || y_Rtype == 0)
		error("SparseArray internal error in "
		      "C_Logic_SVT1_SVT2():\n"
		      "    invalid 'x_type' or 'y_type' value");
	int opcode = _get_Logic_opcode(op);
	int dim0 = INTEGER(x_dim)[0];
	int *nzvals_buf = (int *) R_alloc(dim0, sizeof(int));
	int *nzoffs_buf = (int *) R_alloc(dim0, sizeof(int));
	return REC_Logic_SVT1_SVT2(opcode, x_SVT, x_Rtype, y_SVT, y_Rtype,
				   INTEGER(x_dim), LENGTH(x_dim),
				   nzvals_buf, nzoffs_buf);
}

