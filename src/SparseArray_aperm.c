/****************************************************************************
 *                  Transposition of a SparseArray object                   *
 ****************************************************************************/
#include "SparseArray_aperm.h"

#include "Rvector_utils.h"
#include "leaf_vector_utils.h"
#include "SBT_utils.h"


static SEXP transpose_SVT(SEXP SVT, const int *dim, int ndim, SEXPTYPE Rtype)
{
	int ans_dim[2], nrow, ncol, j, coords0[2], lv_len, k;
	SEXP ans, lv, lv_offs, lv_vals;
	const int *lv_offs_p, *lv_vals_p;

	if (ndim != 2)
		error("object to transpose must have exactly 2 dimensions");

	if (SVT == R_NilValue)
		return SVT;

	ans_dim[1] = nrow = dim[0];
	ans_dim[0] = ncol = dim[1];
	ans = PROTECT(NEW_LIST(nrow));
	for (j = 0; j < ncol; j++) {
		lv = VECTOR_ELT(SVT, j);
		if (lv == R_NilValue)
			continue;
		coords0[0] = j;
		lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
		lv_offs_p = INTEGER(lv_offs);
		lv_vals_p = INTEGER(lv_vals);
		for (k = 0; k < lv_len; k++) {
			coords0[1] = *lv_offs_p;
			_push_int_to_SBT(ans, ans_dim, 2, coords0, *lv_vals_p);
			lv_offs_p++;
			lv_vals_p++;
		}
	}
	_SBT2SVT(ans, dim, ndim, Rtype);
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_transpose_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT)
{
	SEXPTYPE Rtype;

	Rtype = _get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("SparseArray internal error in "
		      "C_transpose_SVT():\n"
		      "    SVT_SparseArray object has invalid type");

	return transpose_SVT(x_SVT, INTEGER(x_dim), LENGTH(x_dim), Rtype);
}

