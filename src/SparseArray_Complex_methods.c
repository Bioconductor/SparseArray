/****************************************************************************
 *                 Complex methods for SparseArray objects                  *
 ****************************************************************************/
#include "SparseArray_Complex_methods.h"

#include "leaf_vector_utils.h"


/* --- .Call ENTRY POINT --- */
SEXP C_Complex_SVT(SEXP z_dim, SEXP z_type, SEXP z_SVT, SEXP op)
{
	SEXP ans, ans_elt;

	error("not implemented yet, sorry!");
	ans = PROTECT(NEW_LIST(2));

	ans_elt = PROTECT(mkChar("complex"));
	SET_VECTOR_ELT(ans, 0, ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
}

