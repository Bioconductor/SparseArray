/****************************************************************************
 *  Low-level stuff shared by COO_SparseArray and SVT_SparseArray objects   *
 ****************************************************************************/
#include "SparseArray_class.h"

#include "Rvector_utils.h"

int _coercion_can_introduce_zeros(SEXPTYPE from_Rtype, SEXPTYPE to_Rtype)
{
	if (to_Rtype == from_Rtype)
		return 0;
	if (to_Rtype == RAWSXP || from_Rtype == STRSXP || from_Rtype == VECSXP)
		return 1;
	if (from_Rtype == REALSXP)
		return to_Rtype == INTSXP;
	if (from_Rtype == CPLXSXP)
		return to_Rtype == INTSXP || to_Rtype == REALSXP;
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_coercion_can_introduce_zeros(SEXP from_type, SEXP to_type)
{
	SEXPTYPE from_Rtype, to_Rtype;

	from_Rtype = _get_Rtype_from_Rstring(from_type);
	to_Rtype = _get_Rtype_from_Rstring(to_type);
	if (from_Rtype == 0 || to_Rtype == 0)
		error("'from_type' and 'to_type' must be valid "
		      "vector types specified\n  as single strings");
	return ScalarLogical(_coercion_can_introduce_zeros(from_Rtype,
							   to_Rtype));
}

