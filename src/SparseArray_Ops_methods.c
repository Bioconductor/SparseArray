/****************************************************************************
 *                   Ops methods for SparseArray objects                    *
 ****************************************************************************/
#include "SparseArray_Ops_methods.h"

#include "Rvector_utils.h"
#include "leaf_vector_Ops.h"


/* --- .Call ENTRY POINT --- */
SEXP C_SVT_Arith(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		 SEXP y_dim, SEXP y_type, SEXP y_SVT,
		 SEXP op)
{
	SEXPTYPE x_Rtype, y_Rtype;
	int opcode;

	x_Rtype = _get_Rtype_from_Rstring(x_type);
	y_Rtype = _get_Rtype_from_Rstring(y_type);
	if (x_Rtype == 0 || y_Rtype == 0)
		error("SparseArray internal error in "
		      "C_SVT_Arith():\n"
                      "    invalid 'x_type' or 'y_type' value");
	opcode = _get_Arith_opcode(op, x_Rtype, y_Rtype);
	error("not implemented yet, sorry!");
	return R_NilValue;
}

/* --- .Call ENTRY POINT --- */
SEXP C_SVT_Compare(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		   SEXP y_dim, SEXP y_type, SEXP y_SVT,
		   SEXP op)
{
	SEXPTYPE x_Rtype, y_Rtype;
	int opcode;

	x_Rtype = _get_Rtype_from_Rstring(x_type);
	y_Rtype = _get_Rtype_from_Rstring(y_type);
	if (x_Rtype == 0 || y_Rtype == 0)
		error("SparseArray internal error in "
		      "C_SVT_Compare():\n"
                      "    invalid 'x_type' or 'y_type' value");
	opcode = _get_Compare_opcode(op, x_Rtype, y_Rtype);
	error("not implemented yet, sorry!");
	return R_NilValue;
}

/* --- .Call ENTRY POINT --- */
SEXP C_SVT_Logic(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		 SEXP y_dim, SEXP y_type, SEXP y_SVT,
		 SEXP op)
{
	SEXPTYPE x_Rtype, y_Rtype;
	int opcode;

	x_Rtype = _get_Rtype_from_Rstring(x_type);
	y_Rtype = _get_Rtype_from_Rstring(y_type);
	if (x_Rtype == 0 || y_Rtype == 0)
		error("SparseArray internal error in "
		      "C_SVT_Logic():\n"
                      "    invalid 'x_type' or 'y_type' value");
	opcode = _get_Logic_opcode(op, x_Rtype, y_Rtype);
	error("not implemented yet, sorry!");
	return R_NilValue;
}

