/****************************************************************************
 *                   Ops methods for SparseArray objects                    *
 ****************************************************************************/
#include "SparseArray_Ops_methods.h"

#include "Rvector_utils.h"
#include "leaf_vector_Ops.h"

#include <string.h>  /* for memcmp() */


static SEXP REC_SVT_Arith(SEXP SVT1, SEXP SVT2,
			  const int *dims, int ndim, int opcode,
			  int *offs_buf, void *vals_buf, int *ovflow)
{
	SEXP ans, ans_elt, subSVT1, subSVT2;
	int ans_len, is_empty, i;

	if (SVT1 == R_NilValue || SVT2 == R_NilValue) {
		if (SVT1 == R_NilValue && SVT2 == R_NilValue)
			return R_NilValue;
		if (opcode == ADD_OPCODE || opcode == SUB_OPCODE)
			return SVT1 == R_NilValue ? SVT2 : SVT1;
	}
	if (ndim == 1) {
		/* 'SVT1' and 'SVT2' are "leaf vectors", but one of them
		   can be NULL if 'opcode' is MULT_OPCODE. */
		return _Arith_leaf_vectors(SVT1, SVT2, opcode,
					   offs_buf, vals_buf, ovflow);
	}

	/* Each of 'SVT1' and 'SVT2' is either a list or a NULL, but they
	   cannot both be NULL. */
	ans_len = dims[ndim - 1];
	ans = PROTECT(NEW_LIST(ans_len));
	subSVT1 = subSVT2 = R_NilValue;
	is_empty = 1;
	for (i = 0; i < ans_len; i++) {
		if (SVT1 != R_NilValue)
			subSVT1 = VECTOR_ELT(SVT1, i);
		if (SVT2 != R_NilValue)
			subSVT2 = VECTOR_ELT(SVT2, i);
		ans_elt = REC_SVT_Arith(subSVT1, subSVT2,
					dims, ndim - 1, opcode,
					offs_buf, vals_buf, ovflow);
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

static void check_array_conformability(SEXP x_dim, SEXP y_dim)
{
	int ndim;

	ndim = LENGTH(x_dim);
	if (ndim != LENGTH(y_dim) ||
	    memcmp(INTEGER(x_dim), INTEGER(y_dim), sizeof(int) * ndim) != 0)
		error("non-conformable arrays");
	return;
}

/* --- .Call ENTRY POINT --- */
SEXP C_SVT_Arith(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		 SEXP y_dim, SEXP y_type, SEXP y_SVT,
		 SEXP op)
{
	SEXPTYPE x_Rtype, y_Rtype;
	int opcode, *offs_buf, ovflow;
	double *vals_buf;
	SEXP ans_SVT;

	check_array_conformability(x_dim, y_dim);
	x_Rtype = _get_Rtype_from_Rstring(x_type);
	y_Rtype = _get_Rtype_from_Rstring(y_type);
	if (x_Rtype == 0 || y_Rtype == 0)
		error("SparseArray internal error in "
		      "C_SVT_Arith():\n"
                      "    invalid 'x_type' or 'y_type' value");
	opcode = _get_Arith_opcode(op, x_Rtype, y_Rtype);
	if (opcode != ADD_OPCODE &&
	    opcode != SUB_OPCODE &&
	    opcode != MULT_OPCODE)
	{
		error("\"%s\" not implemented between SVT_SparseArray "
		      "objects", CHAR(STRING_ELT(op, 0)));
	}
	offs_buf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	/* Must be big enough to contain ints or doubles. */
	vals_buf = (double *) R_alloc(INTEGER(x_dim)[0], sizeof(double));
	ovflow = 0;
	ans_SVT = REC_SVT_Arith(x_SVT, y_SVT,
				INTEGER(x_dim), LENGTH(x_dim), opcode,
				offs_buf, vals_buf, &ovflow);
	if (ans_SVT != R_NilValue)
		PROTECT(ans_SVT);
	if (ovflow)
		warning("NAs produced by integer overflow");
	if (ans_SVT != R_NilValue)
		UNPROTECT(1);
	return ans_SVT;
}

/* --- .Call ENTRY POINT --- */
SEXP C_SVT_Compare(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		   SEXP y_dim, SEXP y_type, SEXP y_SVT,
		   SEXP op)
{
	SEXPTYPE x_Rtype, y_Rtype;
	int opcode;

	check_array_conformability(x_dim, y_dim);
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

	check_array_conformability(x_dim, y_dim);
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

