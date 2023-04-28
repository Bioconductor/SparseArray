/****************************************************************************
 *                  'Ops' methods for SparseArray objects                   *
 ****************************************************************************/
#include "SparseArray_Ops_methods.h"

#include "Rvector_utils.h"
#include "leaf_vector_Arith.h"
#include "leaf_vector_Compare.h"
#include "leaf_vector_Logic.h"
#include "SVT_SparseArray_class.h"  /* for _coerce_SVT() */

#include <string.h>  /* for memcmp() */


static void REC_unary_minus_SVT(SEXP SVT, const int *dims, int ndim)
{
	int SVT_len, i;

	if (SVT == R_NilValue)
		return;
	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		_unary_minus_leaf_vector(SVT, 0);
		return;
	}
	SVT_len = dims[ndim - 1];
	for (i = 0; i < SVT_len; i++)
		REC_unary_minus_SVT(VECTOR_ELT(SVT, i), dims, ndim - 1);
	return;
}

static SEXP REC_Arith_SVT1_v2(SEXP SVT1, SEXP v2,
			      const int *dims, int ndim,
			      int opcode, SEXPTYPE ans_Rtype,
			      int *offs_buf, void *vals_buf, int *ovflow)
{
	int ans_len, is_empty, i;
	SEXP ans, ans_elt, subSVT1;

	if (SVT1 == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* 'SVT1' is a "leaf vector". */
		return _Arith_lv1_v2(SVT1, v2, opcode, ans_Rtype,
				     offs_buf, vals_buf, ovflow);
	}

	/* 'SVT1' is a list. */
	ans_len = dims[ndim - 1];
	ans = PROTECT(NEW_LIST(ans_len));
	is_empty = 1;
	for (i = 0; i < ans_len; i++) {
		subSVT1 = VECTOR_ELT(SVT1, i);
		ans_elt = REC_Arith_SVT1_v2(subSVT1, v2,
					    dims, ndim - 1,
					    opcode, ans_Rtype,
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

static SEXP REC_Compare_SVT1_v2(SEXP SVT1, SEXP v2,
				const int *dims, int ndim,
				int opcode, int *offs_buf, void *vals_buf)
{
	int ans_len, is_empty, i;
	SEXP ans, ans_elt, subSVT1;

	if (SVT1 == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* 'SVT1' is a "leaf vector". */
		return _Compare_lv1_v2(SVT1, v2, opcode, offs_buf, vals_buf);
	}

	/* 'SVT1' is a list. */
	ans_len = dims[ndim - 1];
	ans = PROTECT(NEW_LIST(ans_len));
	is_empty = 1;
	for (i = 0; i < ans_len; i++) {
		subSVT1 = VECTOR_ELT(SVT1, i);
		ans_elt = REC_Compare_SVT1_v2(subSVT1, v2,
					      dims, ndim - 1,
					      opcode, offs_buf, vals_buf);
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

static SEXP REC_Arith_SVT1_SVT2(SEXP SVT1, SEXPTYPE Rtype1,
				SEXP SVT2, SEXPTYPE Rtype2,
				const int *dims, int ndim,
				int opcode, SEXPTYPE ans_Rtype,
				int *offs_buf, void *vals_buf, int *ovflow)
{
	int ans_len, is_empty, i;
	SEXP ans, ans_elt, subSVT1, subSVT2;

	if (SVT1 == R_NilValue) {
		if (SVT2 == R_NilValue)
			return R_NilValue;
		if (opcode == ADD_OPCODE)
			return _coerce_SVT(SVT2, dims, ndim,
					   Rtype2, ans_Rtype, offs_buf);
	} else if (SVT2 == R_NilValue) {
		if (opcode == ADD_OPCODE || opcode == SUB_OPCODE)
			return _coerce_SVT(SVT1, dims, ndim,
					   Rtype1, ans_Rtype, offs_buf);
	}

	if (ndim == 1) {
		/* 'SVT1' and 'SVT2' are "leaf vectors", but:
		   - 'SVT1' can be NULL if 'opcode' is SUB_OPCODE;
		   - either 'SVT1' or 'SVT2' (but not both) can be NULL
		     if 'opcode' is MULT_OPCODE. */
		return _Arith_lv1_lv2(SVT1, SVT2, opcode, ans_Rtype,
				      offs_buf, vals_buf, ovflow);
	}

	/* Each of 'SVT1' and 'SVT2' is either a list or NULL, but they
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
		ans_elt = REC_Arith_SVT1_SVT2(subSVT1, Rtype1, subSVT2, Rtype2,
					      dims, ndim - 1,
					      opcode, ans_Rtype,
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

static SEXP REC_Compare_SVT1_SVT2(SEXP SVT1, SEXP SVT2,
				  const int *dims, int ndim,
				  int opcode, int *offs_buf, int *vals_buf)
{
	int ans_len, is_empty, i;
	SEXP ans, ans_elt, subSVT1, subSVT2;

	if (SVT1 == R_NilValue && SVT2 == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* Each of 'SVT1' and 'SVT2' is either a "leaf vector" or NULL,
		   but they cannot both be NULL. */
		return _Compare_lv1_lv2(SVT1, SVT2, opcode, offs_buf, vals_buf);
	}

	/* Each of 'SVT1' and 'SVT2' is either a list or NULL, but they
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
		ans_elt = REC_Compare_SVT1_SVT2(subSVT1, subSVT2,
						dims, ndim - 1,
						opcode, offs_buf, vals_buf);
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

static SEXP REC_Logic_SVT1_SVT2(SEXP SVT1, SEXP SVT2,
				const int *dims, int ndim,
				int opcode, int *offs_buf, int *vals_buf)
{
	int ans_len, is_empty, i;
	SEXP ans, ans_elt, subSVT1, subSVT2;

	if (SVT1 == R_NilValue || SVT2 == R_NilValue) {
		if (opcode == AND_OPCODE)
			return R_NilValue;
		return SVT1 == R_NilValue ? SVT2 : SVT1;
	}

	if (ndim == 1)
		return _Logic_lv1_lv2(SVT1, SVT2, opcode, offs_buf, vals_buf);

	/* Each of 'SVT1' and 'SVT2' is a list. */
	ans_len = dims[ndim - 1];
	ans = PROTECT(NEW_LIST(ans_len));
	is_empty = 1;
	for (i = 0; i < ans_len; i++) {
		subSVT1 = VECTOR_ELT(SVT1, i);
		subSVT2 = VECTOR_ELT(SVT2, i);
		ans_elt = REC_Logic_SVT1_SVT2(subSVT1, subSVT2,
					      dims, ndim - 1,
					      opcode, offs_buf, vals_buf);
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
	SEXP ans;

	ans = PROTECT(duplicate(x_SVT));
	REC_unary_minus_SVT(ans, INTEGER(x_dim), LENGTH(x_dim));
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_Arith_SVT1_v2(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP v2,
		     SEXP op, SEXP ans_type)
{
	SEXPTYPE x_Rtype, ans_Rtype;
	int opcode, *offs_buf, ovflow;
	double *vals_buf;
	SEXP ans;

	x_Rtype = _get_Rtype_from_Rstring(x_type);
	ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (x_Rtype == 0 || ans_Rtype == 0)
		error("SparseArray internal error in "
		      "C_Arith_SVT1_v2():\n"
                      "    invalid 'x_type' or 'ans_type' value");
	opcode = _get_Arith_opcode(op);
	if (opcode != MULT_OPCODE &&
	    opcode != DIV_OPCODE &&
	    opcode != POW_OPCODE &&
	    opcode != MOD_OPCODE &&
	    opcode != IDIV_OPCODE)
	{
		error("\"%s\" is not supported between an SVT_SparseArray "
		      "object and a numeric vector", CHAR(STRING_ELT(op, 0)));
	}
	offs_buf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	/* Must be big enough to contain ints or doubles. */
	vals_buf = (double *) R_alloc(INTEGER(x_dim)[0], sizeof(double));
	ovflow = 0;
	ans = REC_Arith_SVT1_v2(x_SVT, v2,
				INTEGER(x_dim), LENGTH(x_dim),
				opcode, ans_Rtype,
				offs_buf, vals_buf, &ovflow);
	if (ans != R_NilValue)
		PROTECT(ans);
	if (ovflow)
		warning("NAs produced by integer overflow");
	if (ans != R_NilValue)
		UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_Compare_SVT1_v2(SEXP x_dim, SEXP x_type, SEXP x_SVT, SEXP v2, SEXP op)
{
	SEXPTYPE x_Rtype;
	int opcode, *offs_buf, *vals_buf;

	x_Rtype = _get_Rtype_from_Rstring(x_type);
	if (x_Rtype == 0)
		error("SparseArray internal error in "
		      "C_Compare_SVT1_v2():\n"
                      "    invalid 'x_type'");
	opcode = _get_Compare_opcode(op);
	offs_buf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	vals_buf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	return REC_Compare_SVT1_v2(x_SVT, v2,
				   INTEGER(x_dim), LENGTH(x_dim),
				   opcode, offs_buf, vals_buf);
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
SEXP C_Arith_SVT1_SVT2(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		       SEXP y_dim, SEXP y_type, SEXP y_SVT,
		       SEXP op, SEXP ans_type)
{
	SEXPTYPE x_Rtype, y_Rtype, ans_Rtype;
	int opcode, *offs_buf, ovflow;
	double *vals_buf;
	SEXP ans;

	check_array_conformability(x_dim, y_dim);
	x_Rtype = _get_Rtype_from_Rstring(x_type);
	y_Rtype = _get_Rtype_from_Rstring(y_type);
	ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (x_Rtype == 0 || y_Rtype == 0 || ans_Rtype == 0)
		error("SparseArray internal error in "
		      "C_Arith_SVT1_SVT2():\n"
                      "    invalid 'x_type', 'y_type', or 'ans_type' value");
	opcode = _get_Arith_opcode(op);
	if (opcode != ADD_OPCODE &&
	    opcode != SUB_OPCODE &&
	    opcode != MULT_OPCODE)
	{
		error("\"%s\" is not supported between SVT_SparseArray "
		      "objects", CHAR(STRING_ELT(op, 0)));
	}
	offs_buf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	/* Must be big enough to contain ints or doubles. */
	vals_buf = (double *) R_alloc(INTEGER(x_dim)[0], sizeof(double));
	ovflow = 0;
	ans = REC_Arith_SVT1_SVT2(x_SVT, x_Rtype, y_SVT, y_Rtype,
				  INTEGER(x_dim), LENGTH(x_dim),
				  opcode, ans_Rtype,
				  offs_buf, vals_buf, &ovflow);
	if (ans != R_NilValue)
		PROTECT(ans);
	if (ovflow)
		warning("NAs produced by integer overflow");
	if (ans != R_NilValue)
		UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT ---
   'x_type' and 'y_type' asre ignored. */
SEXP C_Compare_SVT1_SVT2(SEXP x_dim, SEXP x_type, SEXP x_SVT,
			 SEXP y_dim, SEXP y_type, SEXP y_SVT,
			 SEXP op)
{
	int opcode, *offs_buf, *vals_buf;

	check_array_conformability(x_dim, y_dim);
	opcode = _get_Compare_opcode(op);
	if (opcode != NE_OPCODE &&
	    opcode != LT_OPCODE &&
	    opcode != GT_OPCODE)
	{
		error("\"%s\" is not supported between SVT_SparseArray "
		      "objects", CHAR(STRING_ELT(op, 0)));
	}
	offs_buf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	vals_buf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	return REC_Compare_SVT1_SVT2(x_SVT, y_SVT,
				     INTEGER(x_dim), LENGTH(x_dim),
				     opcode, offs_buf, vals_buf);
}

/* --- .Call ENTRY POINT ---
   'x_type' and 'y_type' are ignored. */
SEXP C_Logic_SVT1_SVT2(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		       SEXP y_dim, SEXP y_type, SEXP y_SVT,
		       SEXP op)
{
	int opcode, *offs_buf, *vals_buf;

	check_array_conformability(x_dim, y_dim);
	opcode = _get_Logic_opcode(op);
	offs_buf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	vals_buf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	return REC_Logic_SVT1_SVT2(x_SVT, y_SVT,
				   INTEGER(x_dim), LENGTH(x_dim),
				   opcode, offs_buf, vals_buf);
}

