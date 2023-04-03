/****************************************************************************
 *        crossprod(), tcrossprod(), and %*% of SparseArray objects         *
 ****************************************************************************/
#include "SparseMatrix_mult.h"

#include "Rvector_utils.h"
#include "leaf_vector_utils.h"

#include <string.h>  /* for memset() */


/* --- .Call ENTRY POINT --- */
SEXP C_SVT_crossprod1(SEXP x_dim, SEXP x_SVT, SEXP ans_type, SEXP ans_dimnames)
{
	SEXPTYPE ans_Rtype;
	size_t ans_Rtype_size;
	int x_ncol, i, j;
	SEXP ans, subSVT1, subSVT2;
	double dp;

	/* Check 'ans_type'. */
	ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (ans_Rtype == 0)
		error("SparseArray internal error in "
		      "C_SVT_crossprod1():\n"
		      "    invalid 'ans_type' value");
	if (ans_Rtype != REALSXP)
		error("SparseArray internal error in "
		      "C_SVT_crossprod1():\n"
		      "    type \"%s\" is not supported yet",
		      type2char(ans_Rtype));
	ans_Rtype_size = _get_Rtype_size(ans_Rtype);

	/* Check 'x_dim'. */
	if (LENGTH(x_dim) != 2)
		error("'x' must have 2 dimensions");
	x_ncol = INTEGER(x_dim)[1];

	/* Allocate 'ans' and fill with zeros. */
	ans = PROTECT(allocMatrix(ans_Rtype, x_ncol, x_ncol));
	SET_DIMNAMES(ans, ans_dimnames);
        /* allocMatrix() is just a thin wrapper around allocVector() and
           the latter does NOT initialize the vector elements, except for
           a list or a character vector. */
	if (ans_Rtype != STRSXP && ans_Rtype != VECSXP)
		memset(DATAPTR(ans), 0, ans_Rtype_size * XLENGTH(ans));

	/* The easy case. */
	if (x_SVT == R_NilValue) {
		UNPROTECT(1);
		return(ans);
	}

	/* The general case. */
	for (j = 0; j < x_ncol; j++) {
		subSVT2 = VECTOR_ELT(x_SVT, j);
		if (subSVT2 != R_NilValue) {
			dp = _dotprod_leaf_vectors(subSVT2, subSVT2);
			REAL(ans)[j + j * x_ncol] = dp;
		}
		for (i = j + 1; i < x_ncol; i++) {
			subSVT1 = VECTOR_ELT(x_SVT, i);
			if (subSVT1 == R_NilValue && subSVT2 == R_NilValue)
				continue;
			if (subSVT1 == R_NilValue) {
				dp = _dotprod0_leaf_vector(subSVT2);
			} else if (subSVT2 == R_NilValue) {
				dp = _dotprod0_leaf_vector(subSVT1);
			} else {
				dp = _dotprod_leaf_vectors(subSVT1, subSVT2);
			}
			REAL(ans)[i + j * x_ncol] = dp;
			REAL(ans)[j + i * x_ncol] = dp;
		}
	}

	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_SVT_crossprod2(SEXP x_dim, SEXP x_SVT, SEXP y_dim, SEXP y_SVT,
		      SEXP ans_type, SEXP ans_dimnames)
{
	SEXPTYPE ans_Rtype;
	size_t ans_Rtype_size;
	int x_ncol, y_ncol, i, j;
	SEXP ans, subSVT1, subSVT2;
	double dp;

	/* Check 'ans_type'. */
	ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (ans_Rtype == 0)
		error("SparseArray internal error in "
		      "C_SVT_crossprod2():\n"
		      "    invalid 'ans_type' value");
	if (ans_Rtype != REALSXP)
		error("SparseArray internal error in "
		      "C_SVT_crossprod2():\n"
		      "    type \"%s\" is not supported yet",
		      type2char(ans_Rtype));
	ans_Rtype_size = _get_Rtype_size(ans_Rtype);

	/* Check 'x_dim' and 'y_dim'. */
	if (LENGTH(x_dim) != 2 || LENGTH(y_dim) != 2)
		error("arguments must have 2 dimensions");
	if (INTEGER(x_dim)[0] != INTEGER(y_dim)[0])
		error("non-conformable arguments");
	x_ncol = INTEGER(x_dim)[1];
	y_ncol = INTEGER(y_dim)[1];

	/* Allocate 'ans' and fill with zeros. */
	ans = PROTECT(allocMatrix(ans_Rtype, x_ncol, y_ncol));
	SET_DIMNAMES(ans, ans_dimnames);
        /* allocMatrix() is just a thin wrapper around allocVector() and
           the latter does NOT initialize the vector elements, except for
           a list or a character vector. */
	if (ans_Rtype != STRSXP && ans_Rtype != VECSXP)
		memset(DATAPTR(ans), 0, ans_Rtype_size * XLENGTH(ans));

	/* The easy case. */
	if (x_SVT == R_NilValue && y_SVT == R_NilValue) {
		UNPROTECT(1);
		return(ans);
	}

	/* The general case. */
	subSVT1 = subSVT2 = R_NilValue;
	for (j = 0; j < y_ncol; j++) {
		if (y_SVT != R_NilValue)
			subSVT2 = VECTOR_ELT(y_SVT, j);
		for (i = 0; i < x_ncol; i++) {
			if (x_SVT != R_NilValue)
				subSVT1 = VECTOR_ELT(x_SVT, i);
			if (subSVT1 == R_NilValue && subSVT2 == R_NilValue)
				continue;
			if (subSVT1 == R_NilValue) {
				dp = _dotprod0_leaf_vector(subSVT2);
			} else if (subSVT2 == R_NilValue) {
				dp = _dotprod0_leaf_vector(subSVT1);
			} else {
				dp = _dotprod_leaf_vectors(subSVT1, subSVT2);
			}
			REAL(ans)[i + j * x_ncol] = dp;
		}
	}

	UNPROTECT(1);
	return ans;
}

