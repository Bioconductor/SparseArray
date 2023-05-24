/****************************************************************************
 *                  Transposition of a SparseArray object                   *
 ****************************************************************************/
#include "SparseArray_aperm.h"

#include "Rvector_utils.h"
#include "leaf_vector_utils.h"
#include "SBT_utils.h"


static void push_leaf_vector_to_SBT_row(
		SEXP SBT, const int *dim, int ndim,
		SEXP lv, int *coords0_buf)
{
	int lv_len, k;
	SEXP lv_offs, lv_vals;
	const int *lv_offs_p;

	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	lv_offs_p = INTEGER(lv_offs);
	switch (TYPEOF(lv_vals)) {
	    case LGLSXP: case INTSXP: {
		const int *lv_vals_p = INTEGER(lv_vals);
		for (k = 0; k < lv_len; k++) {
			coords0_buf[ndim - 1] = *lv_offs_p;
			_push_int_to_SBT(SBT, dim, ndim,
					 coords0_buf, *lv_vals_p);
			lv_offs_p++;
			lv_vals_p++;
		}
		return;
	    }
	    case REALSXP: {
		const double *lv_vals_p = REAL(lv_vals);
		for (k = 0; k < lv_len; k++) {
			coords0_buf[ndim - 1] = *lv_offs_p;
			_push_double_to_SBT(SBT, dim, ndim,
					    coords0_buf, *lv_vals_p);
			lv_offs_p++;
			lv_vals_p++;
		}
		return;
	    }
	    case CPLXSXP: {
		const Rcomplex *lv_vals_p = COMPLEX(lv_vals);
		for (k = 0; k < lv_len; k++) {
			coords0_buf[ndim - 1] = *lv_offs_p;
			_push_Rcomplex_to_SBT(SBT, dim, ndim,
					      coords0_buf, *lv_vals_p);
			lv_offs_p++;
			lv_vals_p++;
		}
		return;
	    }
	    case RAWSXP: {
		const Rbyte *lv_vals_p = RAW(lv_vals);
		for (k = 0; k < lv_len; k++) {
			coords0_buf[ndim - 1] = *lv_offs_p;
			_push_Rbyte_to_SBT(SBT, dim, ndim,
					   coords0_buf, *lv_vals_p);
			lv_offs_p++;
			lv_vals_p++;
		}
		return;
	    }
	    case STRSXP: {
		for (k = 0; k < lv_len; k++) {
			coords0_buf[ndim - 1] = *lv_offs_p;
			_push_SEXP_to_SBT(SBT, dim, ndim,
					  coords0_buf, STRING_ELT(lv_vals, k));
			lv_offs_p++;
		}
		return;
	    }
	    case VECSXP: {
		for (k = 0; k < lv_len; k++) {
			coords0_buf[ndim - 1] = *lv_offs_p;
			_push_SEXP_to_SBT(SBT, dim, ndim, coords0_buf,
						VECTOR_ELT(lv_vals, k));
			lv_offs_p++;
		}
		return;
	    }
	}
	error("SparseArray internal error in "
	      "push_leaf_vector_to_SBT_row():\n"
	      "    type \"%s\" is not supported", type2char(TYPEOF(lv_vals)));
	return;  /* will never reach this */
}

/* Recursive. */
static void REC_push_transposed_SVT_to_SBT(
		SEXP SBT, const int *SBT_dim, int SBT_ndim,
		SEXP SVT, const int *SVT_dim, int SVT_ndim,
		int *coords0_buf)
{
	int SVT_len, i;
	SEXP subSVT;

	SVT_len = LENGTH(SVT);
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		if (subSVT == R_NilValue)
			continue;
		coords0_buf[SBT_ndim - SVT_ndim] = i;
		if (SVT_ndim == 2) {
			push_leaf_vector_to_SBT_row(
					SBT, SBT_dim, SBT_ndim,
					subSVT, coords0_buf);
		} else {
			REC_push_transposed_SVT_to_SBT(
					SBT, SBT_dim, SBT_ndim,
					subSVT, SVT_dim, SVT_ndim - 1,
					coords0_buf);
		}
	}
	return;
}

/* --- .Call ENTRY POINT --- */
SEXP C_transpose_SVT(SEXP x_dim, SEXP x_type, SEXP x_SVT)
{
	SEXPTYPE Rtype;
	int ndim, *ans_dim, along, *coords0_buf;
	SEXP ans;

	Rtype = _get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("SparseArray internal error in "
		      "C_transpose_SVT():\n"
		      "    SVT_SparseArray object has invalid type");

	ndim = LENGTH(x_dim);
	if (ndim < 2)
		error("object to transpose must have at least 2 dimensions");

	if (x_SVT == R_NilValue)
		return R_NilValue;

	ans_dim = (int *) R_alloc(ndim, sizeof(int));
	for (along = 0; along < ndim; along++)
		ans_dim[along] = INTEGER(x_dim)[ndim - 1 - along];

	coords0_buf = (int *) R_alloc(ndim, sizeof(int));

	ans = PROTECT(NEW_LIST(ans_dim[ndim - 1]));
	REC_push_transposed_SVT_to_SBT(ans, ans_dim, ndim,
				       x_SVT, INTEGER(x_dim), ndim,
				       coords0_buf);
	_SBT2SVT(ans, ans_dim, ndim, Rtype);
	UNPROTECT(1);
	return ans;
}

