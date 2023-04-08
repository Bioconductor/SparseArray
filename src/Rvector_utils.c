/****************************************************************************
 *       Some low-level convenience utilities to deal with R vectors        *
 *                      Maybe should go to S4Vectors?                       *
 ****************************************************************************/
#include "Rvector_utils.h"

#include <string.h>  /* for memset() */


/* The 7 types of R vectors (6 types of atomic vectors + the "list" type). */
static const SEXPTYPE Rvector_types[] = {
        LGLSXP,   // "logical"
        INTSXP,   // "integer"
        REALSXP,  // "double"
        CPLXSXP,  // "complex"
        RAWSXP,   // "raw"
        STRSXP,   // "character"

        VECSXP    // "list"
};

/* Return 0 if supplied 'type' is invalid string. */
SEXPTYPE _get_Rtype_from_Rstring(SEXP type)
{
	SEXP type0;
	SEXPTYPE Rtype;
	int ntypes, i;

	if (!IS_CHARACTER(type) || LENGTH(type) != 1)
		return 0;
	type0 = STRING_ELT(type, 0);
	if (type0 == NA_STRING)
		return 0;
	Rtype = str2type(CHAR(type0));
	ntypes = sizeof(Rvector_types) / sizeof(SEXPTYPE);
	for (i = 0; i < ntypes; i++)
		if (Rtype == Rvector_types[i])
			return Rtype;
	return 0;
}

size_t _get_Rtype_size(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return sizeof(int);
	    case REALSXP:             return sizeof(double);
	    case CPLXSXP:             return sizeof(Rcomplex);
	    case RAWSXP:              return sizeof(Rbyte);
	}
	return 0;
}

/* Like allocVector() but with initialization of the vector elements. */
SEXP _new_Rvector0(SEXPTYPE Rtype, R_xlen_t len)
{
	SEXP ans;
	size_t Rtype_size;

	ans = PROTECT(allocVector(Rtype, len));
	/* allocVector() does NOT initialize the vector elements, except
	   for a list or a character vector. */
	if (Rtype != STRSXP && Rtype != VECSXP) {
		Rtype_size = _get_Rtype_size(Rtype);
		if (Rtype_size == 0) {
			UNPROTECT(1);
			error("SparseArray internal error in _new_Rvector0():\n"
			      "    type \"%s\" is not supported",
			      type2char(Rtype));
		}
		memset(DATAPTR(ans), 0, Rtype_size * XLENGTH(ans));
	}
	UNPROTECT(1);
	return ans;
}

/* Like allocMatrix() but with initialization of the matrix elements and
   addition of the dimnames. */
SEXP _new_Rmatrix0(SEXPTYPE Rtype, int nrow, int ncol, SEXP dimnames)
{
	SEXP ans;
	size_t Rtype_size;

	ans = PROTECT(allocMatrix(Rtype, nrow, ncol));
	SET_DIMNAMES(ans, dimnames);
	/* allocMatrix() is just a thin wrapper around allocVector() and
	   the latter does NOT initialize the vector elements, except for
	   a list or a character vector. */
	if (Rtype != STRSXP && Rtype != VECSXP) {
		Rtype_size = _get_Rtype_size(Rtype);
		if (Rtype_size == 0) {
			UNPROTECT(1);
			error("SparseArray internal error in _new_Rmatrix0():\n"
			      "    type \"%s\" is not supported",
			      type2char(Rtype));
		}
		memset(DATAPTR(ans), 0, Rtype_size * XLENGTH(ans));
	}
	UNPROTECT(1);
	return ans;
}

/* Like allocArray() but with initialization of the array elements and
   addition of the dimnames. */
SEXP _new_Rarray0(SEXPTYPE Rtype, SEXP dim, SEXP dimnames)
{
	SEXP ans;
	size_t Rtype_size;

	ans = PROTECT(allocArray(Rtype, dim));
	SET_DIMNAMES(ans, dimnames);
	/* allocArray() is just a thin wrapper around allocVector() and
	   the latter does NOT initialize the vector elements, except for
	   a list or a character vector. */
	if (Rtype != STRSXP && Rtype != VECSXP) {
		Rtype_size = _get_Rtype_size(Rtype);
		if (Rtype_size == 0) {
			UNPROTECT(1);
			error("SparseArray internal error in _new_Rarray0():\n"
			      "    type \"%s\" is not supported",
			      type2char(Rtype));
		}
		memset(DATAPTR(ans), 0, Rtype_size * XLENGTH(ans));
	}
	UNPROTECT(1);
	return ans;
}

CopyRVectorElt_FUNType _select_copy_Rvector_elt_FUN(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return _copy_INTEGER_elt;
	    case REALSXP:             return _copy_NUMERIC_elt;
	    case CPLXSXP:             return _copy_COMPLEX_elt;
	    case RAWSXP:              return _copy_RAW_elt;
	    case STRSXP:              return _copy_CHARACTER_elt;
	    case VECSXP:              return _copy_LIST_elt;
	}
	return NULL;
}

CopyRVectorElts_FUNType _select_copy_Rvector_elts_FUN(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return _copy_INTEGER_elts;
	    case REALSXP:             return _copy_NUMERIC_elts;
	    case CPLXSXP:             return _copy_COMPLEX_elts;
	    case RAWSXP:              return _copy_RAW_elts;
	    case STRSXP:              return _copy_CHARACTER_elts;
	    case VECSXP:              return _copy_LIST_elts;
	}
	return NULL;
}


/****************************************************************************
 * _collect_offsets_of_nonzero_Rsubvec_elts()
 */

static int collect_offsets_of_nonzero_int_elts(
		const int *in, int in_len, int *offs_buf)
{
	int *off_p = offs_buf;
	for (int offset = 0; offset < in_len; offset++, in++)
		if (*in != 0)
			*(off_p++) = offset;
	return (int) (off_p - offs_buf);
}

static int collect_offsets_of_nonzero_double_elts(
		const double *in, int in_len, int *offs_buf)
{
	int *off_p = offs_buf;
	for (int offset = 0; offset < in_len; offset++, in++)
		if (*in != 0.0)
			*(off_p++) = offset;
	return (int) (off_p - offs_buf);
}

#define	IS_NONZERO_RCOMPLEX(x) ((x)->r != 0.0 || (x)->i != 0.0)
static int collect_offsets_of_nonzero_Rcomplex_elts(
		const Rcomplex *in, int in_len, int *offs_buf)
{
	int *off_p = offs_buf;
	for (int offset = 0; offset < in_len; offset++, in++)
		if (IS_NONZERO_RCOMPLEX(in))
			*(off_p++) = offset;
	return (int) (off_p - offs_buf);
}

static int collect_offsets_of_nonzero_Rbyte_elts(
		const Rbyte *in, int in_len, int *offs_buf)
{
	int *off_p = offs_buf;
	for (int offset = 0; offset < in_len; offset++, in++)
		if (*in != 0)
			*(off_p++) = offset;
	return (int) (off_p - offs_buf);
}

#define	IS_NONEMPTY_CHARSXP(x) ((x) == NA_STRING || XLENGTH(x) != 0)
static int collect_offsets_of_nonempty_character_elts(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *offs_buf)
{
	int *off_p = offs_buf;
	for (int offset = 0; offset < subvec_len; offset++, subvec_offset++) {
		SEXP Rvector_elt = STRING_ELT(Rvector, subvec_offset);
		if (IS_NONEMPTY_CHARSXP(Rvector_elt))
			*(off_p++) = offset;
	}
	return (int) (off_p - offs_buf);
}

static int collect_offsets_of_nonnull_list_elts(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *offs_buf)
{
	int *off_p = offs_buf;
	for (int offset = 0; offset < subvec_len; offset++, subvec_offset++) {
		SEXP Rvector_elt = VECTOR_ELT(Rvector, subvec_offset);
		if (Rvector_elt != R_NilValue)
			*(off_p++) = offset;
	}
	return (int) (off_p - offs_buf);
}

/* Only looks at the subvector of 'Rvector' made of the range of elements
   defined by 'subvec_offset' and 'subvec_len'.
   Offsets of nonzero elements are collected with respect to this subvector.
   Caller must make sure to supply an 'offs_buf' buffer that is big enough
   to store all the collected offsets. Safe choice is to make the buffer of
   length 'subvec_len'.
   Note that even though 'Rvector' can be a long vector, the subvector
   defined by 'subvec_offset/subvec_len' cannot i.e. 'subvec_len' must be
   supplied as an int.
   Returns the number of collected offsets. */
int _collect_offsets_of_nonzero_Rsubvec_elts(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *offs_buf)
{
	SEXPTYPE Rtype;

	Rtype = TYPEOF(Rvector);
	switch (Rtype) {
	    case LGLSXP: case INTSXP:
		return collect_offsets_of_nonzero_int_elts(
				INTEGER(Rvector) + subvec_offset,
				subvec_len, offs_buf);
	    case REALSXP:
		return collect_offsets_of_nonzero_double_elts(
				REAL(Rvector) + subvec_offset,
				subvec_len, offs_buf);
	    case CPLXSXP:
		return collect_offsets_of_nonzero_Rcomplex_elts(
				COMPLEX(Rvector) + subvec_offset,
				subvec_len, offs_buf);
	    case RAWSXP:
		return collect_offsets_of_nonzero_Rbyte_elts(
				RAW(Rvector) + subvec_offset,
				subvec_len, offs_buf);
	    case STRSXP:
		return collect_offsets_of_nonempty_character_elts(
				Rvector, subvec_offset,
				subvec_len, offs_buf);
	    case VECSXP:
		return collect_offsets_of_nonnull_list_elts(
				Rvector, subvec_offset,
				subvec_len, offs_buf);
	}
	error("SparseArray internal error in "
	      "_collect_offsets_of_nonzero_Rsubvec_elts():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return -1;  /* will never reach this */
}


/****************************************************************************
 * _reset_selected_Rvector_elts()
 */

static void reset_selected_int_elts(int *x,
		const int *selection, int n)
{
	for (int i = 0; i < n; i++, selection++)
		x[*selection] = 0;
	return;
}

static void reset_selected_double_elts(double *x,
		const int *selection, int n)
{
	for (int i = 0; i < n; i++, selection++)
		x[*selection] = 0.0;
	return;
}

static void reset_selected_Rcomplex_elts(Rcomplex *x,
		const int *selection, int n)
{
	Rcomplex *x_elt;

	for (int i = 0; i < n; i++, selection++) {
		x_elt = x + *selection;
		x_elt->r = x_elt->i = 0.0;
	}
	return;
}

static void reset_selected_Rbyte_elts(Rbyte *x,
		const int *selection, int n)
{
	for (int i = 0; i < n; i++, selection++)
		x[*selection] = 0.0;
	return;
}

static void reset_selected_character_elts(SEXP x,
		const int *selection, int n)
{
	SEXP x0;

	x0 = PROTECT(mkChar(""));
	for (int i = 0; i < n; i++, selection++)
		SET_STRING_ELT(x, *selection, x0);
	UNPROTECT(1);
	return;
}

static void reset_selected_list_elts(SEXP x,
		const int *selection, int n)
{
	for (int i = 0; i < n; i++, selection++)
		SET_VECTOR_ELT(x, *selection, R_NilValue);
	return;
}

void _reset_selected_Rvector_elts(SEXP Rvector, const int *selection, int n)
{
	SEXPTYPE Rtype;

	Rtype = TYPEOF(Rvector);
	switch (Rtype) {
	    case LGLSXP: case INTSXP:
		return reset_selected_int_elts(INTEGER(Rvector),
				selection, n);
	    case REALSXP:
		return reset_selected_double_elts(REAL(Rvector),
				selection, n);
	    case CPLXSXP:
		return reset_selected_Rcomplex_elts(COMPLEX(Rvector),
				selection, n);
	    case RAWSXP:
		return reset_selected_Rbyte_elts(RAW(Rvector),
				selection, n);
	    case STRSXP:
		return reset_selected_character_elts(Rvector,
				selection, n);
	    case VECSXP:
		return reset_selected_list_elts(Rvector,
				selection, n);
	}
	error("SparseArray internal error in "
	      "_reset_selected_Rvector_elts():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return;
}


/****************************************************************************
 * _copy_selected_Rsubvec_elts()
 */

void _copy_selected_ints(const int *in,
		const int *selection, int n, int *out)
{
	for (int k = 0; k < n; k++, selection++, out++)
		*out = in[*selection];
	return;
}

void _copy_selected_doubles(const double *in,
		const int *selection, int n, double *out)
{
	for (int k = 0; k < n; k++, selection++, out++)
		*out = in[*selection];
	return;
}

void _copy_selected_Rcomplexes(const Rcomplex *in,
		const int *selection, int n, Rcomplex *out)
{
	for (int k = 0; k < n; k++, selection++, out++)
		*out = in[*selection];
	return;
}

void _copy_selected_Rbytes(const Rbyte *in,
		const int *selection, int n, Rbyte *out)
{
	for (int k = 0; k < n; k++, selection++, out++)
		*out = in[*selection];
	return;
}

/* The selection is assumed to have the same length as 'out_Rvector'.
   Only for a 'selection' and 'out_Rvector' of length <= INT_MAX.
   Do NOT use on a 'selection' or 'out_Rvector' of length > INT_MAX! */
void _copy_selected_Rsubvec_elts(
		SEXP in_Rvector, R_xlen_t in_offset,
		const int *selection, SEXP out_Rvector)
{
	SEXPTYPE Rtype;
	int out_len;
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;

	Rtype = TYPEOF(in_Rvector);
	out_len = LENGTH(out_Rvector);  /* also the length of the selection */

	/* Optimized for LGLSXP, INTSXP, REALSXP, CPLXSXP, and RAWSXP. */
	switch (Rtype) {
	    case LGLSXP: case INTSXP:
		_copy_selected_ints(INTEGER(in_Rvector) + in_offset,
				selection, out_len, INTEGER(out_Rvector));
		return;
	    case REALSXP:
		_copy_selected_doubles(REAL(in_Rvector) + in_offset,
				selection, out_len, REAL(out_Rvector));
		return;
	    case CPLXSXP:
		_copy_selected_Rcomplexes(COMPLEX(in_Rvector) + in_offset,
				selection, out_len, COMPLEX(out_Rvector));
		return;
	    case RAWSXP:
		_copy_selected_Rbytes(RAW(in_Rvector) + in_offset,
				selection, out_len, RAW(out_Rvector));
		return;
	}

	/* STRSXP and VECSXP cases. */
	copy_Rvector_elt_FUN = _select_copy_Rvector_elt_FUN(Rtype);
	if (copy_Rvector_elt_FUN == NULL)
		error("SparseArray internal error in "
		      "_copy_selected_Rsubvec_elts():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));

	for (R_xlen_t k = 0; k < out_len; k++, selection++) {
		R_xlen_t offset = in_offset + *selection;
		copy_Rvector_elt_FUN(in_Rvector, offset, out_Rvector, k);
	}
	return;
}


/****************************************************************************
 * _copy_Rvector_elts_to_offsets()
 */

void _copy_ints_to_offsets(const int *in,
		const int *selection, int n, int *out)
{
	for (int k = 0; k < n; k++, selection++, in++)
		out[*selection] = *in;
	return;
}

void _copy_doubles_to_offsets(const double *in,
		const int *selection, int n, double *out)
{
	for (int k = 0; k < n; k++, selection++, in++)
		out[*selection] = *in;
	return;
}

void _copy_Rcomplexes_to_offsets(const Rcomplex *in,
		const int *selection, int n, Rcomplex *out)
{
	for (int k = 0; k < n; k++, selection++, in++)
		out[*selection] = *in;
	return;
}

void _copy_Rbytes_to_offsets(const Rbyte *in,
		const int *selection, int n, Rbyte *out)
{
	for (int k = 0; k < n; k++, selection++, in++)
		out[*selection] = *in;
	return;
}

/* The selection is assumed to have the same length as 'in_Rvector'.
   Only for a 'selection' and 'in_Rvector' of length <= INT_MAX.
   Do NOT use on a 'selection' or 'in_Rvector' of length > INT_MAX! */
void _copy_Rvector_elts_to_offsets(
		SEXP in_Rvector, const int *selection,
		SEXP out_Rvector, R_xlen_t out_offset)
{
	SEXPTYPE Rtype;
	int in_len;
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;

	Rtype = TYPEOF(in_Rvector);
	in_len = LENGTH(in_Rvector);  /* also the length of the selection */

	/* Optimized for LGLSXP, INTSXP, REALSXP, CPLXSXP, and RAWSXP. */
	switch (Rtype) {
	    case LGLSXP: case INTSXP:
		_copy_ints_to_offsets(INTEGER(in_Rvector),
				selection, in_len,
				INTEGER(out_Rvector) + out_offset);
		return;
	    case REALSXP:
		_copy_doubles_to_offsets(REAL(in_Rvector),
				selection, in_len,
				REAL(out_Rvector) + out_offset);
		return;
	    case CPLXSXP:
		_copy_Rcomplexes_to_offsets(COMPLEX(in_Rvector),
				selection, in_len,
				COMPLEX(out_Rvector) + out_offset);
		return;
	    case RAWSXP:
		_copy_Rbytes_to_offsets(RAW(in_Rvector),
				selection, in_len,
				RAW(out_Rvector) + out_offset);
		return;
	}

	/* STRSXP and VECSXP cases. */
	copy_Rvector_elt_FUN = _select_copy_Rvector_elt_FUN(Rtype);
	if (copy_Rvector_elt_FUN == NULL)
		error("SparseArray internal error in "
		      "_copy_Rvector_elts_to_offsets():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));

	for (R_xlen_t k = 0; k < in_len; k++, selection++) {
		R_xlen_t offset = out_offset + *selection;
		copy_Rvector_elt_FUN(in_Rvector, k, out_Rvector, offset);
	}
	return;
}


/****************************************************************************
 * _copy_Rvector_elts_from_selected_offsets()
 * _copy_Rvector_elts_from_selected_lloffsets()
 */

static void copy_ints_from_selected_offsets(const int *in,
		const int *offsets, const int *offset_selection,
		int n, int *out)
{
	for (int k = 0; k < n; k++, offset_selection++, out++)
		*out = in[offsets[*offset_selection]];
	return;
}
static void copy_ints_from_selected_lloffsets(const int *in,
		const long long *lloffsets, const int *lloffset_selection,
		int n, int *out)
{
	for (int k = 0; k < n; k++, lloffset_selection++, out++)
		*out = in[lloffsets[*lloffset_selection]];
	return;
}

static void copy_doubles_from_selected_offsets(const double *in,
		const int *offsets, const int *offset_selection,
		int n, double *out)
{
	for (int k = 0; k < n; k++, offset_selection++, out++)
		*out = in[offsets[*offset_selection]];
	return;
}
static void copy_doubles_from_selected_lloffsets(const double *in,
		const long long *lloffsets, const int *lloffset_selection,
		int n, double *out)
{
	for (int k = 0; k < n; k++, lloffset_selection++, out++)
		*out = in[lloffsets[*lloffset_selection]];
	return;
}

static void copy_Rcomplexes_from_selected_offsets(const Rcomplex *in,
		const int *offsets, const int *offset_selection,
		int n, Rcomplex *out)
{
	for (int k = 0; k < n; k++, offset_selection++, out++)
		*out = in[offsets[*offset_selection]];
	return;
}
static void copy_Rcomplexes_from_selected_lloffsets(const Rcomplex *in,
		const long long *lloffsets, const int *lloffset_selection,
		int n, Rcomplex *out)
{
	for (int k = 0; k < n; k++, lloffset_selection++, out++)
		*out = in[lloffsets[*lloffset_selection]];
	return;
}

static void copy_Rbytes_from_selected_offsets(const Rbyte *in,
		const int *offsets, const int *offset_selection,
		int n, Rbyte *out)
{
	for (int k = 0; k < n; k++, offset_selection++, out++)
		*out = in[offsets[*offset_selection]];
	return;
}
static void copy_Rbytes_from_selected_lloffsets(const Rbyte *in,
		const long long *lloffsets, const int *lloffset_selection,
		int n, Rbyte *out)
{
	for (int k = 0; k < n; k++, lloffset_selection++, out++)
		*out = in[lloffsets[*lloffset_selection]];
	return;
}

/* The selection is assumed to have the same length as 'out_Rvector'.
   Only for an 'offset_selection' and 'out_Rvector' of length <= INT_MAX.
   Don't use on an 'offset_selection' or 'out_Rvector' of length > INT_MAX! */
void _copy_Rvector_elts_from_selected_offsets(SEXP in_Rvector,
		const int *offsets, const int *offset_selection,
		SEXP out_Rvector)
{
	SEXPTYPE Rtype;
	int out_len;
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;

	Rtype = TYPEOF(in_Rvector);
	out_len = LENGTH(out_Rvector);  /* also the length of the selection */

	/* Optimized for LGLSXP, INTSXP, REALSXP, CPLXSXP, and RAWSXP. */
	switch (Rtype) {
	    case LGLSXP: case INTSXP:
		copy_ints_from_selected_offsets(INTEGER(in_Rvector),
				offsets, offset_selection, out_len,
				INTEGER(out_Rvector));
		return;
	    case REALSXP:
		copy_doubles_from_selected_offsets(REAL(in_Rvector),
				offsets, offset_selection, out_len,
				REAL(out_Rvector));
		return;
	    case CPLXSXP:
		copy_Rcomplexes_from_selected_offsets(COMPLEX(in_Rvector),
				offsets, offset_selection, out_len,
				COMPLEX(out_Rvector));
		return;
	    case RAWSXP:
		copy_Rbytes_from_selected_offsets(RAW(in_Rvector),
				offsets, offset_selection, out_len,
				RAW(out_Rvector));
		return;
	}

	/* STRSXP and VECSXP cases. */
	copy_Rvector_elt_FUN = _select_copy_Rvector_elt_FUN(Rtype);
	if (copy_Rvector_elt_FUN == NULL)
		error("SparseArray internal error in "
		      "_copy_Rvector_elts_from_selected_offsets():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));

	for (R_xlen_t k = 0; k < out_len; k++, offset_selection++)
		copy_Rvector_elt_FUN(
			in_Rvector, (R_xlen_t) offsets[*offset_selection],
			out_Rvector, k);
	return;
}

/* The selection is assumed to have the same length as 'out_Rvector'.
   Only for an 'lloffset_selection' and 'out_Rvector' of length <= INT_MAX.
   Don't use on an 'lloffset_selection' or 'out_Rvector' of length > INT_MAX! */
void _copy_Rvector_elts_from_selected_lloffsets(SEXP in_Rvector,
		const long long *lloffsets, const int *lloffset_selection,
		SEXP out_Rvector)
{
	SEXPTYPE Rtype;
	int out_len;
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;

	Rtype = TYPEOF(in_Rvector);
	out_len = LENGTH(out_Rvector);  /* also the length of the selection */

	/* Optimized for LGLSXP, INTSXP, REALSXP, CPLXSXP, and RAWSXP. */
	switch (Rtype) {
	    case LGLSXP: case INTSXP:
		copy_ints_from_selected_lloffsets(INTEGER(in_Rvector),
				lloffsets, lloffset_selection, out_len,
				INTEGER(out_Rvector));
		return;
	    case REALSXP:
		copy_doubles_from_selected_lloffsets(REAL(in_Rvector),
				lloffsets, lloffset_selection, out_len,
				REAL(out_Rvector));
		return;
	    case CPLXSXP:
		copy_Rcomplexes_from_selected_lloffsets(COMPLEX(in_Rvector),
				lloffsets, lloffset_selection, out_len,
				COMPLEX(out_Rvector));
		return;
	    case RAWSXP:
		copy_Rbytes_from_selected_lloffsets(RAW(in_Rvector),
				lloffsets, lloffset_selection, out_len,
				RAW(out_Rvector));
		return;
	}

	/* STRSXP and VECSXP cases. */
	copy_Rvector_elt_FUN = _select_copy_Rvector_elt_FUN(Rtype);
	if (copy_Rvector_elt_FUN == NULL)
		error("SparseArray internal error in "
		      "_copy_Rvector_elts_from_selected_lloffsets():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));

	for (R_xlen_t k = 0; k < out_len; k++, lloffset_selection++)
		copy_Rvector_elt_FUN(
			in_Rvector, (R_xlen_t) lloffsets[*lloffset_selection],
			out_Rvector, k);
	return;
}

