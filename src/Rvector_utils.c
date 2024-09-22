/****************************************************************************
 *       Some low-level convenience utilities to deal with R vectors        *
 *                      Maybe should go to S4Vectors?                       *
 ****************************************************************************/
#include "Rvector_utils.h"

#include <string.h>  /* for memset() and memcpy() */

/* Initialized in R_init_SparseArray(). */
int intNA;
double doubleNA;
Rcomplex RcomplexNA;

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

/* Returns 0 if supplied 'type' is invalid string. */
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

/* Supported types: "integer", "logical", "double", "complex", or "raw". */
size_t _get_Rtype_size(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case INTSXP: case LGLSXP: return sizeof(int);
	    case REALSXP:             return sizeof(double);
	    case CPLXSXP:             return sizeof(Rcomplex);
	    case RAWSXP:              return sizeof(Rbyte);
	}
	return 0;
}


/****************************************************************************
 * _set_elts_to_val()
 * _set_elts_to_zero()
 * _set_elts_to_one()
 * _set_elts_to_minus_one()
 * _set_elts_to_NA()
 *
 * _set_Rsubvec_elts_to_val()
 * _set_Rsubvec_elts_to_zero()
 * _set_Rsubvec_elts_to_one()
 * _set_Rsubvec_elts_to_minus_one()
 * _set_Rsubvec_elts_to_NA()
 *
 * _set_Rvector_elts_to_val()
 * _set_Rvector_elts_to_zero()
 * _set_Rvector_elts_to_one()
 * _set_Rvector_elts_to_minus_one()
 * _set_Rvector_elts_to_NA()
 */

#define	DEFINE_set_TYPE_elts_to_val_FUN(type)				\
static void set_ ## type ## _elts_to_val				\
	(type *x, R_xlen_t n, type val)					\
{									\
	for (R_xlen_t i = 0; i < n; i++)				\
		x[i] = val;						\
	return;								\
}

DEFINE_set_TYPE_elts_to_val_FUN(int)
DEFINE_set_TYPE_elts_to_val_FUN(double)
DEFINE_set_TYPE_elts_to_val_FUN(Rcomplex)
DEFINE_set_TYPE_elts_to_val_FUN(Rbyte)

#define	set_TYPE_elts_to_val(type, x, offset, n, val) \
	set_ ## type ## _elts_to_val(SHIFT_DATAPTR(type, (x), (offset)), (n), \
				     (*((type *) (val))))

/* Restricted to types "integer", "logical", "double", "complex", and "raw". */
void _set_elts_to_val(SEXPTYPE Rtype, void *x, R_xlen_t offset, R_xlen_t n,
		      const void *val)
{
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		set_TYPE_elts_to_val(int,      x, offset, n, val);
		return;
	    case REALSXP:
		set_TYPE_elts_to_val(double,   x, offset, n, val);
		return;
	    case CPLXSXP:
		set_TYPE_elts_to_val(Rcomplex, x, offset, n, val);
		return;
	    case RAWSXP:
		set_TYPE_elts_to_val(Rbyte,    x, offset, n, val);
		return;
	}
	error("SparseArray internal error in _set_elts_to_val():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return;
}

/* Restricted to types "integer", "logical", "double", "complex", and "raw". */
void _set_elts_to_zero(SEXPTYPE Rtype, void *x, R_xlen_t offset, R_xlen_t n)
{
	size_t Rtype_size = _get_Rtype_size(Rtype);
	if (Rtype_size == 0)
		error("SparseArray internal error in _set_elts_to_zero():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));
	memset(shift_dataptr(Rtype, x, offset), 0, Rtype_size * n);
	return;
}

#define	set_TYPE_elts_to_one(type, x, offset, n) \
	set_ ## type ## _elts_to_val(SHIFT_DATAPTR(type, x, offset), (n), \
				     (type ## 1))

/* Restricted to types "integer", "logical", "double", "complex", and "raw". */
void _set_elts_to_one(SEXPTYPE Rtype, void *x, R_xlen_t offset, R_xlen_t n)
{
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		set_TYPE_elts_to_one(int, x, offset, n);
		return;
	    case REALSXP:
		set_TYPE_elts_to_one(double, x, offset, n);
		return;
	    case CPLXSXP:
		set_TYPE_elts_to_one(Rcomplex, x, offset, n);
		return;
	    case RAWSXP:
		set_TYPE_elts_to_one(Rbyte, x, offset, n);
		return;
	}
	error("SparseArray internal error in _set_elts_to_one():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return;
}

/* Restricted to types "integer", "double", and "complex". */
void _set_elts_to_minus_one(SEXPTYPE Rtype, void *x,
			    R_xlen_t offset, R_xlen_t n)
{
	switch (Rtype) {
	    case INTSXP:
		set_int_elts_to_val(
			SHIFT_DATAPTR(int, x, offset), n, -int1);
		return;
	    case REALSXP:
		set_double_elts_to_val(
			SHIFT_DATAPTR(double, x, offset), n, -double1);
		return;
	    case CPLXSXP: {
		const Rcomplex minus_one = {{-double1, double0}};
		set_Rcomplex_elts_to_val(
			SHIFT_DATAPTR(Rcomplex, x, offset), n, minus_one);
		return;
	    }
	}
	error("SparseArray internal error in _set_elts_to_minus_one():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return;
}

/* Restricted to types "integer", "logical", "double", and "complex". */
void _set_elts_to_NA(SEXPTYPE Rtype, void *x, R_xlen_t offset, R_xlen_t n)
{
	switch (Rtype) {
	    case INTSXP: {
		_set_elts_to_val(Rtype, x, offset, n, &intNA);
		return;
	    }
	    case LGLSXP: {
		const int val = NA_LOGICAL;
		_set_elts_to_val(Rtype, x, offset, n, &val);
		return;
	    }
	    case REALSXP: {
		_set_elts_to_val(Rtype, x, offset, n, &doubleNA);
		return;
	    }
	    case CPLXSXP: {
		_set_elts_to_val(Rtype, x, offset, n, &RcomplexNA);
		return;
	    }
	}
	error("SparseArray internal error in _set_elts_to_NA():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return;
}


static void set_character_elts_to_val(SEXP Rvector,
		R_xlen_t subvec_offset, int subvec_len, const SEXP val)
{
	R_xlen_t i = subvec_offset + subvec_len - 1;
	for ( ; i >= subvec_offset; i--)
		SET_STRING_ELT(Rvector, i, val);
	return;
}

static void set_list_elts_to_val(SEXP Rvector,
		R_xlen_t subvec_offset, int subvec_len, const SEXP val)
{
	R_xlen_t i = subvec_offset + subvec_len - 1;
	for ( ; i >= subvec_offset; i--)
		SET_VECTOR_ELT(Rvector, i, val);
	return;
}

/* When 'Rtype' is STRSXP or VECSXP, 'val' must be an SEXP. Otherwise, it must
   be a pointer to an int, double, Rcomplex, or Rbyte. */
void _set_Rsubvec_elts_to_val(SEXP Rvector,
		R_xlen_t subvec_offset, R_xlen_t subvec_len, const void *val)
{
	SEXPTYPE Rtype = TYPEOF(Rvector);
	if (Rtype == STRSXP) {
		set_character_elts_to_val(Rvector, subvec_offset, subvec_len,
					  (const SEXP) val);
		return;
	}
	if (Rtype == VECSXP) {
		set_list_elts_to_val(Rvector, subvec_offset, subvec_len,
				     (const SEXP) val);
		return;
	}
	_set_elts_to_val(TYPEOF(Rvector), DATAPTR(Rvector),
			 subvec_offset, subvec_len, val);
}

void _set_Rsubvec_elts_to_zero(SEXP Rvector,
		R_xlen_t subvec_offset, R_xlen_t subvec_len)
{
	SEXPTYPE Rtype = TYPEOF(Rvector);
	if (Rtype == STRSXP) {
		set_character_elts_to_val(Rvector, subvec_offset, subvec_len,
					  R_BlankString);
		return;
	}
	if (Rtype == VECSXP) {
		set_list_elts_to_val(Rvector, subvec_offset, subvec_len,
				     R_NilValue);
		return;
	}
	_set_elts_to_zero(Rtype, DATAPTR(Rvector), subvec_offset, subvec_len);
	return;
}

/* Restricted to types "integer", "logical", "double", "complex", and "raw". */
void _set_Rsubvec_elts_to_one(SEXP Rvector,
		R_xlen_t subvec_offset, R_xlen_t subvec_len)
{
	_set_elts_to_one(TYPEOF(Rvector), DATAPTR(Rvector),
			 subvec_offset, subvec_len);
	return;
}

/* Restricted to types "integer", "double", and "complex". */
void _set_Rsubvec_elts_to_minus_one(SEXP Rvector,
		R_xlen_t subvec_offset, R_xlen_t subvec_len)
{
	_set_elts_to_minus_one(TYPEOF(Rvector), DATAPTR(Rvector),
			 subvec_offset, subvec_len);
	return;
}

/* Restricted to types "integer", "logical", "double", "complex",
   and "character". */
void _set_Rsubvec_elts_to_NA(SEXP Rvector,
		R_xlen_t subvec_offset, R_xlen_t subvec_len)
{
	SEXPTYPE Rtype = TYPEOF(Rvector);
	if (Rtype == STRSXP) {
		set_character_elts_to_val(Rvector, subvec_offset, subvec_len,
					  NA_STRING);
		return;
	}
	_set_elts_to_NA(Rtype, DATAPTR(Rvector), subvec_offset, subvec_len);
	return;
}

/* When 'Rtype' is STRSXP or VECSXP, 'val' must be an SEXP. Otherwise, it must
   be a pointer to an int, double, Rcomplex, or Rbyte. */
void _set_Rvector_elts_to_val(SEXP Rvector, const void *val)
{
	_set_Rsubvec_elts_to_val(Rvector, 0, XLENGTH(Rvector), val);
	return;
}

void _set_Rvector_elts_to_zero(SEXP Rvector)
{
	_set_Rsubvec_elts_to_zero(Rvector, 0, XLENGTH(Rvector));
	return;
}

/* Restricted to types "integer", "logical", "double", "complex", and "raw". */
void _set_Rvector_elts_to_one(SEXP Rvector)
{
	_set_Rsubvec_elts_to_one(Rvector, 0, XLENGTH(Rvector));
	return;
}

/* Restricted to types "integer", "double", and "complex". */
void _set_Rvector_elts_to_minus_one(SEXP Rvector)
{
	_set_Rsubvec_elts_to_minus_one(Rvector, 0, XLENGTH(Rvector));
	return;
}

/* Restricted to types "integer", "logical", "double", "complex",
   and "character". */
void _set_Rvector_elts_to_NA(SEXP Rvector)
{
	_set_Rsubvec_elts_to_NA(Rvector, 0, XLENGTH(Rvector));
	return;
}


/****************************************************************************
 * _set_selected_elts_to_zero()
 * _set_selected_elts_to_one()
 * _set_selected_Rvector_elts_to_zero()
 * _set_selected_Rvector_elts_to_one()
 */

#define	DEFINE_set_selected_TYPE_elts_to_val_FUN(type)			\
static void set_selected_ ## type ## _elts_to_val			\
	(type *x, const int *selection, int n, type val)		\
{									\
	for (int i = 0; i < n; i++, selection++)			\
		x[*selection] = val;					\
	return;								\
}

DEFINE_set_selected_TYPE_elts_to_val_FUN(int)
DEFINE_set_selected_TYPE_elts_to_val_FUN(double)
DEFINE_set_selected_TYPE_elts_to_val_FUN(Rcomplex)
DEFINE_set_selected_TYPE_elts_to_val_FUN(Rbyte)

#define	set_selected_TYPE_elts_to_val(type, x, offset, selection, n, val) \
	set_selected_ ## type ## _elts_to_val( \
				SHIFT_DATAPTR(type, (x), (offset)), \
				(selection), (n), (val))

/* Restricted to types "integer", "logical", "double", "complex", and "raw". */
void _set_selected_elts_to_zero(SEXPTYPE Rtype, void *x, R_xlen_t offset,
		const int *selection, int n)
{
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		set_selected_TYPE_elts_to_val(int, x, offset,
					      selection, n, int0);
		return;
	    case REALSXP:
		set_selected_TYPE_elts_to_val(double, x, offset,
					      selection, n, double0);
		return;
	    case CPLXSXP:
		set_selected_TYPE_elts_to_val(Rcomplex, x, offset,
					      selection, n, Rcomplex0);
		return;
	    case RAWSXP:
		set_selected_TYPE_elts_to_val(Rbyte, x, offset,
					      selection, n, Rbyte0);
		return;
	}
	error("SparseArray internal error in _set_selected_elts_to_zero():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return;
}

/* Restricted to types "integer", "logical", "double", "complex", and "raw". */
void _set_selected_elts_to_one(SEXPTYPE Rtype, void *x, R_xlen_t offset,
		const int *selection, int n)
{
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		set_selected_TYPE_elts_to_val(int, x, offset,
					      selection, n, int1);
		return;
	    case REALSXP:
		set_selected_TYPE_elts_to_val(double, x, offset,
					      selection, n, double1);
		return;
	    case CPLXSXP:
		set_selected_TYPE_elts_to_val(Rcomplex, x, offset,
					      selection, n, Rcomplex1);
		return;
	    case RAWSXP:
		set_selected_TYPE_elts_to_val(Rbyte, x, offset,
					      selection, n, Rbyte1);
		return;
	}
	error("SparseArray internal error in _set_selected_elts_to_one():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return;
}

static void set_selected_character_elts(SEXP Rvector, R_xlen_t subvec_offset,
		const int *selection, int n, SEXP val)
{
	for (int i = 0; i < n; i++, selection++)
		SET_STRING_ELT(Rvector, subvec_offset + *selection, val);
	return;
}

static void set_selected_list_elts(SEXP Rvector, R_xlen_t subvec_offset,
		const int *selection, int n, SEXP val)
{
	for (int i = 0; i < n; i++, selection++)
		SET_VECTOR_ELT(Rvector, subvec_offset + *selection, val);
	return;
}

void _set_selected_Rsubvec_elts_to_zero(SEXP Rvector, R_xlen_t subvec_offset,
		const int *selection, int n)
{
	SEXPTYPE Rtype = TYPEOF(Rvector);
	if (Rtype == STRSXP) {
		set_selected_character_elts(Rvector, subvec_offset,
				selection, n, R_BlankString);
		return;
	}
	if (Rtype == VECSXP) {
		set_selected_list_elts(Rvector, subvec_offset,
				selection, n, R_NilValue);
		return;
	}
	_set_selected_elts_to_zero(Rtype, DATAPTR(Rvector), subvec_offset,
				selection, n);
	return;
}

/* Restricted to types "integer", "logical", "double", "complex", and "raw". */
void _set_selected_Rsubvec_elts_to_one(SEXP Rvector, R_xlen_t subvec_offset,
		const int *selection, int n)
{
	_set_selected_elts_to_one(TYPEOF(Rvector), DATAPTR(Rvector),
				subvec_offset,
				selection, n);
	return;
}


/****************************************************************************
 * _new_Rvector0()
 * _new_Rmatrix0()
 * _new_Rarray0()
 * _new_Rvector1()
 * _new_RvectorNA()
 * _new_RarrayNA()
 */

/* Like allocVector() but with initialization of the vector elements. */
SEXP _new_Rvector0(SEXPTYPE Rtype, R_xlen_t len)
{
	SEXP ans = PROTECT(allocVector(Rtype, len));
	/* allocVector() does NOT initialize the vector elements, except for
	   a character vector or a list. */
	if (Rtype != STRSXP && Rtype != VECSXP)
		_set_Rvector_elts_to_zero(ans);
	UNPROTECT(1);
	return ans;
}

/* Like allocMatrix() but with initialization of the matrix elements and
   addition of the dimnames. */
SEXP _new_Rmatrix0(SEXPTYPE Rtype, int nrow, int ncol, SEXP dimnames)
{
	SEXP ans = PROTECT(allocMatrix(Rtype, nrow, ncol));
	/* allocMatrix() is just a thin wrapper around allocVector() and
	   the latter does NOT initialize the vector elements, except for
	   a character vector or a list. */
	if (Rtype != STRSXP && Rtype != VECSXP)
		_set_Rvector_elts_to_zero(ans);
	SET_DIMNAMES(ans, dimnames);
	UNPROTECT(1);
	return ans;
}

/* Like allocArray() but with initialization of the array elements and
   addition of the dimnames. */
SEXP _new_Rarray0(SEXPTYPE Rtype, SEXP dim, SEXP dimnames)
{
	SEXP ans = PROTECT(allocArray(Rtype, dim));
	/* allocArray() is just a thin wrapper around allocVector() and
	   the latter does NOT initialize the vector elements, except for
	   a character vector or a list. */
	if (Rtype != STRSXP && Rtype != VECSXP)
		_set_Rvector_elts_to_zero(ans);
	SET_DIMNAMES(ans, dimnames);
	UNPROTECT(1);
	return ans;
}

/* Like _new_Rvector0() but:
   - initializes the vector elements to 1;
   - restricted to types "integer", "logical", "double", "complex", and "raw";
   - 'len' must be int, not R_xlen_t (i.e. long vectors not supported). */
SEXP _new_Rvector1(SEXPTYPE Rtype, int len)
{
	SEXP ans = PROTECT(allocVector(Rtype, (R_xlen_t) len));
	_set_Rvector_elts_to_one(ans);
	UNPROTECT(1);
	return ans;
}

/* Like _new_Rvector0() but:
   - initializes the array elements to NA;
   - restricted to types "integer", "logical", "double", "complex",
     and "character". */
SEXP _new_RvectorNA(SEXPTYPE Rtype, R_xlen_t len)
{
	SEXP ans = PROTECT(allocVector(Rtype, len));
	_set_Rvector_elts_to_NA(ans);
	UNPROTECT(1);
	return ans;
}

/* Like _new_Rarray0() but:
   - initializes the array elements to NA;
   - restricted to types "integer", "logical", "double", "complex",
     and "character". */
SEXP _new_RarrayNA(SEXPTYPE Rtype, SEXP dim, SEXP dimnames)
{
	SEXP ans = PROTECT(allocArray(Rtype, dim));
	_set_Rvector_elts_to_NA(ans);
	SET_DIMNAMES(ans, dimnames);
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * _collect_offsets_of_nonzero_Rsubvec_elts()
 * _collect_offsets_of_nonNA_Rsubvec_elts()
 */

static int collect_offsets_of_nonzero_int_elts(
		const int *x, int n, int *out)
{
	const int *out0 = out;
	for (int i = 0; i < n; i++)
		if (x[i] != int0)
			*(out++) = i;
	return (int) (out - out0);
}

static int collect_offsets_of_nonNA_int_elts(
		const int *x, int n, int *out)
{
	const int *out0 = out;
	for (int i = 0; i < n; i++)
		if (x[i] != NA_INTEGER)
			*(out++) = i;
	return (int) (out - out0);
}

static int collect_offsets_of_nonzero_double_elts(
		const double *x, int n, int *out)
{
	const int *out0 = out;
	for (int i = 0; i < n; i++)
		if (x[i] != double0)
			*(out++) = i;
	return (int) (out - out0);
}

static int collect_offsets_of_nonNA_double_elts(
		const double *x, int n, int *out)
{
	const int *out0 = out;
	for (int i = 0; i < n; i++)
		if (!R_IsNA(x[i]))  // do NOT use ISNAN()!
			*(out++) = i;
	return (int) (out - out0);
}

static int collect_offsets_of_nonzero_Rcomplex_elts(
		const Rcomplex *x, int n, int *out)
{
	const int *out0 = out;
	for (int i = 0; i < n; i++, x++) {
		if (x->r != Rcomplex0.r || x->i != Rcomplex0.i)
			*(out++) = i;
	}
	return (int) (out - out0);
}

static int collect_offsets_of_nonNA_Rcomplex_elts(
		const Rcomplex *x, int n, int *out)
{
	const int *out0 = out;
	for (int i = 0; i < n; i++, x++) {
		if (!(R_IsNA(x->r) || R_IsNA(x->i)))
			*(out++) = i;
	}
	return (int) (out - out0);
}

static int collect_offsets_of_nonzero_Rbyte_elts(
		const Rbyte *x, int n, int *out)
{
	const int *out0 = out;
	for (int i = 0; i < n; i++)
		if (x[i] != Rbyte0)
			*(out++) = i;
	return (int) (out - out0);
}

static int collect_offsets_of_nonempty_character_elts(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *out)
{
	const int *out0 = out;
	for (int i = 0; i < subvec_len; i++, subvec_offset++) {
		SEXP vec_elt = STRING_ELT(Rvector, subvec_offset);
		if (!IS_EMPTY_CHARSXP(vec_elt))
			*(out++) = i;
	}
	return (int) (out - out0);
}

static int collect_offsets_of_nonNA_character_elts(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *out)
{
	const int *out0 = out;
	for (int i = 0; i < subvec_len; i++, subvec_offset++) {
		SEXP vec_elt = STRING_ELT(Rvector, subvec_offset);
		if (vec_elt != NA_STRING)
			*(out++) = i;
	}
	return (int) (out - out0);
}

static int collect_offsets_of_nonnull_list_elts(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *out)
{
	const int *out0 = out;
	for (int i = 0; i < subvec_len; i++, subvec_offset++) {
		SEXP vec_elt = VECTOR_ELT(Rvector, subvec_offset);
		if (vec_elt != R_NilValue)
			*(out++) = i;
	}
	return (int) (out - out0);
}

/* Only looks at the subvector of 'Rvector' made of the range of elements
   defined by 'subvec_offset' and 'subvec_len'.
   Offsets of nonzero elements are collected with respect to this subvector.
   Caller must make sure that the 'out' array is big enough to store all
   the collected offsets. Safe choice is to allocate an array of 'subvec_len'
   integers.
   Note that even though 'Rvector' can be a long vector, the subvector
   defined by 'subvec_offset/subvec_len' cannot i.e. 'subvec_len' must be
   supplied as an int.
   Returns the number of collected offsets. */
int _collect_offsets_of_nonzero_Rsubvec_elts(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *out)
{
	SEXPTYPE Rtype = TYPEOF(Rvector);
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		return collect_offsets_of_nonzero_int_elts(
				INTEGER(Rvector) + subvec_offset,
				subvec_len, out);
	    case REALSXP:
		return collect_offsets_of_nonzero_double_elts(
				REAL(Rvector) + subvec_offset,
				subvec_len, out);
	    case CPLXSXP:
		return collect_offsets_of_nonzero_Rcomplex_elts(
				COMPLEX(Rvector) + subvec_offset,
				subvec_len, out);
	    case RAWSXP:
		return collect_offsets_of_nonzero_Rbyte_elts(
				RAW(Rvector) + subvec_offset,
				subvec_len, out);
	    case STRSXP:
		return collect_offsets_of_nonempty_character_elts(
				Rvector, subvec_offset,
				subvec_len, out);
	    case VECSXP:
		return collect_offsets_of_nonnull_list_elts(
				Rvector, subvec_offset,
				subvec_len, out);
	}
	error("SparseArray internal error in "
	      "_collect_offsets_of_nonzero_Rsubvec_elts():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return -1;  /* will never reach this */
}

/* Similar to _collect_offsets_of_nonzero_Rsubvec_elts() above but collects
   non-NAs instead of nonzeros. Restricted to types "integer", "logical",
   "double", "complex", and "character". */
int _collect_offsets_of_nonNA_Rsubvec_elts(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *out)
{
	SEXPTYPE Rtype = TYPEOF(Rvector);
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		return collect_offsets_of_nonNA_int_elts(
				INTEGER(Rvector) + subvec_offset,
				subvec_len, out);
	    case REALSXP:
		return collect_offsets_of_nonNA_double_elts(
				REAL(Rvector) + subvec_offset,
				subvec_len, out);
	    case CPLXSXP:
		return collect_offsets_of_nonNA_Rcomplex_elts(
				COMPLEX(Rvector) + subvec_offset,
				subvec_len, out);
	    case STRSXP:
		return collect_offsets_of_nonNA_character_elts(
				Rvector, subvec_offset,
				subvec_len, out);
	}
	error("SparseArray internal error in "
	      "_collect_offsets_of_nonNA_Rsubvec_elts():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return -1;  /* will never reach this */
}


/****************************************************************************
 * _all_elts_equal_one()
 * _all_Rsubvec_elts_equal_one()
 * _all_Rvector_elts_equal_one()
 * _all_selected_Rsubvec_elts_equal_one()
 */

static int all_int_elts_equal_one(const int *x, int n)
{
	for (int i = 0; i < n; i++)
		if (x[i] != int1)
			return 0;
	return 1;
}

static int all_double_elts_equal_one(const double *x, int n)
{
	for (int i = 0; i < n; i++)
		if (x[i] != double1)
			return 0;
	return 1;
}

static int all_Rcomplex_elts_equal_one(const Rcomplex *x, int n)
{
	for (int i = 0; i < n; i++, x++)
		if (x->r != Rcomplex1.r || x->i != Rcomplex1.i)
			return 0;
	return 1;
}

static int all_Rbyte_elts_equal_one(const Rbyte *x, int n)
{
	for (int i = 0; i < n; i++)
		if (x[i] != Rbyte1)
			return 0;
	return 1;
}

/* Always returns 0 on a character vector or list at the moment. */
int _all_elts_equal_one(SEXPTYPE Rtype, const void *x, int n)
{
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		return all_int_elts_equal_one((const int *) x, n);
	    case REALSXP:
		return all_double_elts_equal_one((const double *) x, n);
	    case CPLXSXP:
		return all_Rcomplex_elts_equal_one((const Rcomplex *) x, n);
	    case RAWSXP:
		return all_Rbyte_elts_equal_one((const Rbyte *) x, n);
	    case STRSXP: case VECSXP:
		return 0;
	}
	error("SparseArray internal error in "
	      "_all_elts_equal_one():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
}

/* Always returns 0 on a character vector or list at the moment. */
int _all_Rsubvec_elts_equal_one(SEXP Rvector,
		R_xlen_t subvec_offset, int subvec_len)
{
	SEXPTYPE Rtype = TYPEOF(Rvector);
	switch (Rtype) {
	    case INTSXP: case LGLSXP: {
		const int      *x = INTEGER(Rvector) + subvec_offset;
		return all_int_elts_equal_one(x, subvec_len);
	    }
	    case REALSXP: {
		const double   *x = REAL(Rvector) + subvec_offset;
		return all_double_elts_equal_one(x, subvec_len);
	    }
	    case CPLXSXP: {
		const Rcomplex *x = COMPLEX(Rvector) + subvec_offset;
		return all_Rcomplex_elts_equal_one(x, subvec_len);
	    }
	    case RAWSXP: {
		const Rbyte    *x = RAW(Rvector) + subvec_offset;
		return all_Rbyte_elts_equal_one(x, subvec_len);
	    }
	    case STRSXP: case VECSXP:
		return 0;
	}
	error("SparseArray internal error in "
	      "_all_Rsubvec_elts_equal_one():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
}

/* Always returns 0 on a character vector or list at the moment. */
int _all_Rvector_elts_equal_one(SEXP Rvector)
{
	return _all_Rsubvec_elts_equal_one(Rvector, 0, XLENGTH(Rvector));
}

static int all_selected_int_elts_equal_one(const int *x,
		const int *selection, int n)
{
	for (int k = 0; k < n; k++, selection++)
		if (x[*selection] != int1)
			return 0;
	return 1;
}

static int all_selected_double_elts_equal_one(const double *x,
		const int *selection, int n)
{
	for (int k = 0; k < n; k++, selection++)
		if (x[*selection] != double1)
			return 0;
	return 1;
}

static int all_selected_Rcomplex_elts_equal_one(const Rcomplex *x,
		const int *selection, int n)
{
	for (int k = 0; k < n; k++, selection++) {
		const Rcomplex *z = x + *selection;
		if (z->r != Rcomplex1.r || z->i != Rcomplex1.i)
			return 0;
	}
	return 1;
}

static int all_selected_Rbyte_elts_equal_one(const Rbyte *x,
		const int *selection, int n)
{
	for (int k = 0; k < n; k++, selection++)
		if (x[*selection] != Rbyte1)
			return 0;
	return 1;
}

/* Always returns 0 on a character vector or list at the moment. */
int _all_selected_Rsubvec_elts_equal_one(SEXP Rvector, R_xlen_t subvec_offset,
		const int *selection, int n)
{
	SEXPTYPE Rtype = TYPEOF(Rvector);
	switch (Rtype) {
	    case INTSXP: case LGLSXP: {
		const int      *x = INTEGER(Rvector) + subvec_offset;
		return all_selected_int_elts_equal_one(x, selection, n);
	    }
	    case REALSXP: {
		const double   *x = REAL(Rvector) + subvec_offset;
		return all_selected_double_elts_equal_one(x, selection, n);
	    }
	    case CPLXSXP: {
		const Rcomplex *x = COMPLEX(Rvector) + subvec_offset;
		return all_selected_Rcomplex_elts_equal_one(x, selection, n);
	    }
	    case RAWSXP: {
		const Rbyte    *x = RAW(Rvector) + subvec_offset;
		return all_selected_Rbyte_elts_equal_one(x, selection, n);
	    }
	    case STRSXP: case VECSXP:
		return 0;
	}
	error("SparseArray internal error in "
	      "_all_selected_Rsubvec_elts_equal_one():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return 0;  /* will never reach this */
}


/****************************************************************************
 * _select_copy_Rvector_elt_FUN()
 */

CopyRVectorEltFUN _select_copy_Rvector_elt_FUN(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case INTSXP: case LGLSXP: return copy_INTEGER_elt;
	    case REALSXP:             return copy_NUMERIC_elt;
	    case CPLXSXP:             return copy_COMPLEX_elt;
	    case RAWSXP:              return copy_RAW_elt;
	    case STRSXP:              return copy_CHARACTER_elt;
	    case VECSXP:              return copy_LIST_elt;
	}
	return NULL;
}


/****************************************************************************
 * _copy_Rvector_elts()
 */

#define	DEFINE_copy_TYPE_elts_FUN(type)					\
static void copy_ ## type ## _elts					\
	(const type *src, type *dest, R_xlen_t nelt)			\
{									\
	memcpy(dest, src, sizeof(type) * nelt);				\
	return;								\
}

DEFINE_copy_TYPE_elts_FUN(int)
DEFINE_copy_TYPE_elts_FUN(double)
DEFINE_copy_TYPE_elts_FUN(Rcomplex)
DEFINE_copy_TYPE_elts_FUN(Rbyte)

static void copy_character_elts(SEXP in,  R_xlen_t in_offset,
				SEXP out, R_xlen_t out_offset,
				R_xlen_t nelt)
{
	for (R_xlen_t k = 0; k < nelt; k++, in_offset++, out_offset++)
		SET_STRING_ELT(out, out_offset, STRING_ELT(in, in_offset));
        return;
}

static void copy_list_elts(SEXP in,  R_xlen_t in_offset,
			   SEXP out, R_xlen_t out_offset,
			   R_xlen_t nelt)
{
	for (R_xlen_t k = 0; k < nelt; k++, in_offset++, out_offset++)
		SET_VECTOR_ELT(out, out_offset, VECTOR_ELT(in, in_offset));
        return;
}

void _copy_Rvector_elts(SEXP in,  R_xlen_t in_offset,
			SEXP out, R_xlen_t out_offset,
			R_xlen_t nelt)
{
	SEXPTYPE Rtype = TYPEOF(in);
	if (TYPEOF(out) != Rtype)
		error("SparseArray internal error in _copy_Rvector_elts():\n"
		      "    TYPEOF(in) != TYPEOF(out)");
	if (in_offset + nelt > XLENGTH(in))
		error("SparseArray internal error in _copy_Rvector_elts():\n"
		      "    in_offset + nelt > XLENGTH(in)");
	if (out_offset + nelt > XLENGTH(out))
		error("SparseArray internal error in _copy_Rvector_elts():\n"
		      "    out_offset + nelt > XLENGTH(out)");
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		copy_int_elts(INTEGER(in)  + in_offset,
			INTEGER(out) + out_offset, nelt);
		return;
	    case REALSXP:
		copy_double_elts(REAL(in)  + in_offset,
			REAL(out) + out_offset, nelt);
		return;
	    case CPLXSXP:
		copy_Rcomplex_elts(COMPLEX(in)  + in_offset,
			COMPLEX(out) + out_offset, nelt);
		return;
	    case RAWSXP:
		copy_Rbyte_elts(RAW(in)  + in_offset,
			RAW(out) + out_offset, nelt);
		return;
	    case STRSXP:
		copy_character_elts(in, in_offset, out, out_offset, nelt);
		return;
	    case VECSXP:
		copy_list_elts(in, in_offset, out, out_offset, nelt);
		return;
	}
	error("SparseArray internal error in _copy_Rvector_elts():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return;
}


/****************************************************************************
 * _copy_selected_Rsubvec_elts()
 */

void _copy_selected_int_elts(const int *in,
		const int *selection, int n, int *out)
{
	for (int k = 0; k < n; k++, selection++, out++)
		*out = in[*selection];
	return;
}

void _copy_selected_double_elts(const double *in,
		const int *selection, int n, double *out)
{
	for (int k = 0; k < n; k++, selection++, out++)
		*out = in[*selection];
	return;
}

void _copy_selected_Rcomplex_elts(const Rcomplex *in,
		const int *selection, int n, Rcomplex *out)
{
	for (int k = 0; k < n; k++, selection++, out++)
		*out = in[*selection];
	return;
}

void _copy_selected_Rbyte_elts(const Rbyte *in,
		const int *selection, int n, Rbyte *out)
{
	for (int k = 0; k < n; k++, selection++, out++)
		*out = in[*selection];
	return;
}

void _copy_selected_character_elts(SEXP in, R_xlen_t in_offset,
		const int *selection, int n, SEXP out)
{
	for (int k = 0; k < n; k++, selection++) {
		R_xlen_t i = in_offset + *selection;
		SET_STRING_ELT(out, k, STRING_ELT(in, i));
	}
	return;
}

void _copy_selected_list_elts(SEXP in, R_xlen_t in_offset,
		const int *selection, int n, SEXP out)
{
	for (int k = 0; k < n; k++, selection++) {
		R_xlen_t i = in_offset + *selection;
		SET_VECTOR_ELT(out, k, VECTOR_ELT(in, i));
	}
	return;
}

/* The selection is assumed to have the same length as 'out_Rvector'.
   Only for a 'selection' and 'out_Rvector' of length <= INT_MAX.
   Do NOT use on a 'selection' or 'out_Rvector' of length > INT_MAX! */
void _copy_selected_Rsubvec_elts(
		SEXP Rvector, R_xlen_t subvec_offset,
		const int *selection, SEXP out_Rvector)
{
	SEXPTYPE Rtype = TYPEOF(Rvector);
	int out_len = LENGTH(out_Rvector);  /* assumed to be also the length
					       of the selection */
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		_copy_selected_int_elts(INTEGER(Rvector) + subvec_offset,
				selection, out_len, INTEGER(out_Rvector));
		return;
	    case REALSXP:
		_copy_selected_double_elts(REAL(Rvector) + subvec_offset,
				selection, out_len, REAL(out_Rvector));
		return;
	    case CPLXSXP:
		_copy_selected_Rcomplex_elts(COMPLEX(Rvector) + subvec_offset,
				selection, out_len, COMPLEX(out_Rvector));
		return;
	    case RAWSXP:
		_copy_selected_Rbyte_elts(RAW(Rvector) + subvec_offset,
				selection, out_len, RAW(out_Rvector));
		return;
	    case STRSXP:
		_copy_selected_character_elts(Rvector, subvec_offset,
				selection, out_len, out_Rvector);
		return;
	    case VECSXP:
		_copy_selected_list_elts(Rvector, subvec_offset,
				selection, out_len, out_Rvector);
		return;
	}
	error("SparseArray internal error in _copy_selected_Rsubvec_elts():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
}


/****************************************************************************
 * _subset_Rsubvec()
 */

/* WARNING: 'selection' cannot contain NAs or bad things will happen. This
   is NOT checked! */
SEXP _subset_Rsubvec(SEXP Rvector, R_xlen_t subvec_offset,
		const int *selection, int n)
{
	SEXP ans = PROTECT(allocVector(TYPEOF(Rvector), n));
	_copy_selected_Rsubvec_elts(Rvector, subvec_offset,
				    selection, ans);
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * _copy_Rvector_elts_to_offsets()
 */

void _copy_int_elts_to_offsets(const int *in,
		const int *selection, int n, int *out)
{
	for (int k = 0; k < n; k++, selection++, in++)
		out[*selection] = *in;
	return;
}

void _copy_double_elts_to_offsets(const double *in,
		const int *selection, int n, double *out)
{
	for (int k = 0; k < n; k++, selection++, in++)
		out[*selection] = *in;
	return;
}

void _copy_Rcomplex_elts_to_offsets(const Rcomplex *in,
		const int *selection, int n, Rcomplex *out)
{
	for (int k = 0; k < n; k++, selection++, in++)
		out[*selection] = *in;
	return;
}

void _copy_Rbyte_elts_to_offsets(const Rbyte *in,
		const int *selection, int n, Rbyte *out)
{
	for (int k = 0; k < n; k++, selection++, in++)
		out[*selection] = *in;
	return;
}

void _copy_character_elts_to_offsets(SEXP in,
		const int *selection, int n,
		SEXP out, R_xlen_t out_offset)
{
	for (int k = 0; k < n; k++, selection++) {
		R_xlen_t j = out_offset + *selection;
		SET_STRING_ELT(out, j, STRING_ELT(in, k));
	}
	return;
}

void _copy_list_elts_to_offsets(SEXP in,
		const int *selection, int n,
		SEXP out, R_xlen_t out_offset)
{
	for (int k = 0; k < n; k++, selection++) {
		R_xlen_t j = out_offset + *selection;
		SET_VECTOR_ELT(out, j, VECTOR_ELT(in, k));
	}
	return;
}

/* The selection is assumed to have the same length as 'Rvector'.
   Only for a 'selection' and 'Rvector' of length <= INT_MAX.
   Do NOT use on a 'selection' or 'Rvector' of length > INT_MAX! */
void _copy_Rvector_elts_to_offsets(
		SEXP Rvector, const int *selection,
		SEXP out_Rvector, R_xlen_t out_offset)
{
	SEXPTYPE Rtype = TYPEOF(Rvector);
	int in_len = LENGTH(Rvector);  /* assumed to be also the length
					  of the selection */
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		_copy_int_elts_to_offsets(INTEGER(Rvector),
				selection, in_len,
				INTEGER(out_Rvector) + out_offset);
		return;
	    case REALSXP:
		_copy_double_elts_to_offsets(REAL(Rvector),
				selection, in_len,
				REAL(out_Rvector) + out_offset);
		return;
	    case CPLXSXP:
		_copy_Rcomplex_elts_to_offsets(COMPLEX(Rvector),
				selection, in_len,
				COMPLEX(out_Rvector) + out_offset);
		return;
	    case RAWSXP:
		_copy_Rbyte_elts_to_offsets(RAW(Rvector),
				selection, in_len,
				RAW(out_Rvector) + out_offset);
		return;
	    case STRSXP:
		_copy_character_elts_to_offsets(Rvector,
				selection, in_len,
				out_Rvector, out_offset);
		return;
	    case VECSXP:
		_copy_list_elts_to_offsets(Rvector,
				selection, in_len,
				out_Rvector, out_offset);
		return;
	}
	error("SparseArray internal error in "
	      "_copy_Rvector_elts_to_offsets():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return;
}


/****************************************************************************
 * _copy_Rvector_elts_from_selected_offsets()
 * _copy_Rvector_elts_from_selected_lloffsets()
 */

static void copy_int_elts_from_selected_offsets(const int *in,
		const int *offsets, const int *offset_selection,
		int n, int *out)
{
	for (int k = 0; k < n; k++, offset_selection++, out++)
		*out = in[offsets[*offset_selection]];
	return;
}
static void copy_int_elts_from_selected_lloffsets(const int *in,
		const long long *lloffsets, const int *lloffset_selection,
		int n, int *out)
{
	for (int k = 0; k < n; k++, lloffset_selection++, out++)
		*out = in[lloffsets[*lloffset_selection]];
	return;
}

static void copy_double_elts_from_selected_offsets(const double *in,
		const int *offsets, const int *offset_selection,
		int n, double *out)
{
	for (int k = 0; k < n; k++, offset_selection++, out++)
		*out = in[offsets[*offset_selection]];
	return;
}
static void copy_double_elts_from_selected_lloffsets(const double *in,
		const long long *lloffsets, const int *lloffset_selection,
		int n, double *out)
{
	for (int k = 0; k < n; k++, lloffset_selection++, out++)
		*out = in[lloffsets[*lloffset_selection]];
	return;
}

static void copy_Rcomplex_elts_from_selected_offsets(const Rcomplex *in,
		const int *offsets, const int *offset_selection,
		int n, Rcomplex *out)
{
	for (int k = 0; k < n; k++, offset_selection++, out++)
		*out = in[offsets[*offset_selection]];
	return;
}
static void copy_Rcomplex_elts_from_selected_lloffsets(const Rcomplex *in,
		const long long *lloffsets, const int *lloffset_selection,
		int n, Rcomplex *out)
{
	for (int k = 0; k < n; k++, lloffset_selection++, out++)
		*out = in[lloffsets[*lloffset_selection]];
	return;
}

static void copy_Rbyte_elts_from_selected_offsets(const Rbyte *in,
		const int *offsets, const int *offset_selection,
		int n, Rbyte *out)
{
	for (int k = 0; k < n; k++, offset_selection++, out++)
		*out = in[offsets[*offset_selection]];
	return;
}
static void copy_Rbyte_elts_from_selected_lloffsets(const Rbyte *in,
		const long long *lloffsets, const int *lloffset_selection,
		int n, Rbyte *out)
{
	for (int k = 0; k < n; k++, lloffset_selection++, out++)
		*out = in[lloffsets[*lloffset_selection]];
	return;
}

static void copy_character_elts_from_selected_offsets(SEXP in,
		const int *offsets, const int *offset_selection,
		int n, SEXP out)
{
	for (int k = 0; k < n; k++, offset_selection++)
		SET_STRING_ELT(out, k,
			STRING_ELT(in, offsets[*offset_selection]));
	return;
}
static void copy_character_elts_from_selected_lloffsets(SEXP in,
		const long long *lloffsets, const int *lloffset_selection,
		int n, SEXP out)
{
	for (int k = 0; k < n; k++, lloffset_selection++)
		SET_STRING_ELT(out, k,
			STRING_ELT(in, lloffsets[*lloffset_selection]));
	return;
}

static void copy_list_elts_from_selected_offsets(SEXP in,
		const int *offsets, const int *offset_selection,
		int n, SEXP out)
{
	for (int k = 0; k < n; k++, offset_selection++)
		SET_VECTOR_ELT(out, k,
			VECTOR_ELT(in, offsets[*offset_selection]));
	return;
}
static void copy_list_elts_from_selected_lloffsets(SEXP in,
		const long long *lloffsets, const int *lloffset_selection,
		int n, SEXP out)
{
	for (int k = 0; k < n; k++, lloffset_selection++)
		SET_VECTOR_ELT(out, k,
			VECTOR_ELT(in, lloffsets[*lloffset_selection]));
	return;
}

/* The selection is assumed to have the same length as 'out_Rvector'.
   Only for an 'offset_selection' and 'out_Rvector' of length <= INT_MAX.
   Don't use on an 'offset_selection' or 'out_Rvector' of length > INT_MAX! */
void _copy_Rvector_elts_from_selected_offsets(SEXP Rvector,
		const int *offsets, const int *offset_selection,
		SEXP out_Rvector)
{
	SEXPTYPE Rtype = TYPEOF(Rvector);
	int out_len = LENGTH(out_Rvector);  /* assumed to be also the length
					       of the selection */
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		copy_int_elts_from_selected_offsets(INTEGER(Rvector),
				offsets, offset_selection, out_len,
				INTEGER(out_Rvector));
		return;
	    case REALSXP:
		copy_double_elts_from_selected_offsets(REAL(Rvector),
				offsets, offset_selection, out_len,
				REAL(out_Rvector));
		return;
	    case CPLXSXP:
		copy_Rcomplex_elts_from_selected_offsets(COMPLEX(Rvector),
				offsets, offset_selection, out_len,
				COMPLEX(out_Rvector));
		return;
	    case RAWSXP:
		copy_Rbyte_elts_from_selected_offsets(RAW(Rvector),
				offsets, offset_selection, out_len,
				RAW(out_Rvector));
		return;
	    case STRSXP:
		copy_character_elts_from_selected_offsets(Rvector,
				offsets, offset_selection, out_len,
				out_Rvector);
		return;
	    case VECSXP:
		copy_list_elts_from_selected_offsets(Rvector,
				offsets, offset_selection, out_len,
				out_Rvector);
		return;
	}
	error("SparseArray internal error in "
	      "_copy_Rvector_elts_from_selected_offsets():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
}

/* The selection is assumed to have the same length as 'out_Rvector'.
   Only for an 'lloffset_selection' and 'out_Rvector' of length <= INT_MAX.
   Don't use on an 'lloffset_selection' or 'out_Rvector' of length > INT_MAX! */
void _copy_Rvector_elts_from_selected_lloffsets(SEXP Rvector,
		const long long *lloffsets, const int *lloffset_selection,
		SEXP out_Rvector)
{
	SEXPTYPE Rtype = TYPEOF(Rvector);
	int out_len = LENGTH(out_Rvector);  /* assumed to be also the length
					       of the selection */
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		copy_int_elts_from_selected_lloffsets(INTEGER(Rvector),
				lloffsets, lloffset_selection, out_len,
				INTEGER(out_Rvector));
		return;
	    case REALSXP:
		copy_double_elts_from_selected_lloffsets(REAL(Rvector),
				lloffsets, lloffset_selection, out_len,
				REAL(out_Rvector));
		return;
	    case CPLXSXP:
		copy_Rcomplex_elts_from_selected_lloffsets(COMPLEX(Rvector),
				lloffsets, lloffset_selection, out_len,
				COMPLEX(out_Rvector));
		return;
	    case RAWSXP:
		copy_Rbyte_elts_from_selected_lloffsets(RAW(Rvector),
				lloffsets, lloffset_selection, out_len,
				RAW(out_Rvector));
		return;
	    case STRSXP:
		copy_character_elts_from_selected_lloffsets(Rvector,
				lloffsets, lloffset_selection, out_len,
				out_Rvector);
		return;
	    case VECSXP:
		copy_list_elts_from_selected_lloffsets(Rvector,
				lloffsets, lloffset_selection, out_len,
				out_Rvector);
		return;
	}
	error("SparseArray internal error in "
	      "_copy_Rvector_elts_from_selected_lloffsets():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
}


/****************************************************************************
 * _unary_minus_Rvector()
 */

static void unary_minus_int(const int *x, SEXP out_Rvector)
{
	R_xlen_t n = XLENGTH(out_Rvector);
	SEXPTYPE out_Rtype = TYPEOF(out_Rvector);
	int val;
	switch (out_Rtype) {
	    case INTSXP: {
		int *out = INTEGER(out_Rvector);
		for (R_xlen_t i = 0; i < n; i++) {
			val = x[i];
			if (val != NA_INTEGER)
				val = - val;
			out[i] = val;
		}
		return;
	    }
	    case REALSXP: {
		double *out = REAL(out_Rvector);
		for (R_xlen_t i = 0; i < n; i++) {
			val = x[i];
			if (val == NA_INTEGER) {
				out[i] = NA_REAL;
			} else {
				val = - val;
				out[i] = (double) val;
			}
		}
		return;
	    }
	    case CPLXSXP: {
		Rcomplex *out = COMPLEX(out_Rvector);
		for (R_xlen_t i = 0; i < n; i++, out++) {
			val = x[i];
			if (val == NA_INTEGER) {
				out->r = NA_REAL;
			} else {
				val = - val;
				out->r = (double) val;
			}
			out->i = double0;
		}
		return;
	    }
	}
	error("SparseArray internal error in unary_minus_int():\n"
	      "    output type \"%s\" is not supported", type2char(out_Rtype));
}

static void unary_minus_double(const double *x, SEXP out_Rvector)
{
	R_xlen_t n = XLENGTH(out_Rvector);
	SEXPTYPE out_Rtype = TYPEOF(out_Rvector);
	switch (out_Rtype) {
	    case REALSXP: {
		double *out = REAL(out_Rvector);
		for (R_xlen_t i = 0; i < n; i++)
			out[i] = - x[i];
		return;
	    }
	    case CPLXSXP: {
		Rcomplex *out = COMPLEX(out_Rvector);
		for (R_xlen_t i = 0; i < n; i++, out++) {
			out->r = - x[i];
			out->i = double0;
		}
		return;
	    }
	}
	error("SparseArray internal error in unary_minus_double():\n"
	      "    output type \"%s\" is not supported", type2char(out_Rtype));
}

static void unary_minus_Rcomplex(const Rcomplex *x, SEXP out_Rvector)
{
	R_xlen_t n = XLENGTH(out_Rvector);
	SEXPTYPE out_Rtype = TYPEOF(out_Rvector);
	if (out_Rtype == CPLXSXP) {
		Rcomplex *out = COMPLEX(out_Rvector);
		for (R_xlen_t i = 0; i < n; i++, x++, out++) {
			out->r = - x->r;
			out->i = - x->i;
		}
		return;
	}
	error("SparseArray internal error in unary_minus_Rcomplex():\n"
	      "    output type \"%s\" is not supported", type2char(out_Rtype));
}

/* 'out_Rvector' must be already allocated to the length of 'Rvector' and its
   type must be 'TYPEOF(Rvector)' or a bigger type.
   Use '_unary_minus_Rvector(x, x)' for **in-place** replacement. */
void _unary_minus_Rvector(SEXP Rvector, SEXP out_Rvector)
{
	R_xlen_t in_len = XLENGTH(Rvector);
	if (XLENGTH(out_Rvector) != in_len)
		error("SparseArray internal error in _unary_minus_Rvector():\n"
		      "    XLENGTH(out_Rvector) != in_len");
	SEXPTYPE Rtype = TYPEOF(Rvector);
	switch (Rtype) {
	    case INTSXP:
		unary_minus_int(INTEGER(Rvector), out_Rvector);
		return;
	    case REALSXP:
		unary_minus_double(REAL(Rvector), out_Rvector);
		return;
	    case CPLXSXP:
		unary_minus_Rcomplex(COMPLEX(Rvector), out_Rvector);
		return;
	}
	error("SparseArray internal error in _unary_minus_Rvector():\n"
	      "    input type \"%s\" is not supported", type2char(Rtype));
        return;
}

