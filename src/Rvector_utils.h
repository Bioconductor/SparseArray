#ifndef _RVECTOR_UTILS_H_
#define _RVECTOR_UTILS_H_

#include <Rdefines.h>

static const Rbyte Rbyte0 = 0;
static const int int0 = 0;
static const double double0 = 0.0;
/* Some old versions of gcc choke on this:
static const Rcomplex Rcomplex0 = {{double0, double0}}; */
static const Rcomplex Rcomplex0 = {{0.0, 0.0}};

static const Rbyte Rbyte1 = 1;
static const int int1 = 1;
static const double double1 = 1.0;
/* Some old versions of gcc choke on this:
static const Rcomplex Rcomplex1 = {{double1, double0}}; */
static const Rcomplex Rcomplex1 = {{1.0, 0.0}};

#define	IS_EMPTY_CHARSXP(x) ((x) != NA_STRING && isBlankString(CHAR(x)))

#define RCOMPLEX_IS_NA_OR_NaN(z) (ISNAN((z)->r) || ISNAN((z)->i))

typedef void (*CopyRVectorElt_FUNType)(
	SEXP in,  R_xlen_t in_offset,
	SEXP out, R_xlen_t out_offset);


/****************************************************************************
 * Inline functions
 */

/* Should be the same as doing (char *) (x) + sizeof(type) * (offset). */
#define SHIFT_DATAPTR(type, x, offset) (type *) (x) + (offset)

/* Restricted to types "logical", "integer", "double", "complex", and "raw".
   Should be the same as doing (char *) x + _get_Rtype_size(Rtype) * offset. */
static inline void *shift_dataptr(SEXPTYPE Rtype, void *x, R_xlen_t offset)
{
	switch (Rtype) {
	    case INTSXP: case LGLSXP: return SHIFT_DATAPTR(int, x, offset);
	    case REALSXP:             return SHIFT_DATAPTR(double, x, offset);
	    case CPLXSXP:             return SHIFT_DATAPTR(Rcomplex, x, offset);
	    case RAWSXP:              return SHIFT_DATAPTR(Rbyte, x, offset);
	}
	error("SparseArray internal error in shift_dataptr():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return NULL;  /* will never reach this */
}

static inline void copy_INTEGER_elt(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset)
{
	INTEGER(out)[out_offset] =
		in == R_NilValue ? int1 : INTEGER(in)[in_offset];
	return;
}

static inline void copy_NUMERIC_elt(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset)
{
	REAL(out)[out_offset] =
		in == R_NilValue ? double1 : REAL(in)[in_offset];
	return;
}

static inline void copy_COMPLEX_elt(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset)
{
	COMPLEX(out)[out_offset] =
		in == R_NilValue ? Rcomplex1 : COMPLEX(in)[in_offset];
	return;
}

static inline void copy_RAW_elt(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset)
{
	RAW(out)[out_offset] =
		in == R_NilValue ? Rbyte1 : RAW(in)[in_offset];
	return;
}

static inline void copy_CHARACTER_elt(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset)
{
	if (in == R_NilValue)
		error("SparseArray internal error in copy_CHARACTER_elt():\n"
		      "    lacunar leaf found in an SVT_SparseArray object "
		      "of type \"character\"");
	SET_STRING_ELT(out, out_offset, STRING_ELT(in, in_offset));
	return;
}

static inline void copy_LIST_elt(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset)
{
	if (in == R_NilValue)
		error("SparseArray internal error in copy_LIST_elt():\n"
		      "    lacunar leaf found in an SVT_SparseArray object "
		      "of type \"list\"");
	SET_VECTOR_ELT(out, out_offset, VECTOR_ELT(in, in_offset));
	return;
}


/****************************************************************************
 * Function prototypes
 */

SEXPTYPE _get_Rtype_from_Rstring(SEXP type);

size_t _get_Rtype_size(SEXPTYPE Rtype);

void _set_elts_to_val(
	SEXPTYPE Rtype,
	void *x,
	R_xlen_t offset,
	R_xlen_t n,
	const void *val
);

void _set_elts_to_zero(
	SEXPTYPE Rtype,
	void *x,
	R_xlen_t offset,
	R_xlen_t n
);

void _set_elts_to_one(
	SEXPTYPE Rtype,
	void *x,
	R_xlen_t offset,
	R_xlen_t n
);

void _set_elts_to_minus_one(
	SEXPTYPE Rtype,
	void *x,
	R_xlen_t offset,
	R_xlen_t n
);

void _set_elts_to_NA(
	SEXPTYPE Rtype,
	void *x,
	R_xlen_t offset,
	R_xlen_t n
);

void _set_Rsubvec_elts_to_val(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	R_xlen_t subvec_len,
	const void *val
);

void _set_Rsubvec_elts_to_zero(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	R_xlen_t subvec_len
);

void _set_Rsubvec_elts_to_one(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	R_xlen_t subvec_len
);

void _set_Rsubvec_elts_to_minus_one(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	R_xlen_t subvec_len
);

void _set_Rsubvec_elts_to_NA(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	R_xlen_t subvec_len
);

void _set_Rvector_elts_to_val(SEXP Rvector, const void *val);

void _set_Rvector_elts_to_zero(SEXP Rvector);

void _set_Rvector_elts_to_one(SEXP Rvector);

void _set_Rvector_elts_to_minus_one(SEXP Rvector);

void _set_Rvector_elts_to_NA(SEXP Rvector);

void _set_selected_elts_to_zero(
	SEXPTYPE Rtype,
	void *x,
	R_xlen_t offset,
	const int *selection,
	int n
);

void _set_selected_elts_to_one(
	SEXPTYPE Rtype,
	void *x,
	R_xlen_t offset,
	const int *selection,
	int n
);

void _set_selected_Rsubvec_elts_to_zero(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	const int *selection,
	int n
);

void _set_selected_Rsubvec_elts_to_one(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	const int *selection,
	int n
);

SEXP _new_Rvector0(
	SEXPTYPE Rtype,
	R_xlen_t len
);

SEXP _new_Rmatrix0(
	SEXPTYPE Rtype,
	int nrow,
	int ncol,
	SEXP dimnames
);

SEXP _new_Rarray0(
	SEXPTYPE Rtype,
	SEXP dim,
	SEXP dimnames
);

SEXP _new_Rvector1(
	SEXPTYPE Rtype,
	int len
);

SEXP _new_RarrayNA(
	SEXPTYPE Rtype,
	SEXP dim,
	SEXP dimnames
);

int _collect_offsets_of_nonzero_Rsubvec_elts(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	int subvec_len,
	int *out
);

int _collect_offsets_of_nonNA_Rsubvec_elts(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	int subvec_len,
	int *out
);

int _all_elts_equal_one(
	SEXPTYPE Rtype,
	const void *x,
	int n
);

int _all_Rsubvec_elts_equal_one(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	int subvec_len
);

int _all_Rvector_elts_equal_one(SEXP Rvector);

int _all_selected_Rsubvec_elts_equal_one(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	const int *selection,
	int n
);

CopyRVectorElt_FUNType _select_copy_Rvector_elt_FUN(SEXPTYPE Rtype);

void _copy_Rvector_elts(
	SEXP in,
	R_xlen_t in_offset,
	SEXP out,
	R_xlen_t out_offset,
	R_xlen_t nelt
);

void _copy_selected_int_elts(
	const int *in,
	const int *selection,
	int n,
	int *out
);

void _copy_selected_double_elts(
	const double *in,
	const int *selection,
	int n,
	double *out
);

void _copy_selected_Rcomplex_elts(
	const Rcomplex *in,
	const int *selection,
	int n,
	Rcomplex *out
);

void _copy_selected_Rbyte_elts(
	const Rbyte *in,
	const int *selection,
	int n,
	Rbyte *out
);

void _copy_selected_character_elts(
	SEXP in,
	R_xlen_t in_offset,
	const int *selection,
	int n,
	SEXP out
);

void _copy_selected_list_elts(
	SEXP in,
	R_xlen_t in_offset,
	const int *selection,
	int n,
	SEXP out
);

void _copy_selected_Rsubvec_elts(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	const int *selection,
	SEXP out_Rvector
);

SEXP _subset_Rsubvec(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	const int *selection,
	int n
);

void _copy_int_elts_to_offsets(
	const int *in,
	const int *selection,
	int n,
	int *out
);

void _copy_double_elts_to_offsets(
	const double *in,
	const int *selection,
	int n,
	double *out
);

void _copy_Rcomplex_elts_to_offsets(
	const Rcomplex *in,
	const int *selection,
	int n,
	Rcomplex *out
);

void _copy_Rbyte_elts_to_offsets(
	const Rbyte *in,
	const int *selection,
	int n,
	Rbyte *out
);

void _copy_character_elts_to_offsets(
	SEXP in,
	const int *selection,
	int n,
	SEXP out,
	R_xlen_t out_offset
);

void _copy_list_elts_to_offsets(
	SEXP in,
	const int *selection,
	int n,
	SEXP out,
	R_xlen_t out_offset
);

void _copy_Rvector_elts_to_offsets(
	SEXP in_Rvector,
	const int *selection,
	SEXP out_Rvector,
	R_xlen_t out_offset
);

void _copy_Rvector_elts_from_selected_offsets(
	SEXP in_Rvector,
	const int *offsets,
	const int *offset_selection,
	SEXP out_Rvector
);

void _copy_Rvector_elts_from_selected_lloffsets(
	SEXP in_Rvector,
	const long long *lloffsets,
	const int *lloffset_selection,
	SEXP out_Rvector
);

void _unary_minus_Rvector(
	SEXP in_Rvector,
	SEXP out_Rvector
);

#endif  /* _RVECTOR_UTILS_H_ */

