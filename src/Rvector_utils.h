#ifndef _RVECTOR_UTILS_H_
#define _RVECTOR_UTILS_H_

#include <Rdefines.h>
#include <Rversion.h>
#if R_VERSION < R_Version(4, 6, 0)
# define DATAPTR_RW(x) DATAPTR(x)
#endif


/* Note that R does not define NA_INTEGER, NA_REAL, or NA_STRING as const
   variables so we can't define intNA, doubleNA, RcomplexNA, or characterNA
   as const variable either. Instead we need to define them as regular (i.e.
   non-const) variables and initialize them in R_init_SparseArray(). */

static const int int0 = 0, int1 = 1;
extern int intNA;            /* initialized in R_init_SparseArray() */

static const double double0 = 0.0, double1 = 1.0;
extern double doubleNA;      /* initialized in R_init_SparseArray() */

/* Some old versions of gcc choke on this:
static const Rcomplex Rcomplex0 = {{double0, double0}}, {{double1, double0}}; */
static const Rcomplex Rcomplex0 = {{0.0, 0.0}}, Rcomplex1 = {{1.0, 0.0}};
extern Rcomplex RcomplexNA;  /* initialized in R_init_SparseArray() */

static const Rbyte Rbyte0 = 0, Rbyte1 = 1;

/* Initialized in R_init_SparseArray() */
extern SEXP character0, character1, characterNA;

extern SEXP list0;           /* initialized in R_init_SparseArray() */


#define IS_STRSXP_OR_VECSXP(Rtype) ((Rtype) == STRSXP || (Rtype) == VECSXP)

#define RCOMPLEX_IS_NA_OR_NaN(z) (ISNAN((z)->r) || ISNAN((z)->i))

typedef void (*CopyRVectorEltFUN)(
	SEXP in,  R_xlen_t in_offset,
	SEXP out, R_xlen_t out_offset);


/****************************************************************************
 * Inline functions
 */

static inline int is_int0(int x)
{
	return x == int0;
}
static inline int is_int1(int x)
{
	return x == int1;
}
static inline int is_intNA(int x)
{
	return x == intNA;
}

static inline int is_double0(double x)
{
	return x == double0;
}
static inline int is_double1(double x)
{
	return x == double1;
}
static inline int is_doubleNA(double x)
{
	return R_IsNA(x);  // do NOT use ISNAN()!
}

static inline int is_Rcomplex0(Rcomplex x)
{
	return x.r == double0 && x.i == double0;
}
static inline int is_Rcomplex1(Rcomplex x)
{
	return x.r == double1 && x.i == double0;
}
static inline int is_RcomplexNA(Rcomplex x)
{
	return R_IsNA(x.r) || R_IsNA(x.i);  // do NOT use ISNAN()!
}

static inline int is_Rbyte0(Rbyte x)
{
	return x == Rbyte0;
}
static inline int is_Rbyte1(Rbyte x)
{
	return x == Rbyte1;
}

/* 'x' must be a CHARSXP.
   Note that we're comparing the CHARSXPs' addresses, not their values.
   However, this is much faster, but also, and most importantly, it's
   equivalent to comparing their values. That's because CHARSXPs with the
   same value are expected to have the same address, thanks to R's global
   CHARSXP cache. */
static inline int is_character0(SEXP x)
{
	return x == character0;   /* comparing CHARSXPs' addresses */
}
static inline int is_character1(SEXP x)
{
	return x == character1;   /* comparing CHARSXPs' addresses */
}
static inline int is_characterNA(SEXP x)
{
	return x == characterNA;  /* comparing CHARSXPs' addresses */
}

static inline int is_list0(SEXP x)
{
	return x == list0;
}


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

void _fill_Rvector_block_with_val(
	SEXP Rvector,
	R_xlen_t block_offset,
	R_xlen_t block_len,
	const void *val
);

void _fill_Rvector_block_with_zeros(
	SEXP Rvector,
	R_xlen_t block_offset,
	R_xlen_t block_len
);

void _fill_Rvector_block_with_ones(
	SEXP Rvector,
	R_xlen_t block_offset,
	R_xlen_t block_len
);

void _fill_Rvector_block_with_minus_one(
	SEXP Rvector,
	R_xlen_t block_offset,
	R_xlen_t block_len
);

void _fill_Rvector_block_with_NA(
	SEXP Rvector,
	R_xlen_t block_offset,
	R_xlen_t block_len
);

void _fill_Rvector_with_val(SEXP Rvector, const void *val);

void _fill_Rvector_with_zeros(SEXP Rvector);

void _fill_Rvector_with_ones(SEXP Rvector);

void _fill_Rvector_with_minus_one(SEXP Rvector);

void _fill_Rvector_with_NA(SEXP Rvector);

void _set_selected_elts_to_zero(
	SEXPTYPE Rtype,
	void *x,
	const int *selection,
	int selection_len,
	R_xlen_t selection_offset
);

void _set_selected_elts_to_one(
	SEXPTYPE Rtype,
	void *x,
	const int *selection,
	int selection_len,
	R_xlen_t selection_offset
);

void _fill_Rvector_subset_with_zeros(
	SEXP Rvector,
	const int *selection,
	int selection_len,
	R_xlen_t selection_offset
);

void _fill_Rvector_subset_with_ones(
	SEXP Rvector,
	const int *selection,
	int selection_len,
	R_xlen_t selection_offset
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

SEXP _new_RvectorNA(
	SEXPTYPE Rtype,
	R_xlen_t len
);

SEXP _new_RarrayNA(
	SEXPTYPE Rtype,
	SEXP dim,
	SEXP dimnames
);

int _collect_offsets_of_nonzero_elts_in_Rvector_block(
	SEXP Rvector,
	R_xlen_t block_offset,
	int block_len,
	int *out
);

int _collect_offsets_of_nonNA_elts_in_Rvector_block(
	SEXP Rvector,
	R_xlen_t block_offset,
	int block_len,
	int *out
);

int _all_elts_equal_one(
	SEXPTYPE Rtype,
	const void *x,
	int n
);

int _Rvector_block_is_filled_with_ones(
	SEXP Rvector,
	R_xlen_t block_offset,
	int block_len
);

int _Rvector_is_filled_with_ones(SEXP Rvector);

int _Rvector_subset_is_filled_with_ones(
	SEXP Rvector,
	const int *selection,
	int selection_len,
	R_xlen_t selection_offset
);

CopyRVectorEltFUN _select_copy_Rvector_elt_FUN(SEXPTYPE Rtype);

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

void _copy_Rvector_subset(
	SEXP Rvector,
	const int *selection,
	R_xlen_t selection_offset,
	SEXP out_Rvector
);

SEXP _subset_Rvector(
	SEXP Rvector,
	const int *selection,
	int selection_len,
	R_xlen_t selection_offset
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

