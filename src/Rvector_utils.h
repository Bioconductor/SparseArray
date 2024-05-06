#ifndef _RVECTOR_UTILS_H_
#define _RVECTOR_UTILS_H_

#include <Rdefines.h>
#include <string.h>  /* for memcpy() */


static const Rbyte Rbyte0 = 0;
static const int int0 = 0;
static const double double0 = 0.0;
static const Rcomplex Rcomplex0 = {{0.0, 0.0}};

static const Rbyte Rbyte1 = 1;
static const int int1 = 1;
static const double double1 = 1.0;
static const Rcomplex Rcomplex1 = {{1.0, 0.0}};


/****************************************************************************
 * typedefs
 */

typedef void (*CopyRVectorElt_FUNType)(
	SEXP in,  R_xlen_t in_offset,
	SEXP out, R_xlen_t out_offset);

typedef void (*CopyRVectorElts_FUNType)(
	SEXP in,  R_xlen_t in_offset,
	SEXP out, R_xlen_t out_offset,
	R_xlen_t nelt);


/****************************************************************************
 * Inline functions
 */

static inline void _copy_INTEGER_elt(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset)
{
	INTEGER(out)[out_offset] = INTEGER(in)[in_offset];
	return;
}

static inline void _copy_INTEGER_elts(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset,
		R_xlen_t nelt)
{
	void *dest, *src;

	dest = INTEGER(out) + out_offset;
	src  = INTEGER(in)  + in_offset;
	memcpy(dest, src, sizeof(int) * nelt);
	return;
}

static inline void _copy_NUMERIC_elt(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset)
{
	REAL(out)[out_offset] = REAL(in)[in_offset];
	return;
}

static inline void _copy_NUMERIC_elts(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset,
		R_xlen_t nelt)
{
	void *dest = REAL(out) + out_offset;
	void *src  = REAL(in)  + in_offset;
	memcpy(dest, src, sizeof(double) * nelt);
	return;
}

static inline void _copy_COMPLEX_elt(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset)
{
	COMPLEX(out)[out_offset] = COMPLEX(in)[in_offset];
	return;
}

static inline void _copy_COMPLEX_elts(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset,
		R_xlen_t nelt)
{
	void *dest = COMPLEX(out) + out_offset;
	void *src  = COMPLEX(in)  + in_offset;
	memcpy(dest, src, sizeof(Rcomplex) * nelt);
	return;
}

static inline void _copy_RAW_elt(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset)
{
	RAW(out)[out_offset] = RAW(in)[in_offset];
	return;
}

static inline void _copy_RAW_elts(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset,
		R_xlen_t nelt)
{
	void *dest = RAW(out) + out_offset;
	void *src  = RAW(in)  + in_offset;
	memcpy(dest, src, sizeof(Rbyte) * nelt);
	return;
}

static inline void _copy_CHARACTER_elt(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset)
{
	SET_STRING_ELT(out, out_offset, STRING_ELT(in, in_offset));
	return;
}

static inline void _copy_CHARACTER_elts(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset,
		R_xlen_t nelt)
{
	for (R_xlen_t k = 0; k < nelt; k++)
		_copy_CHARACTER_elt(in, in_offset + k, out, out_offset + k);
	return;
}

static inline void _copy_LIST_elt(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset)
{
	SET_VECTOR_ELT(out, out_offset, VECTOR_ELT(in, in_offset));
	return;
}

static inline void _copy_LIST_elts(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset,
		R_xlen_t nelt)
{
	for (R_xlen_t k = 0; k < nelt; k++)
		_copy_LIST_elt(in, in_offset + k, out, out_offset + k);
	return;
}


/****************************************************************************
 * Function prototypes
 */

SEXPTYPE _get_Rtype_from_Rstring(SEXP type);

size_t _get_Rtype_size(SEXPTYPE Rtype);

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

void _set_Rsubvec_to_one(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	int subvec_len
);

SEXP _new_Rvector1(
	SEXPTYPE Rtype,
	int len
);

CopyRVectorElt_FUNType _select_copy_Rvector_elt_FUN(SEXPTYPE Rtype);

CopyRVectorElts_FUNType _select_copy_Rvector_elts_FUN(SEXPTYPE Rtype);

int _collect_offsets_of_nonzero_Rsubvec_elts(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	int subvec_len,
	int *out
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

int _all_Rsubvec_elts_equal_one(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	int subvec_len
);

int _all_selected_Rsubvec_elts_equal_one(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	const int *selection,
	int n
);

void _copy_selected_ints(
	const int *in,
	const int *selection,
	int n,
	int *out
);

void _copy_selected_doubles(
	const double *in,
	const int *selection,
	int n,
	double *out
);

void _copy_selected_Rcomplexes(
	const Rcomplex *in,
	const int *selection,
	int n,
	Rcomplex *out
);

void _copy_selected_Rbytes(
	const Rbyte *in,
	const int *selection,
	int n,
	Rbyte *out
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

void _copy_ints_to_offsets(
	const int *in,
	const int *selection,
	int n,
	int *out
);

void _copy_doubles_to_offsets(
	const double *in,
	const int *selection,
	int n,
	double *out
);

void _copy_Rcomplexes_to_offsets(
	const Rcomplex *in,
	const int *selection,
	int n,
	Rcomplex *out
);

void _copy_Rbytes_to_offsets(
	const Rbyte *in,
	const int *selection,
	int n,
	Rbyte *out
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

const char *_unary_minus_Rvector(
	SEXP in_Rvector,
	SEXP out_Rvector
);

#endif  /* _RVECTOR_UTILS_H_ */

