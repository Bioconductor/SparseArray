#ifndef _RVECTOR_UTILS_H_
#define _RVECTOR_UTILS_H_

#include <Rdefines.h>

typedef void (*CopyRVectorElt_FUNType)(
	SEXP in,  R_xlen_t in_offset,
	SEXP out, R_xlen_t out_offset);

typedef void (*CopyRVectorElts_FUNType)(
	SEXP in,  R_xlen_t in_offset,
	SEXP out, R_xlen_t out_offset,
	R_xlen_t nelt);

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
	void *dest, *src;

	dest = REAL(out) + out_offset;
	src  = REAL(in)  + in_offset;
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
	void *dest, *src;

	dest = COMPLEX(out) + out_offset;
	src  = COMPLEX(in)  + in_offset;
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
	void *dest, *src;

	dest = RAW(out) + out_offset;
	src  = RAW(in)  + in_offset;
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
	R_xlen_t k;

	for (k = 0; k < nelt; k++)
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
	R_xlen_t k;

	for (k = 0; k < nelt; k++)
		_copy_LIST_elt(in, in_offset + k, out, out_offset + k);
	return;
}

SEXPTYPE _get_Rtype_from_Rstring(SEXP type);

SEXP _new_Rvector(
	SEXPTYPE Rtype,
	R_xlen_t len
);

SEXP _new_Rarray(
	SEXPTYPE Rtype,
	SEXP dim,
	SEXP dimnames
);

CopyRVectorElt_FUNType _select_copy_Rvector_elt_FUN(SEXPTYPE Rtype);

CopyRVectorElts_FUNType _select_copy_Rvector_elts_FUN(SEXPTYPE Rtype);

int _collect_offsets_of_nonzero_Rsubvec_elts(
	SEXP Rvector,
	R_xlen_t subvec_offset,
	int subvec_len,
	int *offs_buf
);

void _reset_selected_Rvector_elts(
	SEXP Rvector,
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
	SEXP in_Rvector,
	R_xlen_t in_offset,
	const int *selection,
	SEXP out_Rvector
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

#endif  /* _RVECTOR_UTILS_H_ */

