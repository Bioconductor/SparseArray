/****************************************************************************
 *                 Basic manipulation of SparseVec structs                  *
 ****************************************************************************/
#include "SparseVec.h"

#include "Rvector_utils.h"


/****************************************************************************
 * _alloc_buf_SparseVec()
 */

/* IMPORTANT: The caller must immediately call 'PROTECT(sv.nzvals)' on
   the returned SparseVec struct when 'Rtype' is STRSXP or VECSXP. */
SparseVec _alloc_buf_SparseVec(SEXPTYPE Rtype, int len, int na_background)
{
	if (na_background && (Rtype == RAWSXP || Rtype == VECSXP))
		error("SparseArray internal error in "
		      "_alloc_buf_SparseVec():\n    NaArray objects "
		      "of type \"%s\" are not supported", type2char(Rtype));
	SparseVec sv;
	sv.Rtype = Rtype;
	if (IS_STRSXP_OR_VECSXP(Rtype)) {
		sv.nzvals = PROTECT(allocVector(Rtype, (R_xlen_t) len));
	} else {
		size_t Rtype_size = _get_Rtype_size(Rtype);
		if (Rtype_size == 0)
			error("SparseArray internal error in "
			      "_alloc_buf_SparseVec():\n    type \"%s\" is "
			      "not supported", type2char(Rtype));
		sv.nzvals = R_alloc(len, Rtype_size);
	}
	sv.nzoffs = (int *) R_alloc(len, sizeof(int));
	sv.nzcount = 0;
	sv.len = len;
	sv.na_background = na_background;
	if (IS_STRSXP_OR_VECSXP(Rtype))
		UNPROTECT(1);
	return sv;
}


/****************************************************************************
 * _write_Rvector_block_to_SV()
 *
 * TODO: Maybe _make_leaf_from_Rvector_block() and
 * _make_naleaf_from_Rvector_block() should be based on this?
 */

static SEXPTYPE get_SV_write_Rtype(SEXP Rvector,
		const int *out_offs, int n, const SparseVec *out_sv)
{
	if (Rvector == R_NilValue)
		error("SparseArray internal error in "
		      "get_SV_write_Rtype():\n"
		      "    'Rvector' cannot be NULL");
	SEXPTYPE Rtype = get_SV_Rtype(out_sv);
	if (TYPEOF(Rvector) != Rtype)
		error("SparseArray internal error in "
		      "get_SV_write_Rtype():\n"
		      "    'Rvector' and 'out_sv' don't have the same type");
	if (out_offs == NULL) {
		if (n != out_sv->len)
			error("SparseArray internal error in "
			      "get_SV_write_Rtype():\n"
			      "    'out_offs == NULL' and 'n != out_sv->len'");
	} else {
		if (n > out_sv->len)
			error("SparseArray internal error in "
			      "get_SV_write_Rtype():\n"
			      "    'out_offs != NULL' and 'n > out_sv->len'");
	}
	return Rtype;
}

/* The only reason that we define 'RbyteNA' is to make
   'DEFINE_write_Rvector_block_to_typeSV_FUN(Rbyte)' and
   'DEFINE_write_Rvector_subset_to_typeSV_FUN(Rbyte)' work.
   Note that it absolutely doesn't matter what value we set 'RbyteNA'
   to because we don't support NaArray objects of type raw.
   This means that if SparseVec 'out_sv' is of type raw
   then 'out_sv->na_background' is guaranteed to be FALSE.
   In other words, 'RbyteNA' will **never** be used! */

#define RbyteNA Rbyte0  /* exact value doesn't matter, see above */

#define DEFINE_write_Rvector_block_to_typeSV_FUN(type)			  \
static void write_Rvector_block_to_ ## type ## SV(			  \
		const type *vals, int cycle_len,			  \
		const int *out_offs, int n, SparseVec *out_sv)		  \
{									  \
	type *out_nzvals = (type *) out_sv->nzvals;			  \
	type bg_val = out_sv->na_background ? type ## NA : type ## 0;	  \
	out_sv->nzcount = 0;						  \
	for (int k = 0; k < n; k++) {					  \
		type out_val = vals[cycle_len ? k % cycle_len : k];	  \
		if (type ## _equal(out_val, bg_val))			  \
			continue;					  \
		int out_off = out_offs == NULL ? k : out_offs[k];	  \
		APPEND_TO_NZVALS_NZOFFS(out_val, out_off,		  \
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);	  \
	}								  \
	return;								  \
}

DEFINE_write_Rvector_block_to_typeSV_FUN(int)
DEFINE_write_Rvector_block_to_typeSV_FUN(double)
DEFINE_write_Rvector_block_to_typeSV_FUN(Rcomplex)
DEFINE_write_Rvector_block_to_typeSV_FUN(Rbyte)

static void write_Rvector_block_to_characterSV(
		SEXP Rvector, R_xlen_t block_offset, int cycle_len,
		const int *out_offs, int n, SparseVec *out_sv)
{
	SEXP out_nzvals = (SEXP) out_sv->nzvals;  /* STRSXP */
	out_sv->nzcount = 0;
	for (int k = 0; k < n; k++) {
		int i = cycle_len ? k % cycle_len : block_offset + k;
		SEXP out_val = STRING_ELT(Rvector, i);
		if (IS_BG_CHARSXP(out_val, out_sv->na_background))
			continue;
		int out_off = out_offs == NULL ? k : out_offs[k];
		SET_STRING_ELT(out_nzvals, out_sv->nzcount, out_val);
		out_sv->nzoffs[out_sv->nzcount] = out_off;
		out_sv->nzcount++;
	}
	return;
}

static void write_Rvector_block_to_listSV(
		SEXP Rvector, R_xlen_t block_offset, int cycle_len,
		const int *out_offs, int n, SparseVec *out_sv)
{
	SEXP out_nzvals = (SEXP) out_sv->nzvals;  /* VECSXP */
	out_sv->nzcount = 0;
	for (int k = 0; k < n; k++) {
		int i = cycle_len ? k % cycle_len : block_offset + k;
		SEXP out_val = VECTOR_ELT(Rvector, i);
		if (out_val == R_NilValue)
			continue;
		int out_off = out_offs == NULL ? k : out_offs[k];
		SET_VECTOR_ELT(out_nzvals, out_sv->nzcount, out_val);
		out_sv->nzoffs[out_sv->nzcount] = out_off;
		out_sv->nzcount++;
	}
	return;
}

/* Fills 'out_sv' with the nonzero values from the block of 'n' elements
   in 'Rvector' that starts at offset 'block_offset'.
   'out_offs' must be NULL or an array of 'n' offsets (non-negative integers)
   that are strictly sorted (in ascending order). The last offset in the
   array must be < 'out_sv->len'.
   About the length of 'Rvector':
   - If 'block_offset != 0' then its length must be >= 'block_offset + n'.
   - If 'block_offset == 0' then 'Rvector' can be of any length (except 0
     unless 'n' is also 0) and will be recycled if its length is < 'n'. */
void _write_Rvector_block_to_SV(SEXP Rvector, R_xlen_t block_offset,
		const int *out_offs, int n, SparseVec *out_sv)
{
	SEXPTYPE Rtype = get_SV_write_Rtype(Rvector, out_offs, n, out_sv);
	int Rvector_len = LENGTH(Rvector), cycle_len = 0;
	if (block_offset != 0) {
		if (Rvector_len < block_offset + n)
			error("SparseArray internal error in "
			      "_write_Rvector_block_to_SV():\n"
			      "    LENGTH(Rvector) < block_offset + n");
	} else if (Rvector_len < n) {
		if (Rvector_len == 0)
			error("SparseArray internal error in "
			      "_write_Rvector_block_to_SV():\n"
			      "    LENGTH(Rvector) == 0");
		cycle_len = Rvector_len;
	}
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		write_Rvector_block_to_intSV(
			INTEGER(Rvector) + block_offset, cycle_len,
			out_offs, n, out_sv);
		return;
	    case REALSXP:
		write_Rvector_block_to_doubleSV(
			REAL(Rvector) + block_offset, cycle_len,
			out_offs, n, out_sv);
		return;
	    case CPLXSXP:
		write_Rvector_block_to_RcomplexSV(
			COMPLEX(Rvector) + block_offset, cycle_len,
			out_offs, n, out_sv);
		return;
	    case RAWSXP:
		write_Rvector_block_to_RbyteSV(
			RAW(Rvector) + block_offset, cycle_len,
			out_offs, n, out_sv);
		return;
	    case STRSXP:
		write_Rvector_block_to_characterSV(
			Rvector, block_offset, cycle_len,
			out_offs, n, out_sv);
		return;
	    case VECSXP:
		write_Rvector_block_to_listSV(
			Rvector, block_offset, cycle_len,
			out_offs, n, out_sv);
		return;
	}
	error("SparseArray internal error in "
	      "_write_Rvector_block_to_SV():\n"
	      "    'out_sv' of type \"%s\" not supported", type2char(Rtype));
	return;  /* will never reach this */
}


/****************************************************************************
 * _write_Rvector_subset_to_SV()
 */

#define DEFINE_write_Rvector_subset_to_typeSV_FUN(type)			  \
static void write_Rvector_subset_to_ ## type ## SV(			  \
		const type *vals, const int *selection,			  \
		const int *out_offs, int n, SparseVec *out_sv)		  \
{									  \
	type *out_nzvals = (type *) out_sv->nzvals;			  \
	type bg_val = out_sv->na_background ? type ## NA : type ## 0;	  \
	out_sv->nzcount = 0;						  \
	for (int k = 0; k < n; k++) {					  \
		type out_val = vals[selection[k]];			  \
		if (type ## _equal(out_val, bg_val))			  \
			continue;					  \
		int out_off = out_offs == NULL ? k : out_offs[k];	  \
		APPEND_TO_NZVALS_NZOFFS(out_val, out_off,		  \
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);	  \
	}								  \
	return;								  \
}

DEFINE_write_Rvector_subset_to_typeSV_FUN(int)
DEFINE_write_Rvector_subset_to_typeSV_FUN(double)
DEFINE_write_Rvector_subset_to_typeSV_FUN(Rcomplex)
DEFINE_write_Rvector_subset_to_typeSV_FUN(Rbyte)

static void write_Rvector_subset_to_characterSV(
		SEXP Rvector, const int *selection,
		const int *out_offs, int n, SparseVec *out_sv)
{
	SEXP out_nzvals = (SEXP) out_sv->nzvals;  /* STRSXP */
	out_sv->nzcount = 0;
	for (int k = 0; k < n; k++) {
		SEXP out_val = STRING_ELT(Rvector, selection[k]);
		if (IS_BG_CHARSXP(out_val, out_sv->na_background))
			continue;
		int out_off = out_offs == NULL ? k : out_offs[k];
		SET_STRING_ELT(out_nzvals, out_sv->nzcount, out_val);
		out_sv->nzoffs[out_sv->nzcount] = out_off;
		out_sv->nzcount++;
	}
	return;
}

static void write_Rvector_subset_to_listSV(
		SEXP Rvector, const int *selection,
		const int *out_offs, int n, SparseVec *out_sv)
{
	SEXP out_nzvals = (SEXP) out_sv->nzvals;  /* VECSXP */
	out_sv->nzcount = 0;
	for (int k = 0; k < n; k++) {
		SEXP out_val = VECTOR_ELT(Rvector, selection[k]);
		if (out_val == R_NilValue)
			continue;
		int out_off = out_offs == NULL ? k : out_offs[k];
		SET_VECTOR_ELT(out_nzvals, out_sv->nzcount, out_val);
		out_sv->nzoffs[out_sv->nzcount] = out_off;
		out_sv->nzcount++;
	}
	return;
}

void _write_Rvector_subset_to_SV(SEXP Rvector, const int *selection,
		const int *out_offs, int n, SparseVec *out_sv)
{
	SEXPTYPE Rtype = get_SV_write_Rtype(Rvector, out_offs, n, out_sv);
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		write_Rvector_subset_to_intSV(INTEGER(Rvector), selection,
			out_offs, n, out_sv);
		return;
	    case REALSXP:
		write_Rvector_subset_to_doubleSV(REAL(Rvector), selection,
			out_offs, n, out_sv);
		return;
	    case CPLXSXP:
		write_Rvector_subset_to_RcomplexSV(COMPLEX(Rvector), selection,
			out_offs, n, out_sv);
		return;
	    case RAWSXP:
		write_Rvector_subset_to_RbyteSV(RAW(Rvector), selection,
			out_offs, n, out_sv);
		return;
	    case STRSXP:
		write_Rvector_subset_to_characterSV(Rvector, selection,
			out_offs, n, out_sv);
		return;
	    case VECSXP:
		write_Rvector_subset_to_listSV(Rvector, selection,
			out_offs, n, out_sv);
		return;
	}
	error("SparseArray internal error in "
	      "_write_Rvector_subset_to_SV():\n"
	      "    'out_sv' of type \"%s\" not supported", type2char(Rtype));
	return;  /* will never reach this */
}


/****************************************************************************
 * _expand_intSV()
 * _expand_doubleSV()
 */

void _expand_intSV(const SparseVec *sv, int *out, int set_background)
{
	if (set_background) {
		if (sv->na_background) {
			_set_elts_to_NA(INTSXP, out, 0, sv->len);
		} else {
			_set_elts_to_zero(INTSXP, out, 0, sv->len);
		}
	}
	const int *nzvals_p = get_intSV_nzvals_p(sv);
	if (nzvals_p == NULL) {  /* lacunar SparseVec */
		_set_selected_elts_to_one(INTSXP, out,
				sv->nzoffs, get_SV_nzcount(sv), 0);
	} else {  /* regular SparseVec */
		_copy_int_elts_to_offsets(nzvals_p,
				sv->nzoffs, get_SV_nzcount(sv), out);
	}
	return;
}

void _expand_doubleSV(const SparseVec *sv, double *out, int set_background)
{
	if (set_background) {
		if (sv->na_background) {
			_set_elts_to_NA(REALSXP, out, 0, sv->len);
		} else {
			_set_elts_to_zero(REALSXP, out, 0, sv->len);
		}
	}
	const double *nzvals_p = get_doubleSV_nzvals_p(sv);
	if (nzvals_p == NULL) {  /* lacunar SparseVec */
		_set_selected_elts_to_one(REALSXP, out,
				sv->nzoffs, get_SV_nzcount(sv), 0);
	} else {  /* regular SparseVec */
		_copy_double_elts_to_offsets(nzvals_p,
				sv->nzoffs, get_SV_nzcount(sv), out);
	}
	return;
}

