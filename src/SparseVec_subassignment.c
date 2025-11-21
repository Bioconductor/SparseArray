/****************************************************************************
 *                     Subassignment of a sparse vector                     *
 ****************************************************************************/
#include "SparseVec_subassignment.h"


/* The only reason that we define 'RbyteNA' and 'listNA' is to make
   'DEFINE_subassign_typeSV_FUN(Rbyte)' and other macros defined in
   this file work. Note that these macros contain the following line:

       type bg_val = out_sv->na_background ? type ## NA : type ## 0;

   So it absolutely doesn't matter what value we set 'RbyteNA' and 'listNA'
   to because we don't support NaArray objects of type raw or list.
   This means that if SparseVec 'out_sv' is of type raw or list
   then 'out_sv->na_background' is guaranteed to be FALSE.
   In other words, 'RbyteNA' and 'listNA' will **never** be used! */

#define RbyteNA Rbyte0  /* exact value doesn't matter, see above */
#define listNA list0    /* exact value doesn't matter, see above */


static SEXPTYPE get_SV_subassign_Rtype(SEXPTYPE expected_Rtype,
		const SparseVec *sv, const SparseVec *out_sv)
{
	SEXPTYPE out_Rtype = get_SV_Rtype(out_sv);
	if (out_Rtype != expected_Rtype)
		error("SparseArray internal error in "
		      "get_SV_subassign_Rtype():\n"
		      "    'out_sv' does not have the expected type");
	if (sv->len != out_sv->len || get_SV_Rtype(sv) != out_Rtype)
		error("SparseArray internal error in "
		      "get_SV_subassign_Rtype():\n"
		      "    'sv' and 'out_sv' are incompatible");
	return expected_Rtype;
}


/****************************************************************************
 * Inline functions next_subassign_<type>SV_out_val()
 */

#define DEFINE_next_subassign_typeSV_out_val_FUN(type)			   \
static inline int next_subassign_ ## type ## SV_out_val(		   \
		const SparseVec *sv1, const int *offs, int n,		   \
		const type *vals2, int cycle_len, const int *selection2,   \
		int *k1, int *k, int *out_off, type *out_val)		   \
{									   \
	int ret = next_offset(sv1->nzoffs, get_SV_nzcount(sv1),		   \
			      offs, n, *k1, *k, out_off);		   \
	switch (ret) {							   \
	    case 1:							   \
		*out_val = get_ ## type ## SV_nzval(sv1, *k1);		   \
		(*k1)++;						   \
		break;							   \
	    case 3:							   \
		(*k1)++;						   \
	    case 2: {							   \
		int i2;							   \
		if (selection2 == NULL) {				   \
			i2 = cycle_len ? *k % cycle_len : *k;		   \
		} else {						   \
			i2 = selection2[*k];				   \
		}							   \
		*out_val = vals2[i2];					   \
		(*k)++;							   \
		break;							   \
	    }								   \
	}								   \
	return ret;							   \
}

DEFINE_next_subassign_typeSV_out_val_FUN(int)
DEFINE_next_subassign_typeSV_out_val_FUN(double)
DEFINE_next_subassign_typeSV_out_val_FUN(Rcomplex)
DEFINE_next_subassign_typeSV_out_val_FUN(Rbyte)

static inline int next_subassign_characterSV_out_val(
		const SparseVec *sv1, const int *offs, int n,
		SEXP Rvector, R_xlen_t block_offset, int cycle_len,
		const int *selection,
		int *k1, int *k, int *out_off, SEXP *out_val)
{
	int ret = next_offset(sv1->nzoffs, get_SV_nzcount(sv1),
			      offs, n, *k1, *k, out_off);
	switch (ret) {
	    case 1:
		*out_val = get_characterSV_nzval(sv1, *k1);
		(*k1)++;
		break;
	    case 3:
		(*k1)++;
	    case 2: {
		int i;
		if (selection == NULL) {
			i = block_offset + (cycle_len ? *k % cycle_len : *k);
		} else {
			i = selection[*k];
		}
		*out_val = STRING_ELT(Rvector, i);
		(*k)++;
		break;
	    }
	}
	return ret;
}

static inline int next_subassign_listSV_out_val(
		const SparseVec *sv1, const int *offs, int n,
		SEXP Rvector, R_xlen_t block_offset, int cycle_len,
		const int *selection,
		int *k1, int *k, int *out_off, SEXP *out_val)
{
	int ret = next_offset(sv1->nzoffs, get_SV_nzcount(sv1),
			      offs, n, *k1, *k, out_off);
	switch (ret) {
	    case 1:
		*out_val = get_listSV_nzval(sv1, *k1);
		(*k1)++;
		break;
	    case 3:
		(*k1)++;
	    case 2: {
		int i;
		if (selection == NULL) {
			i = block_offset + (cycle_len ? *k % cycle_len : *k);
		} else {
			i = selection[*k];
		}
		*out_val = VECTOR_ELT(Rvector, i);
		(*k)++;
		break;
	    }
	}
	return ret;
}


/****************************************************************************
 * Inline functions process_subassign_<type>SV_out_val()
 */

#define DEFINE_process_subassign_typeSV_out_val_FUN(type)		   \
static inline int process_subassign_ ## type ## SV_out_val(int ret,	   \
		type out_val, int out_off, SparseVec *out_sv,		   \
		type bg_val, const SparseVec *sv1, int k1)		   \
{									   \
	int is_effrep = 0;						   \
	if (ret == 2) {							   \
		if (type ## _equal(out_val, bg_val))			   \
			return 0;  /* zero replaces zero */		   \
		is_effrep = 1;     /* nonzero replaces zero */		   \
	} else if (ret == 3) {						   \
		if (type ## _equal(out_val, bg_val))			   \
			return 1;  /* zero replaces nonzero */		   \
		/* nonzero replaces nonzero */				   \
		type v1 = get_ ## type ## SV_nzval(sv1, k1 - 1);	   \
		if (!type ## _equal(v1, out_val))			   \
			is_effrep = 1;					   \
	}								   \
	APPEND_TO_NZVALS_NZOFFS(out_val, out_off,			   \
		(type *) out_sv->nzvals, out_sv->nzoffs, out_sv->nzcount); \
	return is_effrep;						   \
}

DEFINE_process_subassign_typeSV_out_val_FUN(int)
DEFINE_process_subassign_typeSV_out_val_FUN(double)
DEFINE_process_subassign_typeSV_out_val_FUN(Rcomplex)
DEFINE_process_subassign_typeSV_out_val_FUN(Rbyte)

/* Note that when comparing CHARSXPs 'v1' and 'out_val' below (v1 != out_val),
   we compare their **addresses**, not their **values**.
   However, this is much faster, but also, and most importantly, it's
   equivalent to comparing their values. That's because CHARSXPs with the
   same value are expected to have the same address, thanks to R's global
   CHARSXP cache.
   In any case, even if 'v1 != out_val' were to produce false positives,
   it would not be such a big deal because the main reason for counting the
   number of **effective** replacements (neffrep) is to avoid copying an SVT
   leaf when a subassignment does not modify it (i.e. when 'neffrep == 0').
   So in the worst case, these false positives simply mean that we would
   still copy a leaf touched by the subassignment operation, even if the
   leaf is not modified by the subassignment. */
static inline int process_subassign_characterSV_out_val(int ret,
		SEXP out_val, int out_off, SparseVec *out_sv,
		SEXP bg_val, const SparseVec *sv1, int k1)
{
	int is_effrep = 1;
	if (ret == 2) {
		if (out_val == bg_val)
			return 0;
	} else if (ret == 3) {
		if (out_val == bg_val)
			return 1;
		SEXP v1 = get_characterSV_nzval(sv1, k1 - 1);
		/* See note above about this comparison. */
		if (v1 == out_val)
			is_effrep = 0;
	}
	SET_STRING_ELT((SEXP) out_sv->nzvals, out_sv->nzcount, out_val);
	out_sv->nzoffs[out_sv->nzcount] = out_off;
	out_sv->nzcount++;
	return is_effrep;
}

/* Note that when comparing VECSXP elements 'v1' and 'out_val' below
   (v1 != out_val), we compare their **addresses**, not their **values**.
   Comparing the values would be too costly. So yes, 'v1 != out_val' can
   produce false positives, but it's not a big deal because the main
   reason for counting the number of **effective** replacements (neffrep)
   is to avoid copying an SVT leaf when a subassignment does not modify
   it (i.e. when 'neffrep == 0').
   So in the worst case, these false positives simply mean that we will
   still copy a leaf touched by the subassignment operation, even if the
   leaf is not modified by the subassignment.
   However, comparing the addresses will still do a good job in a situation
   like:

       value <- subset_Array_by_Nindex(svt1, Nindex)
       svt2 <- subassign_Array_by_Nindex(svt1, Nindex, as.array(value))

   where no copy will be triggered ('svt2@SVT' will have the same address
   as 'svt1@SVT'). */
static inline int process_subassign_listSV_out_val(int ret,
		SEXP out_val, int out_off, SparseVec *out_sv,
		SEXP bg_val, const SparseVec *sv1, int k1)
{
	int is_effrep = 1;
	if (ret == 2) {
		if (out_val == bg_val)
			return 0;
	} else if (ret == 3) {
		if (out_val == bg_val)
			return 1;
		SEXP v1 = get_listSV_nzval(sv1, k1 - 1);
		/* See note above about this comparison. */
		if (v1 == out_val)
			is_effrep = 0;
	}
	SET_VECTOR_ELT((SEXP) out_sv->nzvals, out_sv->nzcount, out_val);
	out_sv->nzoffs[out_sv->nzcount] = out_off;
	out_sv->nzcount++;
	return is_effrep;
}


/****************************************************************************
 * subassign_SV()
 */

#define DEFINE_subassign_typeSV_FUN(type)				   \
static int subassign_ ## type ## SV(					   \
		const SparseVec *sv1, const int *offs, int n,		   \
		const type *vals2, int cycle_len, const int *selection2,   \
		SparseVec *out_sv)					   \
{									   \
	type bg_val = out_sv->na_background ? type ## NA : type ## 0;	   \
	out_sv->nzcount = 0;						   \
	int neffrep = 0, ret, k1 = 0, k = 0, out_off;			   \
	type out_val;							   \
	while ((ret = next_subassign_ ## type ## SV_out_val(sv1, offs, n,  \
				vals2, cycle_len, selection2,		   \
				&k1, &k, &out_off, &out_val)))		   \
	{								   \
		neffrep += process_subassign_ ## type ## SV_out_val(ret,   \
				out_val, out_off, out_sv,		   \
				bg_val, sv1, k1);			   \
	}								   \
	return neffrep;							   \
}

DEFINE_subassign_typeSV_FUN(int)
DEFINE_subassign_typeSV_FUN(double)
DEFINE_subassign_typeSV_FUN(Rcomplex)
DEFINE_subassign_typeSV_FUN(Rbyte)

static int subassign_characterSV(
		const SparseVec *sv1, const int *offs, int n,
		SEXP Rvector, R_xlen_t block_offset, int cycle_len,
		const int *selection,
		SparseVec *out_sv)
{
	SEXP bg_val = out_sv->na_background ? characterNA : character0;
	out_sv->nzcount = 0;
	int neffrep = 0, ret, k1 = 0, k = 0, out_off;
	SEXP out_val;
	while ((ret = next_subassign_characterSV_out_val(sv1, offs, n,
				Rvector, block_offset, cycle_len, selection,
				&k1, &k, &out_off, &out_val)))
	{
		neffrep += process_subassign_characterSV_out_val(ret,
				out_val, out_off, out_sv,
				bg_val, sv1, k1);
	}
	return neffrep;
}

static int subassign_listSV(
		const SparseVec *sv1, const int *offs, int n,
		SEXP Rvector, R_xlen_t block_offset, int cycle_len,
		const int *selection,
		SparseVec *out_sv)
{
	out_sv->nzcount = 0;
	int neffrep = 0, ret, k1 = 0, k = 0, out_off;
	SEXP out_val;
	while ((ret = next_subassign_listSV_out_val(sv1, offs, n,
				Rvector, block_offset, cycle_len, selection,
				&k1, &k, &out_off, &out_val)))
	{
		neffrep += process_subassign_listSV_out_val(ret,
				out_val, out_off, out_sv,
				list0, sv1, k1);
	}
	return neffrep;
}

static int subassign_SV(
		const SparseVec *sv1, const int *offs, int n,
		SEXP Rvector, R_xlen_t block_offset, const int *selection,
		SparseVec *out_sv)
{
	if (Rvector == R_NilValue)
		error("SparseArray internal error in subassign_SV():\n"
		      "    'Rvector' cannot be NULL");
	SEXPTYPE Rtype = get_SV_subassign_Rtype(TYPEOF(Rvector), sv1, out_sv);
	int Rvector_len = LENGTH(Rvector), cycle_len = 0;
	if (block_offset != 0) {
		if (selection != NULL)
			error("SparseArray internal error in subassign_SV():\n"
			      "    block_offset != 0 && selection != NULL");
		if (Rvector_len < block_offset + n)
			error("SparseArray internal error in subassign_SV():\n"
			      "    LENGTH(Rvector) < block_offset + n");
	} else if (selection == NULL && Rvector_len < n) {
		if (Rvector_len == 0)
			error("SparseArray internal error in subassign_SV():\n"
			      "    LENGTH(Rvector) == 0");
		cycle_len = Rvector_len;
	}
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		return subassign_intSV(sv1, offs, n,
				INTEGER(Rvector) + block_offset, cycle_len,
				selection, out_sv);
	    case REALSXP:
		return subassign_doubleSV(sv1, offs, n,
				REAL(Rvector) + block_offset, cycle_len,
				selection, out_sv);
	    case CPLXSXP:
		return subassign_RcomplexSV(sv1, offs, n,
				COMPLEX(Rvector) + block_offset, cycle_len,
				selection, out_sv);
	    case RAWSXP:
		return subassign_RbyteSV(sv1, offs, n,
				RAW(Rvector) + block_offset, cycle_len,
				selection, out_sv);
	    case STRSXP:
		return subassign_characterSV(sv1, offs, n,
				Rvector, block_offset, cycle_len,
				selection, out_sv);
	    case VECSXP:
		return subassign_listSV(sv1, offs, n,
				Rvector, block_offset, cycle_len,
				selection, out_sv);
	}
	error("SparseArray internal error in subassign_SV():\n"
	      "    'out_sv' of type \"%s\" not supported", type2char(Rtype));
	return 0;  /* will never reach this */
}


/****************************************************************************
 * _subassign_SV_with_Rvector_block()
 * _subassign_SV_with_Rvector_subset()
 */

/* 'sv->len' and 'out_sv->len' must be the same. 'sv' can be lacunar.
   'offs' must be an array of 'n' offsets (non-negative integers) that
   are sorted in strictly ascending order. The last offset in the array
   must be < 'sv->len'.
   The replacement value (a.k.a. right value) is the block of 'n' elements
   in 'Rvector' that starts at offset 'block_offset'. It can contain zeros.
   About the length of 'Rvector':
   - If 'block_offset != 0' then its length must be >= 'block_offset + n'.
   - If 'block_offset == 0' then 'Rvector' can be of any length (except 0
     unless 'n' is also 0) and will be recycled if its length is < 'n'.
   Returns the number of **effective** replacements, that is, the number of
   offsets for which the subassignment operation effectively modifies the
   original value. */
int _subassign_SV_with_Rvector_block(
		const SparseVec *sv, const int *offs, int n,
		SEXP Rvector, R_xlen_t block_offset, SparseVec *out_sv)
{
	return subassign_SV(sv, offs, n,
			    Rvector, block_offset, NULL, out_sv);
}

/* Same as _subassign_SV_with_Rvector_block() above except that the
   replacement value is the subset of 'Rvector' obtained by extracting
   the selected elements. These are the 'n' elements at the positions
   indicated by 'selection', an array of zero-based indices into 'Rvector'. */
int _subassign_SV_with_Rvector_subset(
		const SparseVec *sv, const int *offs, int n,
		SEXP Rvector, const int *selection, SparseVec *out_sv)
{
	return subassign_SV(sv, offs, n, Rvector, 0, selection, out_sv);
}


/****************************************************************************
 * Inline functions next_subassign_<type>SV_with_SV_out_val()
 */

#define DEFINE_next_subassign_typeSV_with_SV_out_val_FUN(type, out_type)   \
static inline int next_subassign_ ## type ## SV_with_SV_out_val(	   \
		const SparseVec *sv1, const int *offs,			   \
		const SparseVec *sv2, out_type bg_val,			   \
		int *k1, int *k, int *k2,				   \
		int *out_off, out_type *out_val)			   \
{									   \
	int ret = next_offset(sv1->nzoffs, get_SV_nzcount(sv1),		   \
			      offs, sv2->len, *k1, *k, out_off);	   \
	switch (ret) {							   \
	    case 1:							   \
		*out_val = get_ ## type ## SV_nzval(sv1, *k1);		   \
		(*k1)++;						   \
		break;							   \
	    case 3:							   \
		(*k1)++;						   \
	    case 2:							   \
		if (*k2 < get_SV_nzcount(sv2) && sv2->nzoffs[*k2] == *k) { \
			*out_val = get_ ## type ## SV_nzval(sv2, *k2);	   \
			(*k2)++;					   \
		} else {						   \
			*out_val = bg_val;				   \
		}							   \
		(*k)++;							   \
		break;							   \
	}								   \
	return ret;							   \
}

DEFINE_next_subassign_typeSV_with_SV_out_val_FUN(int, int)
DEFINE_next_subassign_typeSV_with_SV_out_val_FUN(double, double)
DEFINE_next_subassign_typeSV_with_SV_out_val_FUN(Rcomplex, Rcomplex)
DEFINE_next_subassign_typeSV_with_SV_out_val_FUN(Rbyte, Rbyte)
DEFINE_next_subassign_typeSV_with_SV_out_val_FUN(character, SEXP)
DEFINE_next_subassign_typeSV_with_SV_out_val_FUN(list, SEXP)


/****************************************************************************
 * _subassign_SV_with_SV()
 */

#define DEFINE_subassign_typeSV_with_SV_FUN(type, out_type)		   \
static int subassign_ ## type ## SV_with_SV(				   \
		const SparseVec *sv1, const int *offs,			   \
		const SparseVec *sv2, SparseVec *out_sv)		   \
{									   \
	out_type bg_val = out_sv->na_background ? type ## NA : type ## 0;  \
	out_sv->nzcount = 0;						   \
	int neffrep = 0, ret, k1 = 0, k = 0, k2 = 0, out_off;		   \
	out_type out_val;						   \
	while ((ret = next_subassign_ ## type ## SV_with_SV_out_val(	   \
				sv1, offs, sv2, bg_val,			   \
				&k1, &k, &k2, &out_off, &out_val)))	   \
	{								   \
		neffrep += process_subassign_ ## type ## SV_out_val(ret,   \
				out_val, out_off, out_sv,		   \
				bg_val, sv1, k1);			   \
	}								   \
	return neffrep;							   \
}

DEFINE_subassign_typeSV_with_SV_FUN(int, int)
DEFINE_subassign_typeSV_with_SV_FUN(double, double)
DEFINE_subassign_typeSV_with_SV_FUN(Rcomplex, Rcomplex)
DEFINE_subassign_typeSV_with_SV_FUN(Rbyte, Rbyte)
DEFINE_subassign_typeSV_with_SV_FUN(character, SEXP)
DEFINE_subassign_typeSV_with_SV_FUN(list, SEXP)

/* 'sv1->len' and 'out_sv->len' must be the same.
   'sv1' and/or 'sv2' can be lacunar.
   'offs' must be an array of 'sv2->len' offsets (non-negative integers)
   that are sorted in strictly ascending order. The last offset in the
   array must be < 'sv1->len'.
   Returns the number of **effective** replacements, that is, the number of
   offsets for which the subassignment operation effectively modifies the
   original value. */
int _subassign_SV_with_SV(
		const SparseVec *sv1, const int *offs,
		const SparseVec *sv2, SparseVec *out_sv)
{
	SEXPTYPE Rtype = get_SV_subassign_Rtype(get_SV_Rtype(sv2), sv1, out_sv);
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		return subassign_intSV_with_SV(sv1, offs, sv2, out_sv);
	    case REALSXP:
		return subassign_doubleSV_with_SV(sv1, offs, sv2, out_sv);
	    case CPLXSXP:
		return subassign_RcomplexSV_with_SV(sv1, offs, sv2, out_sv);
	    case RAWSXP:
		return subassign_RbyteSV_with_SV(sv1, offs, sv2, out_sv);
	    case STRSXP:
		return subassign_characterSV_with_SV(sv1, offs, sv2, out_sv);
	    case VECSXP:
		return subassign_listSV_with_SV(sv1, offs, sv2, out_sv);
	}
	error("SparseArray internal error in _subassign_SV_with_SV():\n"
	      "    'out_sv' of type \"%s\" not supported", type2char(Rtype));
	return 0;  /* will never reach this */
}


/****************************************************************************
 * subassign_full_<type>SV_with_Rvector_block()
 */

static inline int next_k1(const SparseVec *sv1, int *k1, int out_off)
{
	if (*k1 < get_SV_nzcount(sv1) && sv1->nzoffs[*k1] == out_off) {
		(*k1)++;
		return 3;
	}
	return 2;
}

#define DEFINE_subassign_full_typeSV_with_Rvector_block_FUN(type)	  \
static int subassign_full_ ## type ## SV_with_Rvector_block(		  \
		const SparseVec *sv1,					  \
		const type *vals2, int cycle_len, SparseVec *out_sv)	  \
{									  \
	type bg_val = out_sv->na_background ? type ## NA : type ## 0;	  \
	out_sv->nzcount = 0;						  \
	int neffrep = 0, k1 = 0;					  \
	for (int out_off = 0; out_off < out_sv->len; out_off++) {	  \
		int ret = next_k1(sv1, &k1, out_off);			  \
		int i = cycle_len ? out_off % cycle_len : out_off;	  \
		type out_val = vals2[i];				  \
		neffrep += process_subassign_ ## type ## SV_out_val(ret,  \
				out_val, out_off, out_sv,		  \
				bg_val, sv1, k1);			  \
	}								  \
	return neffrep;							  \
}

DEFINE_subassign_full_typeSV_with_Rvector_block_FUN(int)
DEFINE_subassign_full_typeSV_with_Rvector_block_FUN(double)
DEFINE_subassign_full_typeSV_with_Rvector_block_FUN(Rcomplex)
DEFINE_subassign_full_typeSV_with_Rvector_block_FUN(Rbyte)

static int subassign_full_characterSV_with_Rvector_block(const SparseVec *sv1,
		SEXP Rvector, R_xlen_t block_offset, int cycle_len,
		SparseVec *out_sv)
{
	SEXP bg_val = out_sv->na_background ? characterNA : character0;
	out_sv->nzcount = 0;
	int neffrep = 0, k1 = 0;
	for (int out_off = 0; out_off < out_sv->len; out_off++) {
		int ret = next_k1(sv1, &k1, out_off);
		int i = cycle_len ? out_off % cycle_len :
				    block_offset + out_off;
		SEXP out_val = STRING_ELT(Rvector, i);
		neffrep += process_subassign_characterSV_out_val(ret,
				out_val, out_off, out_sv,
				bg_val, sv1, k1);
	}
	return neffrep;
}

static int subassign_full_listSV_with_Rvector_block(const SparseVec *sv1,
		SEXP Rvector, R_xlen_t block_offset, int cycle_len,
		SparseVec *out_sv)
{
	out_sv->nzcount = 0;
	int neffrep = 0, k1 = 0;
	for (int out_off = 0; out_off < out_sv->len; out_off++) {
		int ret = next_k1(sv1, &k1, out_off);
		int i = cycle_len ? out_off % cycle_len :
				    block_offset + out_off;
		SEXP out_val = VECTOR_ELT(Rvector, i);
		neffrep += process_subassign_listSV_out_val(ret,
				out_val, out_off, out_sv,
				list0, sv1, k1);
	}
	return neffrep;
}


/****************************************************************************
 * _subassign_full_SV_with_Rvector_block()
 */

/* Note that the content of input SparseVec 'sv' is used **only** to compute
   the number of **effective** replacements. In particular, it has NO impact
   whatsoever on the content that gets written to 'out_sv'.
   In other words, _subassign_full_SV_with_Rvector_block() is equivalent to:

     _write_Rvector_block_to_SV(Rvector, block_offset,
                                NULL, out_sv->len, out_sv)

   except that it uses the content of 'sv' to compute the number
   of **effective** replacements (which it returns).

   'sv->len' and 'out_sv->len' must be the same. 'sv' can be lacunar.
   The block of 'sv->len' (or 'out_sv->len') elements in 'Rvector' that
   starts at offset 'block_offset' forms the replacement value (a.k.a.
   right value) of the subassignment operation. It can contain zeros.
   About the length of 'Rvector':
   - If 'block_offset != 0' then its length must be >= 'block_offset + sv->len'.
   - If 'block_offset == 0' then 'Rvector' can be of any length (except 0
     unless 'sv->len' is also 0) and will be recycled if its length
     is < 'sv->len'.
   Returns the number of **effective** replacements, that is, the number of
   offsets for which the subassignment operation effectively modifies the
   original value. */
int _subassign_full_SV_with_Rvector_block(const SparseVec *sv,
		SEXP Rvector, R_xlen_t block_offset, SparseVec *out_sv)
{
	if (Rvector == R_NilValue)
		error("SparseArray internal error in "
		      "_subassign_full_SV_with_Rvector_block():\n"
		      "    'Rvector' cannot be NULL");
	SEXPTYPE Rtype = get_SV_subassign_Rtype(TYPEOF(Rvector), sv, out_sv);
	int Rvector_len = LENGTH(Rvector), cycle_len = 0;
	if (block_offset != 0) {
		if (Rvector_len < block_offset + sv->len)
			error("SparseArray internal error in "
			      "_subassign_full_SV_with_Rvector_block():\n"
			      "    LENGTH(Rvector) < block_offset + sv->len");
	} else if (Rvector_len < sv->len) {
		if (Rvector_len == 0)
			error("SparseArray internal error in "
			      "_subassign_full_SV_with_Rvector_block():\n"
			      "    LENGTH(Rvector) == 0");
		cycle_len = Rvector_len;
	}
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		return subassign_full_intSV_with_Rvector_block(sv,
				INTEGER(Rvector) + block_offset, cycle_len,
				out_sv);
	    case REALSXP:
		return subassign_full_doubleSV_with_Rvector_block(sv,
				REAL(Rvector) + block_offset, cycle_len,
				out_sv);
	    case CPLXSXP:
		return subassign_full_RcomplexSV_with_Rvector_block(sv,
				COMPLEX(Rvector) + block_offset, cycle_len,
				out_sv);
	    case RAWSXP:
		return subassign_full_RbyteSV_with_Rvector_block(sv,
				RAW(Rvector) + block_offset, cycle_len,
				out_sv);
	    case STRSXP:
		return subassign_full_characterSV_with_Rvector_block(sv,
				Rvector, block_offset, cycle_len,
				out_sv);
	    case VECSXP:
		return subassign_full_listSV_with_Rvector_block(sv,
				Rvector, block_offset, cycle_len,
				out_sv);
	}
	error("SparseArray internal error in "
	      "_subassign_full_SV_with_Rvector_block():\n"
	      "    'out_sv' of type \"%s\" not supported", type2char(Rtype));
	return 0;  /* will never reach this */
}

