/****************************************************************************
 *                      Subsetting of a sparse vector                       *
 ****************************************************************************/
#include "SparseVec_subsetting.h"


/* Returns a value >= 0 and < 'nzcount' if success, or -1 if failure. */
static inline int bsearch_off_to_k1(int off, const int *nzoffs, int nzcount)
{
	/* Compare with first offset. */
	int k1 = 0;
	int nzoff = nzoffs[k1];
	if (off < nzoff)
		return -1;
	if (off == nzoff)
		return k1;

	/* Compare with last offset. */
	int k2 = nzcount - 1;
	nzoff = nzoffs[k2];
	if (off > nzoff)
		return -1;
	if (off == nzoff)
		return k2;

	/* Binary search.
	   Seems that using >> 1 instead of / 2 is faster, even when compiling
	   with 'gcc -O2' (one would hope that the optimizer is able to do that
	   kind of optimization). */
	int k;
	while ((k = (k1 + k2) >> 1) != k1) {
		nzoff = nzoffs[k];
		if (off == nzoff)
			return k;
		if (off > nzoff)
			k1 = k;
		else
			k2 = k;
	}
	return -1;
}

static void build_lookup_table(int *lookup_table,
		const int *nzoffs, int nzcount)
{
	for (int k1 = 0; k1 < nzcount; k1++)
		lookup_table[*(nzoffs++)] = k1;
	return;
}

static void reset_lookup_table(int *lookup_table,
		const int *nzoffs, int nzcount)
{
	for (int k = 0; k < nzcount; k++)
		lookup_table[*(nzoffs)++] = -1;
	return;
}

/* Mapping 'off' to k1 thru the lookup table is always faster than using a
   binary search. However, allocating and building the lookup table has a
   small cost (that is proportional to 'sv->len'), so is only worth it if
   we're going to use it to map more than one 'off' value.
   See C_subset_SVT_by_Noffs() in SparseArray_subsetting.c for the (naive)
   strategy that is used at the moment to decide whether to use a lookup
   table or not. */
#define MAP_OFF_TO_K1(off, nzoffs, nzcount) \
	lookup_table == NULL ? bsearch_off_to_k1((off), (nzoffs), (nzcount)) \
			     : lookup_table[(off)]

#define DEFINE_subset_typeSV_FUN(type)					   \
static void subset_ ## type ## SV(const SparseVec *sv, const int *offs,	   \
		SparseVec *out_sv, int *lookup_table)			   \
{									   \
	out_sv->nzcount = 0;						   \
	type *out_nzvals = (type *) out_sv->nzvals;			   \
	if (lookup_table != NULL)					   \
		build_lookup_table(lookup_table, sv->nzoffs, sv->nzcount); \
	for (int out_off = 0; out_off < out_sv->len; out_off++) {	   \
		int off = offs[out_off];				   \
		int k1 = MAP_OFF_TO_K1(off, sv->nzoffs, sv->nzcount);	   \
		if (k1 < 0)						   \
			continue;					   \
		type out_val = get_ ## type ## SV_nzval(sv, k1);	   \
		APPEND_TO_NZVALS_NZOFFS(out_val, out_off,		   \
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);	   \
	}								   \
        if (lookup_table != NULL)					   \
		reset_lookup_table(lookup_table, sv->nzoffs, sv->nzcount); \
	return;								   \
}

DEFINE_subset_typeSV_FUN(int)
DEFINE_subset_typeSV_FUN(double)
DEFINE_subset_typeSV_FUN(Rcomplex)
DEFINE_subset_typeSV_FUN(Rbyte)

static void subset_characterSV(const SparseVec *sv, const int *offs,
		SparseVec *out_sv, int *lookup_table)
{
	out_sv->nzcount = 0;
	SEXP out_nzvals = (SEXP) out_sv->nzvals;
	if (lookup_table != NULL)
		build_lookup_table(lookup_table, sv->nzoffs, sv->nzcount);
	for (int out_off = 0; out_off < out_sv->len; out_off++) {
		int off = offs[out_off];
		int k1 = MAP_OFF_TO_K1(off, sv->nzoffs, sv->nzcount);
		if (k1 < 0)
			continue;
		SEXP out_val = get_characterSV_nzval(sv, k1);
		SET_STRING_ELT(out_nzvals, out_sv->nzcount, out_val);
		out_sv->nzoffs[out_sv->nzcount] = out_off;
		out_sv->nzcount++;
	}
        if (lookup_table != NULL)
		reset_lookup_table(lookup_table, sv->nzoffs, sv->nzcount);
	return;
}

static void subset_listSV(const SparseVec *sv, const int *offs,
		SparseVec *out_sv, int *lookup_table)
{
	out_sv->nzcount = 0;
	SEXP out_nzvals = (SEXP) out_sv->nzvals;
	if (lookup_table != NULL)
		build_lookup_table(lookup_table, sv->nzoffs, sv->nzcount);
	for (int out_off = 0; out_off < out_sv->len; out_off++) {
		int off = offs[out_off];
		int k1 = MAP_OFF_TO_K1(off, sv->nzoffs, sv->nzcount);
		if (k1 < 0)
			continue;
		SEXP out_val = get_listSV_nzval(sv, k1);
		SET_VECTOR_ELT(out_nzvals, out_sv->nzcount, out_val);
		out_sv->nzoffs[out_sv->nzcount] = out_off;
		out_sv->nzcount++;
	}
        if (lookup_table != NULL)
		reset_lookup_table(lookup_table, sv->nzoffs, sv->nzcount);
	return;
}

/* 'offs' must be an array of 'out_sv->len' offsets that are >= 0
   and < 'sv->len'.
   'lookup_table' must be NULL or an array of 'sv->len' integers
   filled with -1.
   Note that background values 'sv->na_background' and 'out_sv->na_background'
   should be the same. However, we never make any use of them (background
   does not matter in the context of N-index subsetting), so we don't bother
   checking them.
 */
void _subset_SV(const SparseVec *sv, const int *offs,
		SparseVec *out_sv, int *lookup_table)
{
	SEXPTYPE Rtype = get_SV_Rtype(sv);
	if (get_SV_Rtype(out_sv) != Rtype)
		error("SparseArray internal error in _subset_SV():\n"
		      "    'sv' and 'out_sv' must have the same type");
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		return subset_intSV(sv, offs, out_sv, lookup_table);
	    case REALSXP:
		return subset_doubleSV(sv, offs, out_sv, lookup_table);
	    case CPLXSXP:
		return subset_RcomplexSV(sv, offs, out_sv, lookup_table);
	    case RAWSXP:
		return subset_RbyteSV(sv, offs, out_sv, lookup_table);
	    case STRSXP:
		return subset_characterSV(sv, offs, out_sv, lookup_table);
	    case VECSXP:
		return subset_listSV(sv, offs, out_sv, lookup_table);
	}
	error("SparseArray internal error in _subset_SV():\n"
	      "    'sv' of type \"%s\" not supported", type2char(Rtype));
        return;  /* will never reach this */
}

