/****************************************************************************
 *                     Subassignment of a sparse vector                     *
 ****************************************************************************/
#include "SparseVec_subassignment.h"


static inline int int_equal(int x, int y)
{
	return x == y;
}
static inline int double_equal(double x, double y)
{
	return x == y;
}
static inline int Rcomplex_equal(Rcomplex x, Rcomplex y)
{
	return x.r == y.r && x.i == y.i;
}


/****************************************************************************
 * _fill_SV_with_Rsubvec()
 *
 * TODO: Move this to its own file. Maybe _make_leaf_from_Rsubvec() and
 * _make_naleaf_from_Rsubvec() should be based on _fill_SV_with_Rsubvec().
 */

static void fill_intSV(const int *vals,
		const int *offs, int n, SparseVec *out_sv)
{
	int *out_nzvals = (int *) out_sv->nzvals;
	out_sv->nzcount = 0;
	int out_bg_val = out_sv->na_background ? intNA : int0;
	for (int k = 0; k < n; k++) {
		int out_val = vals[k];
		if (int_equal(out_val, out_bg_val))
			continue;
		int off = offs == NULL ? k : offs[k];
		APPEND_TO_NZVALS_NZOFFS(out_val, off,
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);
	}
	return;
}

static void fill_doubleSV(const double *vals,
		const int *offs, int n, SparseVec *out_sv)
{
	double *out_nzvals = (double *) out_sv->nzvals;
	out_sv->nzcount = 0;
	double out_bg_val = out_sv->na_background ? doubleNA : double0;
	for (int k = 0; k < n; k++) {
		double out_val = vals[k];
		if (double_equal(out_val, out_bg_val))
			continue;
		int off = offs == NULL ? k : offs[k];
		APPEND_TO_NZVALS_NZOFFS(out_val, off,
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);
	}
	return;
}

static void fill_RcomplexSV(const Rcomplex *vals,
		const int *offs, int n, SparseVec *out_sv)
{
	Rcomplex *out_nzvals = (Rcomplex *) out_sv->nzvals;
	out_sv->nzcount = 0;
	Rcomplex out_bg_val = out_sv->na_background ? RcomplexNA : Rcomplex0;
	for (int k = 0; k < n; k++) {
		Rcomplex out_val = vals[k];
		if (Rcomplex_equal(out_val, out_bg_val))
			continue;
		int off = offs == NULL ? k : offs[k];
		APPEND_TO_NZVALS_NZOFFS(out_val, off,
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);
	}
	return;
}

static void fill_RbyteSV(const Rbyte *vals,
		const int *offs, int n, SparseVec *out_sv)
{
	Rbyte *out_nzvals = (Rbyte *) out_sv->nzvals;
	out_sv->nzcount = 0;
	for (int k = 0; k < n; k++) {
		Rbyte out_val = vals[k];
		if (out_val == Rbyte0)
			continue;
		int off = offs == NULL ? k : offs[k];
		APPEND_TO_NZVALS_NZOFFS(out_val, off,
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);
	}
	return;
}

#define	IS_BG_CHARSXP(x, na_background) \
	(na_background) ? ((x) == NA_STRING) \
			: ((x) != NA_STRING && LENGTH(x) == 0)

static void fill_characterSV(SEXP Rvector, R_xlen_t subvec_offset,
		const int *offs, int n, SparseVec *out_sv)
{
	SEXP out_nzvals = (SEXP) out_sv->nzvals;
	out_sv->nzcount = 0;
	for (int k = 0; k < n; k++) {
		SEXP out_val = STRING_ELT(Rvector, subvec_offset + k);
		if (IS_BG_CHARSXP(out_val, out_sv->na_background))
			continue;
		int off = offs == NULL ? k : offs[k];
		SET_STRING_ELT(out_nzvals, out_sv->nzcount, out_val);
		out_sv->nzoffs[out_sv->nzcount] = off;
		out_sv->nzcount++;
	}
	return;
}

/* Fills 'out_sv' with the nonzero elements of 'Rvector' that have an
   index 'i' that is >= 'subvec_offset' and < 'subvec_offset + n'.
   'offs' must be NULL or an array of 'n' offsets (non-negative integers)
   that are strictly sorted (in ascending order). The last offset in the
   array must be < 'out_sv->len'. */
void _fill_SV_with_Rsubvec(SEXP Rvector, R_xlen_t subvec_offset,
		const int *offs, int n, SparseVec *out_sv)
{
	SEXPTYPE Rtype = get_SV_Rtype(out_sv);
	if (TYPEOF(Rvector) != Rtype)
		error("SparseArray internal error in "
		      "_fill_SV_with_Rsubvec():\n"
		      "    'Rvector' and 'out_sv' don't have the same type");
	if (offs == NULL) {
		if (n != out_sv->len)
			error("SparseArray internal error in "
			      "_fill_SV_with_Rsubvec():\n"
			      "    'offs == NULL' and 'n != out_sv->len'");
	} else {
		if (n > out_sv->len)
			error("SparseArray internal error in "
			      "_fill_SV_with_Rsubvec():\n"
			      "    'offs != NULL' and 'n > out_sv->len'");
	}
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		fill_intSV(INTEGER(Rvector) + subvec_offset,
			offs, n, out_sv);
		return;
	    case REALSXP:
		fill_doubleSV(REAL(Rvector) + subvec_offset,
			offs, n, out_sv);
		return;
	    case CPLXSXP:
		fill_RcomplexSV(COMPLEX(Rvector) + subvec_offset,
			offs, n, out_sv);
		return;
	    case RAWSXP:
		fill_RbyteSV(RAW(Rvector) + subvec_offset,
			offs, n, out_sv);
		return;
	    case STRSXP:
		fill_characterSV(Rvector, subvec_offset,
			offs, n, out_sv);
		return;
	}
	error("SparseArray internal error in "
	      "_fill_SV_with_Rsubvec():\n"
	      "    'out_sv' of type \"%s\" not supported yet",
	      type2char(Rtype));
	return;
}


/****************************************************************************
 * subassign_intSV()
 * subassign_doubleSV()
 * subassign_RcomplexSV()
 * subassign_RbyteSV()
 * subassign_characterSV()
 */

#define DEFINE_next_type_out_val_FUN(type)				\
static inline int next_ ## type ## _out_val(const SparseVec *sv1,	\
		const int *offs2, const type *vals2, int n2,		\
		int *k1, int *k2, int *off, type *out_val)		\
{									\
	int ret = next_offset(sv1->nzoffs, get_SV_nzcount(sv1),		\
			      offs2, n2, *k1, *k2, off);		\
	if (ret == 1) {							\
		*out_val = get_ ## type ## SV_nzval(sv1, *k1);		\
		(*k1)++;						\
	} else {							\
		*out_val = vals2[*k2];					\
		(*k2)++;						\
		if (ret == 3)						\
			(*k1)++;					\
	}								\
	return ret;							\
}

DEFINE_next_type_out_val_FUN(int)
DEFINE_next_type_out_val_FUN(double)
DEFINE_next_type_out_val_FUN(Rcomplex)
DEFINE_next_type_out_val_FUN(Rbyte)

static inline int next_character_out_val(
		const SparseVec *sv1,
		const int *offs2, SEXP Rvector, R_xlen_t subvec_offset, int n2,
		int *k1, int *k2, int *off, SEXP *out_val)
{
	int ret = next_offset(sv1->nzoffs, get_SV_nzcount(sv1),
			      offs2, n2, *k1, *k2, off);
	if (ret == 1) {
		*out_val = get_characterSV_nzval(sv1, *k1);
		(*k1)++;
	} else {
		*out_val = STRING_ELT(Rvector, subvec_offset + *k2);
		(*k2)++;
		if (ret == 3)
			(*k1)++;
	}
	return ret;
}

#define DEFINE_subassign_typeSV_FUN(type)				    \
static int subassign_ ## type ## SV(const SparseVec *sv1,		    \
		const int *offs2, const type *vals2, int n2,		    \
		SparseVec *out_sv)					    \
{									    \
	type *out_nzvals = (type *) out_sv->nzvals;			    \
	type out_bg_val = out_sv->na_background ? type ## NA : type ## 0;   \
	out_sv->nzcount = 0;						    \
	int ret, k1 = 0, k2 = 0, off, neffrep = 0;			    \
	type out_val;							    \
	while ((ret = next_ ## type ## _out_val(sv1, offs2, vals2, n2,	    \
						&k1, &k2, &off, &out_val))) \
	{								    \
		if (ret == 2) {						    \
			if (type ## _equal(out_val, out_bg_val)) 	    \
				continue;				    \
			neffrep++;					    \
		} else if (ret == 3) {					    \
			if (type ## _equal(out_val, out_bg_val)) {	    \
				neffrep++;				    \
				continue;				    \
			}						    \
			type v1 = get_ ## type ## SV_nzval(sv1, k1 - 1);    \
			if (!type ## _equal(v1, out_val))		    \
				neffrep++;				    \
		}							    \
		APPEND_TO_NZVALS_NZOFFS(out_val, off,			    \
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);	    \
	}								    \
	return neffrep;							    \
}

DEFINE_subassign_typeSV_FUN(int)
DEFINE_subassign_typeSV_FUN(double)
DEFINE_subassign_typeSV_FUN(Rcomplex)

static int subassign_RbyteSV(const SparseVec *sv1,
		const int *offs2, const Rbyte *vals2, int n2,
		SparseVec *out_sv)
{
	Rbyte *out_nzvals = (Rbyte *) out_sv->nzvals;
	out_sv->nzcount = 0;
	int ret, k1 = 0, k2 = 0, off, neffrep = 0;
	Rbyte out_val;
	while ((ret = next_Rbyte_out_val(sv1, offs2, vals2, n2,
					 &k1, &k2, &off, &out_val)))
	{
		if (ret == 2) {
			if (out_val == Rbyte0)
				continue;  /* zero replaces zero */
			/* nonzero replaces zero */
			neffrep++;
		} else if (ret == 3) {
			if (out_val == Rbyte0) {
				/* zero replaces nonzero */
				neffrep++;
				continue;
			}
			/* nonzero replaces nonzero */
			Rbyte v1 = get_RbyteSV_nzval(sv1, k1 - 1);
			if (v1 != out_val)
				neffrep++;
		}
		APPEND_TO_NZVALS_NZOFFS(out_val, off,
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);
	}
	return neffrep;
}

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
   still copy a leaf touched by the subassignment operation, even when the
   leaf has not changed. */
static int subassign_characterSV(const SparseVec *sv1,
		const int *offs2, SEXP Rvector, R_xlen_t subvec_offset, int n2,
		SparseVec *out_sv)
{
	SEXP out_nzvals = (SEXP) out_sv->nzvals;
	out_sv->nzcount = 0;
	int ret, k1 = 0, k2 = 0, off, neffrep = 0;
	SEXP out_val;
	while ((ret = next_character_out_val(sv1, offs2,
					     Rvector, subvec_offset, n2,
					     &k1, &k2, &off, &out_val)))
	{
		if (ret == 2) {
			if (IS_BG_CHARSXP(out_val, out_sv->na_background))
				continue;
			neffrep++;
		} else if (ret == 3) {
			if (IS_BG_CHARSXP(out_val, out_sv->na_background)) {
				neffrep++;
				continue;
			}
			SEXP v1 = get_characterSV_nzval(sv1, k1 - 1);
			/* See note above about this comparison. */
			if (v1 != out_val)
				neffrep++;
		}
		SET_STRING_ELT(out_nzvals, out_sv->nzcount, out_val);
		out_sv->nzoffs[out_sv->nzcount] = off;
		out_sv->nzcount++;
	}
	return neffrep;
}


/****************************************************************************
 * _subassign_SV_with_Rsubvec()
 */

/* 'sv->len' and 'out_sv->len' must be the same. 'sv' can be lacunar.
   'offs' must be an array of 'n' offsets (non-negative integers) that are
   strictly sorted (in ascending order). The last offset in the array must
   be < 'sv->len'.
   Elements of 'Rvector' with an index 'i' that is >= 'subvec_offset'
   and < 'subvec_offset + n' form the replacement value (a.k.a. right value)
   of the subassignment operation. It can contain zeros.
   Returns the number of **effective** replacements, that is, the number of
   offsets for which the subassignment operation effectively modifies the
   original value. */
int _subassign_SV_with_Rsubvec(const SparseVec *sv, const int *offs, int n,
		SEXP Rvector, R_xlen_t subvec_offset, SparseVec *out_sv)
{
	SEXPTYPE Rtype = get_SV_Rtype(sv);
	if (out_sv->len != sv->len || get_SV_Rtype(out_sv) != Rtype)
		error("SparseArray internal error in "
		      "_subassign_SV_with_Rsubvec():\n"
		      "    'sv' and 'out_sv' are incompatible");
	if (TYPEOF(Rvector) != Rtype)
		error("SparseArray internal error in "
		      "_subassign_SV_with_Rsubvec():\n"
		      "    'sv' and 'Rvector' don't have the same type");
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		return subassign_intSV(sv, offs,
				INTEGER(Rvector) + subvec_offset, n,
				out_sv);
	    case REALSXP:
		return subassign_doubleSV(sv, offs,
				REAL(Rvector) + subvec_offset, n,
				out_sv);
	    case CPLXSXP:
		return subassign_RcomplexSV(sv, offs,
				COMPLEX(Rvector) + subvec_offset, n,
				out_sv);
	    case RAWSXP:
		return subassign_RbyteSV(sv, offs,
				RAW(Rvector) + subvec_offset, n,
				out_sv);
	    case STRSXP:
		return subassign_characterSV(sv, offs,
				Rvector, subvec_offset, n,
				out_sv);
	}
	error("SparseArray internal error in "
	      "_subassign_SV_with_Rsubvec():\n"
	      "    type \"%s\" is not supported at the moment",
	      type2char(Rtype));
	return 0;  /* will never reach this */
}


/****************************************************************************
 * subassign_full_intSV()
 * subassign_full_doubleSV()
 * subassign_full_RcomplexSV()
 * subassign_full_RbyteSV()
 * subassign_full_characterSV()
 */

#define DEFINE_subassign_full_typeSV_FUN(type)				  \
static int subassign_full_ ## type ## SV(const SparseVec *sv1,		  \
		const type *vals2, SparseVec *out_sv)			  \
{									  \
	type *out_nzvals = (type *) out_sv->nzvals;			  \
	type out_bg_val = out_sv->na_background ? type ## NA : type ## 0; \
	out_sv->nzcount = 0;						  \
	int k1 = 0, neffrep = 0;					  \
	for (int i = 0; i < out_sv->len; i++) {				  \
		type v2 = vals2[i];					  \
		if (k1 < get_SV_nzcount(sv1) && sv1->nzoffs[k1] == i) {	  \
			type v1 = get_ ## type ## SV_nzval(sv1, k1);	  \
			k1++;						  \
			if (type ## _equal(v2, out_bg_val)) {		  \
				neffrep++;				  \
				continue;				  \
			}						  \
			if (!type ## _equal(v1, v2))			  \
				neffrep++;				  \
		} else {						  \
			if (type ## _equal(v2, out_bg_val))		  \
				continue;				  \
			neffrep++;					  \
		}							  \
		APPEND_TO_NZVALS_NZOFFS(v2, i,				  \
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);	  \
	}								  \
	return neffrep;							  \
}

DEFINE_subassign_full_typeSV_FUN(int)
DEFINE_subassign_full_typeSV_FUN(double)
DEFINE_subassign_full_typeSV_FUN(Rcomplex)

static int subassign_full_RbyteSV(const SparseVec *sv1,
		const Rbyte *vals2, SparseVec *out_sv)
{
	Rbyte *out_nzvals = (Rbyte *) out_sv->nzvals;
	out_sv->nzcount = 0;
	int k1 = 0, neffrep = 0;
	for (int i = 0; i < out_sv->len; i++) {
		Rbyte v2 = vals2[i];
		if (k1 < get_SV_nzcount(sv1) && sv1->nzoffs[k1] == i) {
			Rbyte v1 = get_RbyteSV_nzval(sv1, k1);
			k1++;
			if (v2 == Rbyte0) {
				/* zero replaces nonzero */
				neffrep++;
				continue;
			}
			/* nonzero replaces nonzero */
			if (v1 != v2)
				neffrep++;
		} else {
			if (v2 == Rbyte0)
				continue;  /* zero replaces zero */
			/* nonzero replaces zero */
			neffrep++;
		}
		APPEND_TO_NZVALS_NZOFFS(v2, i,
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);
	}
	return neffrep;
}

/* Note that when comparing CHARSXPs 'v1' and 'v2' below (v1 != v2), we
   compare their **addresses**, not their **values**.
   See note for subassign_characterSV() above for more information. */
static int subassign_full_characterSV(const SparseVec *sv1,
		SEXP Rvector, R_xlen_t subvec_offset, SparseVec *out_sv)
{
	SEXP out_nzvals = (SEXP) out_sv->nzvals;
	out_sv->nzcount = 0;
	int k1 = 0, neffrep = 0;
	for (int i = 0; i < out_sv->len; i++) {
		SEXP v2 = STRING_ELT(Rvector, subvec_offset + i);
		if (k1 < get_SV_nzcount(sv1) && sv1->nzoffs[k1] == i) {
			SEXP v1 = get_characterSV_nzval(sv1, k1);
			k1++;
			if (IS_BG_CHARSXP(v2, out_sv->na_background)) {
				neffrep++;
				continue;
			}
			/* See note above about this comparison. */
			if (v1 != v2)
				neffrep++;
		} else {
			if (IS_BG_CHARSXP(v2, out_sv->na_background))
				continue;
			neffrep++;
		}
		SET_STRING_ELT(out_nzvals, out_sv->nzcount, v2);
		out_sv->nzoffs[out_sv->nzcount] = i;
		out_sv->nzcount++;
	}
	return neffrep;
}


/****************************************************************************
 * _subassign_full_SV_with_Rsubvec()
 *
 */

/* Note that the content of input SparseVec 'sv' is used only to compute the
   number of **effective** replacements. In particular, it has NO impact on
   what content gets written to 'out_sv'.
   In other words, _subassign_full_SV_with_Rsubvec() is equivalent to:

     _fill_SV_with_Rsubvec(Rvector, subvec_offset, NULL, out_sv->len, out_sv)

   except that the former uses the content of 'sv' to compute the number
   of **effective** replacements and to return it.

   'sv->len' and 'out_sv->len' must be the same. 'sv' can be lacunar.
   Elements of 'Rvector' with an index 'i' that is >= 'subvec_offset'
   and < 'subvec_offset + out_sv->len' form the replacement value (a.k.a.
   right value) of the subassignment operation. It can contain zeros.
   Returns the number of **effective** replacements, that is, the number of
   offsets for which the subassignment operation effectively modifies the
   original value. */
int _subassign_full_SV_with_Rsubvec(const SparseVec *sv,
		SEXP Rvector, R_xlen_t subvec_offset, SparseVec *out_sv)
{
	SEXPTYPE Rtype = get_SV_Rtype(sv);
	if (out_sv->len != sv->len || get_SV_Rtype(out_sv) != Rtype)
		error("SparseArray internal error in "
		      "_subassign_full_SV_with_Rsubvec():\n"
		      "    'sv' and 'out_sv' are incompatible");
	if (TYPEOF(Rvector) != Rtype)
		error("SparseArray internal error in "
		      "_subassign_full_SV_with_Rsubvec():\n"
		      "    'sv' and 'Rvector' don't have the same type");
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		return subassign_full_intSV(sv,
				INTEGER(Rvector) + subvec_offset, out_sv);
	    case REALSXP:
		return subassign_full_doubleSV(sv,
				REAL(Rvector) + subvec_offset, out_sv);
	    case CPLXSXP:
		return subassign_full_RcomplexSV(sv,
				COMPLEX(Rvector) + subvec_offset, out_sv);
	    case RAWSXP:
		return subassign_full_RbyteSV(sv,
				RAW(Rvector) + subvec_offset, out_sv);
	    case STRSXP:
		return subassign_full_characterSV(sv,
				Rvector, subvec_offset, out_sv);
	}
	error("SparseArray internal error in "
	      "_subassign_full_SV_with_Rsubvec():\n"
	      "    type \"%s\" is not supported at the moment",
	      type2char(Rtype));
	return 0;  /* will never reach this */
}

