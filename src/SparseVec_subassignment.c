/****************************************************************************
 *                     Subassignment of a sparse vector                     *
 ****************************************************************************/
#include "SparseVec_subassignment.h"


/****************************************************************************
 * _fill_SV_with_vals()
 *
 * TODO: Move this to its own file. Maybe _make_leaf_from_Rsubvec() and
 * _make_naleaf_from_Rsubvec() should be based on _fill_SV_with_vals().
 */

static void fill_SV_with_ints(const int *vals, const int *offs,
		int n, SparseVec *out_sv)
{
	int *out_nzvals = (int *) out_sv->nzvals;
	out_sv->nzcount = 0;
	int out_bg_val = out_sv->na_background ? intNA : int0;
	for (int k = 0; k < n; k++) {
		int out_val = vals[k];
		if (out_val == out_bg_val)
			continue;
		int off = offs == NULL ? k : offs[k];
		APPEND_TO_NZVALS_NZOFFS(out_val, off,
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);
	}
	return;
}

static void fill_SV_with_doubles(const double *vals, const int *offs,
		int n, SparseVec *out_sv)
{
	double *out_nzvals = (double *) out_sv->nzvals;
	out_sv->nzcount = 0;
	double out_bg_val = out_sv->na_background ? doubleNA : double0;
	for (int k = 0; k < n; k++) {
		double out_val = vals[k];
		if (out_val == out_bg_val)
			continue;
		int off = offs == NULL ? k : offs[k];
		APPEND_TO_NZVALS_NZOFFS(out_val, off,
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);
	}
	return;
}

static void fill_SV_with_Rcomplexes(const Rcomplex *vals, const int *offs,
		int n, SparseVec *out_sv)
{
	Rcomplex *out_nzvals = (Rcomplex *) out_sv->nzvals;
	out_sv->nzcount = 0;
	Rcomplex out_bg_val = out_sv->na_background ? RcomplexNA : Rcomplex0;
	for (int k = 0; k < n; k++) {
		Rcomplex out_val = vals[k];
		if (out_val.r == out_bg_val.r && out_val.i == out_bg_val.i)
			continue;
		int off = offs == NULL ? k : offs[k];
		APPEND_TO_NZVALS_NZOFFS(out_val, off,
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);
	}
	return;
}

static void fill_SV_with_Rbytes(const Rbyte *vals, const int *offs,
		int n, SparseVec *out_sv)
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

/* 'vals' must be an array of 'n' values of the same type as 'out_sv'.
   Zeros in 'vals' are ok.
   'offs' must be NULL or an array of 'n' integers. The values in 'offs'
   must be valid zero-based indices into 'out_sv' sorted in strictly
   ascending order. */
void _fill_SV_with_vals(const void *vals, const int *offs, int n,
		SparseVec *out_sv)
{
	if (offs == NULL) {
		if (n != out_sv->len)
			error("SparseArray internal error in "
			      "_fill_SV_with_vals():\n"
			      "    'offs == NULL' and 'n != out_sv->len'");
	} else {
		if (n > out_sv->len)
			error("SparseArray internal error in "
			      "_fill_SV_with_vals():\n"
			      "    'offs != NULL' and 'n > out_sv->len'");
	}
	SEXPTYPE Rtype = get_SV_Rtype(out_sv);
	switch (Rtype) {
	    case INTSXP: case LGLSXP:
		fill_SV_with_ints((const int *) vals, offs,
				n, out_sv);
		return;
	    case REALSXP:
		fill_SV_with_doubles((const double *) vals, offs,
				n, out_sv);
		return;
	    case CPLXSXP:
		fill_SV_with_Rcomplexes((const Rcomplex *) vals, offs,
				n, out_sv);
		return;
	    case RAWSXP:
		fill_SV_with_Rbytes((const Rbyte *) vals, offs,
				n, out_sv);
		return;
	}
	error("SparseArray internal error in "
	      "_fill_SV_with_vals():\n"
	      "    'out_sv' of type \"%s\" not supported yet",
	      type2char(Rtype));
	return;
}


/****************************************************************************
 * subassign_intSV_with_ints()
 * subassign_doubleSV_with_doubles()
 * subassign_RcomplexSV_with_Rcomplexes()
 * subassign_RbyteSV_with_Rbytes()
 */

#define DEFINE_next_type_out_val_FUN(type)				\
static inline int next_ ## type ## _out_val(				\
	const SparseVec *sv1,						\
	const int *offs2, const type *vals2, int n2,			\
	int *k1, int *k2, int *off, type *out_val)			\
{									\
	int ret = next_offset(sv1->nzoffs, get_SV_nzcount(sv1),		\
			      offs2, n2,				\
			      *k1, *k2, off);				\
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

static void subassign_intSV_with_ints(const SparseVec *sv1,
		const int *offs2, const int *vals2, int n2,
		SparseVec *out_sv)
{
	int *out_nzvals = (int *) out_sv->nzvals;
	out_sv->nzcount = 0;
	int out_bg_val = out_sv->na_background ? intNA : int0;
	int k1 = 0, k2 = 0, off;
	int out_val;
	while ((next_int_out_val(sv1, offs2, vals2, n2,
				 &k1, &k2, &off, &out_val)))
	{
		if (out_val == out_bg_val)
			continue;
		APPEND_TO_NZVALS_NZOFFS(out_val, off,
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);
	}
	return;
}

static void subassign_doubleSV_with_doubles(const SparseVec *sv1,
		const int *offs2, const double *vals2, int n2,
		SparseVec *out_sv)
{
	double *out_nzvals = (double *) out_sv->nzvals;
	out_sv->nzcount = 0;
	double out_bg_val = out_sv->na_background ? doubleNA : double0;
	int k1 = 0, k2 = 0, off;
	double out_val;
	while ((next_double_out_val(sv1, offs2, vals2, n2,
				    &k1, &k2, &off, &out_val)))
	{
		if (out_val == out_bg_val)
			continue;
		APPEND_TO_NZVALS_NZOFFS(out_val, off,
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);
	}
	return;
}

static void subassign_RcomplexSV_with_Rcomplexes(const SparseVec *sv1,
		const int *offs2, const Rcomplex *vals2, int n2,
		SparseVec *out_sv)
{
	Rcomplex *out_nzvals = (Rcomplex *) out_sv->nzvals;
	out_sv->nzcount = 0;
	Rcomplex out_bg_val = out_sv->na_background ? RcomplexNA : Rcomplex0;
	int k1 = 0, k2 = 0, off;
	Rcomplex out_val;
	while ((next_Rcomplex_out_val(sv1, offs2, vals2, n2,
				      &k1, &k2, &off, &out_val)))
	{
		if (out_val.r == out_bg_val.r && out_val.i == out_bg_val.i)
			continue;
		APPEND_TO_NZVALS_NZOFFS(out_val, off,
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);
	}
	return;
}

static void subassign_RbyteSV_with_Rbytes(const SparseVec *sv1,
		const int *offs2, const Rbyte *vals2, int n2,
		SparseVec *out_sv)
{
	Rbyte *out_nzvals = (Rbyte *) out_sv->nzvals;
	out_sv->nzcount = 0;
	int k1 = 0, k2 = 0, off;
	Rbyte out_val;
	while ((next_Rbyte_out_val(sv1, offs2, vals2, n2,
				   &k1, &k2, &off, &out_val)))
	{
		if (out_val == Rbyte0)
			continue;
		APPEND_TO_NZVALS_NZOFFS(out_val, off,
			out_nzvals, out_sv->nzoffs, out_sv->nzcount);
	}
	return;
}


/****************************************************************************
 * _subassign_SV1_with_v2()
 */

/* 'offs2' must be an array of 'n2' integers. The values in 'offs2' must be
   valid zero-based indices into 'sv1' sorted in strictly ascending order.
   'vals2' must be an array of 'n2' values of the same type as the values
   in 'sv1'. Zeros in 'vals2' are ok.
   Ok to use on a lacunar SparseVec. */
void _subassign_SV1_with_v2(const SparseVec *sv1,
		const int *offs2, const void *vals2, int n2,
		SparseVec *out_sv)
{
	SEXPTYPE Rtype1 = get_SV_Rtype(sv1);
	if (get_SV_Rtype(out_sv) != Rtype1 || out_sv->len != sv1->len)
		error("SparseArray internal error in "
		      "_subassign_SV1_with_v2():\n"
		      "    'sv1' and 'out_sv' are incompatible");
	switch (Rtype1) {
	    case INTSXP: case LGLSXP:
		subassign_intSV_with_ints(sv1,
				offs2, (const int *)      vals2, n2, out_sv);
		return;
	    case REALSXP:
		subassign_doubleSV_with_doubles(sv1,
				offs2, (const double *)   vals2, n2, out_sv);
		return;
	    case CPLXSXP:
		subassign_RcomplexSV_with_Rcomplexes(sv1,
				offs2, (const Rcomplex *) vals2, n2, out_sv);
		return;
	    case RAWSXP:
		subassign_RbyteSV_with_Rbytes(sv1,
				offs2, (const Rbyte *)    vals2, n2, out_sv);
		return;
	}
	error("SparseArray internal error in "
	      "_subassign_SV1_with_v2():\n"
	      "    'sv1' of type \"%s\" not supported yet",
	      type2char(Rtype1));
	return;
}

