/****************************************************************************
 *                 Basic manipulation of SparseVec structs                  *
 ****************************************************************************/
#include "SparseVec.h"

#include "Rvector_utils.h"


/* Does NOT work at the moment if 'Rtype' is VECSXP.
   IMPORTANT: The caller must immediately call 'PROTECT(sv.nzvals)' on
   the returned SparseVec struct when 'Rtype' is STRSXP. */
SparseVec _alloc_buf_SparseVec(SEXPTYPE Rtype, int len, int na_background)
{
	SparseVec sv;
	sv.Rtype = Rtype;
	if (Rtype == VECSXP)
		error("SparseArray internal error in "
		      "_alloc_buf_SparseVec():\n    type \"%s\" is "
		      "not supported at the moment", type2char(Rtype));
	if (Rtype == STRSXP) {
		sv.nzvals = PROTECT(NEW_CHARACTER(len));
	} else {
		size_t Rtype_size = _get_Rtype_size(Rtype);
		if (Rtype_size == 0)
			error("SparseArray internal error in "
			      "_alloc_buf_SparseVec():\n    type \"%s\" is "
			      "not supported", type2char(Rtype));
		if (na_background && Rtype == RAWSXP)
			error("SparseArray internal error in "
			      "_alloc_buf_SparseVec():\n    NaArray "
			      "objects of type \"raw\" are not supported");
		sv.nzvals = R_alloc(len, Rtype_size);
	}
	sv.nzoffs = (int *) R_alloc(len, sizeof(int));
	sv.nzcount = 0;
	sv.len = len;
	sv.na_background = na_background;
	if (Rtype == STRSXP)
		UNPROTECT(1);
	return sv;
}

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
		_set_selected_elts_to_one(INTSXP, out, 0,
				sv->nzoffs, get_SV_nzcount(sv));
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
		_set_selected_elts_to_one(REALSXP, out, 0,
				sv->nzoffs, get_SV_nzcount(sv));
	} else {  /* regular SparseVec */
		_copy_double_elts_to_offsets(nzvals_p,
				sv->nzoffs, get_SV_nzcount(sv), out);
	}
	return;
}

