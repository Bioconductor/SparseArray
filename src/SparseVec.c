/****************************************************************************
 *                 Basic manipulation of SparseVec structs                  *
 ****************************************************************************/
#include "SparseVec.h"

#include "Rvector_utils.h"


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

