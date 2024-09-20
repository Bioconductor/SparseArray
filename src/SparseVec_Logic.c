/****************************************************************************
 *                   'Logic' operations on sparse vectors                   *
 ****************************************************************************/
#include "SparseVec_Logic.h"

#include "SparseVec.h"


int _get_Logic_opcode(SEXP op)
{
	const char *s;

	if (!IS_CHARACTER(op) || LENGTH(op) != 1)
		error("SparseArray internal error in _get_Logic_opcode():\n"
		      "    'op' must be a single string");
	op = STRING_ELT(op, 0);
	if (op == NA_STRING)
		error("SparseArray internal error in _get_Logic_opcode():\n"
		      "    'op' cannot be NA");
	s = CHAR(op);
	if (strcmp(s, "&") == 0)
		return AND_OPCODE;
	if (strcmp(s, "|") == 0)
		return OR_OPCODE;
	error("SparseArray internal error in _get_Logic_opcode():\n"
	      "    invalid op: \"%s\"", s);
	return 0;  /* will never reach this */
}

static inline int Logic_int_int(int opcode, int x, int y)
{
	switch (opcode) {
	    case AND_OPCODE:
		if (x == 0 || y == 0)
			return 0;
		if (x == NA_INTEGER || y == NA_INTEGER)
			return NA_INTEGER;
		return 1;
	    case OR_OPCODE:
		if (x == 1 || y == 1)
			return 1;
		if (x == NA_INTEGER || y == NA_INTEGER)
			return NA_INTEGER;
		return 0;
	}
	error("SparseArray internal error in Logic_int_int():\n"
	      "    unsupported 'opcode'");
	return 0;  /* will never reach this */
}

void _Logic_intSV_na(int opcode,
		const SparseVec *sv1, SEXPTYPE Rtype2, SparseVec *out_sv)
{
	if (out_sv->len != sv1->len)
		error("SparseArray internal error in "
		      "_Logic_intSV_na():\n"
		      "    'sv1' and 'out_sv' are incompatible");
	int *out_nzvals = (int *) out_sv->nzvals;
	out_sv->nzcount = 0;
	int out_background = out_sv->na_background ? intNA : int0;
	const int *nzvals1_p = get_intSV_nzvals_p(sv1);
	if (nzvals1_p == NULL) {  /* lacunar SparseVec */
		int out_val = Logic_int_int(opcode, int1, intNA);
		if (out_val == out_background)
			return;
		out_nzvals[0] = out_val;
		out_sv->nzcount = PROPAGATE_NZOFFS;
		return;
	}
	/* regular SparseVec */
	int nzcount1 = get_SV_nzcount(sv1);
	for (int k = 0; k < nzcount1; k++) {
		int out_val = Logic_int_int(opcode, nzvals1_p[k], intNA);
		if (out_val == out_background)
			continue;
		APPEND_TO_NZVALS_NZOFFS(out_val, sv1->nzoffs[k],
				out_nzvals, out_sv->nzoffs, out_sv->nzcount);
	}
	return;
}

void _Logic_intSV_intSV(int opcode,
		const SparseVec *sv1, const SparseVec *sv2, SparseVec *out_sv)
{
	if (out_sv->len != sv1->len || out_sv->len != sv2->len)
		error("SparseArray internal error in "
		      "_Logic_intSV_intSV():\n"
		      "    'sv1', 'sv2', and 'out_sv' are incompatible");
	int *out_nzvals = (int *) out_sv->nzvals;
	out_sv->nzcount = 0;
	int out_background = out_sv->na_background ? intNA : int0;
	int k1 = 0, k2 = 0, off, x, y;
	while (next_int_int_vals(sv1, sv2, &k1, &k2, &off, &x, &y)) {
		int out_val = Logic_int_int(opcode, x, y);
		if (out_val == out_background)
			continue;
		APPEND_TO_NZVALS_NZOFFS(out_val, off,
				out_nzvals, out_sv->nzoffs, out_sv->nzcount);
	}
	return;
}

