/****************************************************************************
 *                 matrixStats operations on "leaf vectors"                 *
 ****************************************************************************/
#include "leaf_vector_matrixStats.h"

#include "leaf_vector_utils.h"

int _get_matrixStats_opcode(SEXP op)
{
	const char *s;

	if (!IS_CHARACTER(op) || LENGTH(op) != 1)
		error("SparseArray internal error in "
		      "_get_matrixStats_opcode():\n"
		      "    'op' must be a single string");
	op = STRING_ELT(op, 0);
	if (op == NA_STRING)
		error("SparseArray internal error in "
		      "_get_matrixStats_opcode():\n"
		      "    'op' cannot be NA");
	s = CHAR(op);
	if (strcmp(s, "CountNAs") == 0)
		return COUNTNAS_OPCODE;
	if (strcmp(s, "AnyNAs") == 0)
		return ANYNAS_OPCODE;
	if (strcmp(s, "Sums") == 0)
		return SUMS_OPCODE;
	if (strcmp(s, "Sums2") == 0)
		return SUMS2_OPCODE;
	if (strcmp(s, "Means") == 0)
		return MEANS_OPCODE;
	if (strcmp(s, "Means2") == 0)
		return MEANS2_OPCODE;
	if (strcmp(s, "Mins") == 0)
		return MINS_OPCODE;
	if (strcmp(s, "Maxs") == 0)
		return MAXS_OPCODE;
	if (strcmp(s, "Ranges") == 0)
		return RANGES_OPCODE;
	if (strcmp(s, "Medians") == 0)
		return MEDIANS_OPCODE;
	if (strcmp(s, "Alls") == 0)
		return ALLS_OPCODE;
	if (strcmp(s, "Anys") == 0)
		return ANYS_OPCODE;
	if (strcmp(s, "Vars") == 0)
		return VARS_OPCODE;
	if (strcmp(s, "Sds") == 0)
		return SDS_OPCODE;
	error("SparseArray internal error in _get_matrixStats_opcode():\n"
	      "    invalid op: \"%s\"", s);
	return 0;  /* will never reach this */
}

