/****************************************************************************
 *                   'Logic' operations on "leaf vectors"                   *
 ****************************************************************************/
#include "leaf_vector_Logic.h"

#include "leaf_vector_utils.h"

int _get_Logic_opcode(SEXP op)
{
	const char *s;

	if (!IS_CHARACTER(op) || LENGTH(op) != 1)
		error("'op' must be a single string");
	op = STRING_ELT(op, 0);
	if (op == NA_STRING)
		error("'op' cannot be NA");
	s = CHAR(op);
	if (strcmp(s, "&") == 0)
		return AND_OPCODE;
	if (strcmp(s, "|") == 0)
		return OR_OPCODE;
	error("SparseArray internal error in "
	      "_get_Logic_opcode():\n"
	      "    invalid op: \"%s\"", s);
	return 0;  /* will never reach this */
}

