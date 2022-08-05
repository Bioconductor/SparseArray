#ifndef _RVECTOR_SUMMARIZE_H_
#define _RVECTOR_SUMMARIZE_H_

#include <Rdefines.h>

/* "Summarize" operations from Summary group generic */
#define	MIN_OPCODE      1
#define	MAX_OPCODE      2
#define	RANGE_OPCODE    3
#define	SUM_OPCODE      4
#define	PROD_OPCODE     5
#define	ANY_OPCODE      6
#define	ALL_OPCODE      7

/* Other "summarize" operations */
#define	SUM_SHIFTED_X2_OPCODE 8  /* to support var1() */
#define	SUM_X_X2_OPCODE       9  /* to support var2() */

int _get_summarize_opcode(SEXP op, SEXPTYPE Rtype);

typedef int (*SummarizeInts_FUNType)(
	void *init, const int *x, int n, int na_rm, R_xlen_t *na_rm_count,
	int status);

typedef int (*SummarizeDoubles_FUNType)(
	void *init, const double *x, int n, int na_rm, R_xlen_t *na_rm_count,
	int status);

typedef struct summarize_op_t {
	int opcode;
	SEXPTYPE Rtype;  /* only INTSXP or REALSXP at the moment */
	int na_rm;
	double shift;
	SummarizeInts_FUNType summarize_ints_FUN;
	SummarizeDoubles_FUNType summarize_doubles_FUN;
} SummarizeOp;

SummarizeOp _init_SummarizeOp(
	int opcode,
	SEXPTYPE Rtype,
	int na_rm,
	double shift,
	void *init
);

int _apply_summarize_op(
	const SummarizeOp *summarize_op,
	void *init,
	const void *x,
	int n,
	R_xlen_t *na_rm_count,
	int status
);

SEXP _make_SEXP_from_summarize_result(
	int opcode,
	SEXPTYPE Rtype,
	void *init,
	int na_rm,
	R_xlen_t na_rm_count,
	int status
);

int _count_Rvector_NAs(SEXP Rvector);

int _Rvector_has_any_NA(SEXP Rvector);

#endif  /* _RVECTOR_SUMMARIZE_H_ */

