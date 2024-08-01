#ifndef _RVECTOR_SUMMARIZATION_H_
#define _RVECTOR_SUMMARIZATION_H_

#include <Rdefines.h>

/* 3 interfaces to the summarization functions in R:
     o Interface 1: FUN(x)
     o Interface 2: FUN(x, na.rm)
     o Interface 3: FUN(x, na.rm, center) */

/* Interface 1: FUN(x) */
#define	ANYNA_OPCODE             1
#define	COUNTNAS_OPCODE          2

/* Operations from 'Summary' group generic.
   Interface 2: FUN(x, na.rm) */
#define	ANY_OPCODE               3
#define	ALL_OPCODE               4
#define	MIN_OPCODE               5
#define	MAX_OPCODE               6
#define	RANGE_OPCODE             7
#define	SUM_OPCODE               8  /* supports MEAN_OPCODE */
#define	PROD_OPCODE              9

/* Other "summarize" operations */
#define	MEAN_OPCODE             10  /* Interface 2 */
#define	CENTERED_X2_SUM_OPCODE  11  /* Interface 3, supports VAR1_OPCODE */
#define	SUM_X_X2_OPCODE         12  /* Interface 2, supports VAR2_OPCODE  */
#define	VAR1_OPCODE             13  /* Interface 3, supports SD1_OPCODE  */
#define	VAR2_OPCODE             14  /* Interface 2, supports SD2_OPCODE  */
#define	SD1_OPCODE              15  /* Interface 3 */
#define	SD2_OPCODE              16  /* Interface 2 */

#define RCOMPLEX_IS_NA(z) (ISNAN((z)->r) || ISNAN((z)->i))

typedef struct summarize_op_t {
	int opcode;
	SEXPTYPE in_Rtype;
	int na_rm;
	double center;
} SummarizeOp;

typedef union summarize_outbuf_t {
	int one_int[1];
	double one_double[1];
	int two_ints[2];
	double two_doubles[2];
	Rcomplex one_Rcomplex[1];  // not used yet
} SummarizeOutbuf;

/* Possible values for the 'outbuf_status' member below. */
#define	OUTBUF_IS_NOT_SET                  1
#define	OUTBUF_IS_SET                      2
#define	OUTBUF_IS_SET_WITH_BREAKING_VALUE  3

typedef struct summarize_result_t {
  /* 'in_length' is the length of the virtual vector we're summarizing.
     We must have 0 <= in_nacount <= in_nzcount <= in_length at any time. */
	R_xlen_t in_length;
	R_xlen_t in_nzcount;
  /* 'in_nacount' is used only when 'summarize_op->na_rm' is True. */
	R_xlen_t in_nacount;
	SEXPTYPE out_Rtype;  // only LGLSXP/INTSXP/REALSXP are supported
	int outbuf_status;   // see OUTBUF_* macros above for possible values
	SummarizeOutbuf outbuf;
	int postprocess_one_zero;
	int warn;
} SummarizeResult;

int _get_summarize_opcode(SEXP op, SEXPTYPE Rtype);

SummarizeOp _make_SummarizeOp(
	int opcode,
	SEXPTYPE in_Rtype,
	int na_rm,
	double center
);

void _init_SummarizeResult(
	const SummarizeOp *summarize_op,
	SummarizeResult *res
);

void _summarize_ones(
	int x_len,
	const SummarizeOp *summarize_op,
	SummarizeResult *res
);

void _summarize_Rvector(
	SEXP x,
	const SummarizeOp *summarize_op,
	SummarizeResult *res
);

void _postprocess_SummarizeResult(
	const SummarizeOp *summarize_op,
	SummarizeResult *res
);

SEXP _make_SEXP_from_summarize_result(
	const SummarizeOp *summarize_op,
	const SummarizeResult *res
);

#endif  /* _RVECTOR_SUMMARIZATION_H_ */

