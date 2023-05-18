#ifndef _RVECTOR_SUMMARIZATION_H_
#define _RVECTOR_SUMMARIZATION_H_

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

typedef struct summarize_op_t {
	int opcode;
	SEXPTYPE in_Rtype;   // only INTSXP/REALSXP supported for now
	int na_rm;
	double shift;
} SummarizeOp;

typedef union summarize_outbuf_t {
	int one_int[1];
	double one_double[1];
	int two_ints[2];
	double two_doubles[2];
	Rcomplex one_complex[1];  // not used yet
} SummarizeOutbuf;

typedef struct summarize_result_t {
	/* 'totalcount' is the length of the virtual vector we're summarizing.
	   We must have 0 <= nacount <= nzcount <= totalcount at any time. */
	R_xlen_t totalcount;
	R_xlen_t nzcount;
	/* 'nacount' is used only when 'summarize_op->na_rm' is True. */
	R_xlen_t nacount;
	/* 'outbuf_is_set' is used only when 'summarize_op->opcode' is
	   MIN_OPCODE, MAX_OPCODE, or RANGE_OPCODE, and 'summarize_op->in_Rtype'
	   is INTSXP. */
	int outbuf_is_set;
	SEXPTYPE out_Rtype;  // only LGLSXP/INTSXP/REALSXP supported for now
	SummarizeOutbuf outbuf;
	int warn;
} SummarizeResult;

int _get_summarize_opcode(SEXP op, SEXPTYPE Rtype);

SummarizeOp _make_SummarizeOp(
	int opcode,
	SEXPTYPE in_Rtype,
	int na_rm,
	double shift
);

void _init_SummarizeResult(
	const SummarizeOp *summarize_op,
	SummarizeResult *res
);

int _summarize_Rvector(
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

int _count_Rvector_NAs(SEXP Rvector);

int _Rvector_has_any_NA(SEXP Rvector);

#endif  /* _RVECTOR_SUMMARIZATION_H_ */

