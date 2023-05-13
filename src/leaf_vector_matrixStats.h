#ifndef _LEAF_VECTOR_MATRIXSTATS_H_
#define _LEAF_VECTOR_MATRIXSTATS_H_

#include <Rdefines.h>

/* matrixStats operations, grouped by interface */

/* Interface 1: FUN(x, dims) */
#define	COUNTNAS_OPCODE	 1	/* "CountNAs" (not a matrixStats op) */
#define	ANYNAS_OPCODE	 2	/* "AnyNAs" */

/* Interface 2: FUN(x, na.rm, dims) */
#define	SUMS_OPCODE	 3	/* "Sums" */
#define	SUMS2_OPCODE	 4	/* "Sums2" */
#define	MEANS_OPCODE	 5	/* "Means" */
#define	MEANS2_OPCODE	 6	/* "Means2" */
#define	MINS_OPCODE	 7	/* "Mins" */
#define	MAXS_OPCODE	 8	/* "Maxs" */
#define	RANGES_OPCODE	 9	/* "Ranges" */
#define	MEDIANS_OPCODE	10	/* "Medians" */
#define	ALLS_OPCODE	11	/* "Alls" */
#define	ANYS_OPCODE	12	/* "Anys" */

/* Interface 3: FUN(x, na.rm, center, dims) */
#define	VARS_OPCODE	13	/* "Vars" */
#define	SDS_OPCODE	14	/* "Sds" */

int _get_matrixStats_opcode(SEXP op);

#endif  /* _LEAF_VECTOR_MATRIXSTATS_H_ */

