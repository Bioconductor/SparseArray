#include <R_ext/Rdynload.h>

#include "sparseMatrix_utils.h"
#include "SparseArray_class.h"
#include "SVT_SparseArray_class.h"
#include "SparseArray_dim_tuning.h"
#include "SparseArray_subsetting.h"
#include "SparseArray_subassignment.h"
#include "SparseArray_combine.h"
#include "SparseArray_summarization.h"
#include "SparseArray_Ops_methods.h"
#include "SparseArray_Math_methods.h"
#include "SparseArray_Complex_methods.h"
#include "SparseMatrix_mult.h"
#include "matrixStats_methods.h"
#include "rowsum_methods.h"
#include "randomSparseArray.h"
#include "readSparseCSV.h"
#include "test.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* sparseMatrix_utils.c */
	CALLMETHOD_DEF(C_colMins_dgCMatrix, 2),
	CALLMETHOD_DEF(C_colMaxs_dgCMatrix, 2),
	CALLMETHOD_DEF(C_colRanges_dgCMatrix, 2),
	CALLMETHOD_DEF(C_colVars_dgCMatrix, 2),

/* SparseArray_class.c */
	CALLMETHOD_DEF(C_coercion_can_introduce_zeros, 2),

/* SVT_SparseArray_class.c */
	CALLMETHOD_DEF(C_get_SVT_SparseArray_nzcount, 2),
	CALLMETHOD_DEF(C_set_SVT_SparseArray_type, 4),
	CALLMETHOD_DEF(C_from_SVT_SparseArray_to_Rarray, 4),
	CALLMETHOD_DEF(C_build_SVT_from_Rarray, 2),
	CALLMETHOD_DEF(C_from_SVT_SparseMatrix_to_CsparseMatrix, 3),
	CALLMETHOD_DEF(C_build_SVT_from_CsparseMatrix, 2),
	CALLMETHOD_DEF(C_from_SVT_SparseArray_to_COO_SparseArray, 3),
	CALLMETHOD_DEF(C_transpose_SVT_SparseMatrix, 3),

/* SparseArray_dim_tuning.c */
	CALLMETHOD_DEF(C_tune_SVT_dims, 4),

/* SparseArray_subsetting.c */
	CALLMETHOD_DEF(C_subset_SVT_SparseArray, 4),

/* SparseArray_subassignment.c */
	CALLMETHOD_DEF(C_subassign_SVT_by_Mindex, 5),
	CALLMETHOD_DEF(C_subassign_SVT_by_Lindex, 5),
	CALLMETHOD_DEF(C_subassign_SVT_with_short_Rvector, 5),
	CALLMETHOD_DEF(C_subassign_SVT_with_Rarray, 5),
	CALLMETHOD_DEF(C_subassign_SVT_with_SVT, 7),

/* SparseArray_combine.c */
	CALLMETHOD_DEF(C_abind_SVT_SparseArray_objects, 3),

/* SparseArray_summarization.c */
	CALLMETHOD_DEF(C_summarize_SVT_SparseArray, 6),
	CALLMETHOD_DEF(C_count_SVT_SparseArray_NAs, 3),
	CALLMETHOD_DEF(C_anyNA_SVT_SparseArray, 3),

/* SparseArray_Ops_methods.c */
	CALLMETHOD_DEF(C_unary_minus_SVT, 3),
	CALLMETHOD_DEF(C_Arith_SVT1_v2, 6),
	CALLMETHOD_DEF(C_Compare_SVT1_v2, 5),
	CALLMETHOD_DEF(C_Arith_SVT1_SVT2, 8),
	CALLMETHOD_DEF(C_Compare_SVT1_SVT2, 7),
	CALLMETHOD_DEF(C_Logic_SVT1_SVT2, 7),

/* SparseArray_Math_methods.c */
	CALLMETHOD_DEF(C_Math_SVT, 4),
	CALLMETHOD_DEF(C_Math2_SVT, 5),

/* SparseArray_Complex_methods.c */
	CALLMETHOD_DEF(C_Complex_SVT, 4),

/* SparseMatrix_mult.c */
	CALLMETHOD_DEF(C_crossprod2_SVT_mat, 7),
	CALLMETHOD_DEF(C_crossprod2_mat_SVT, 7),
	CALLMETHOD_DEF(C_crossprod2_SVT_SVT, 8),
	CALLMETHOD_DEF(C_crossprod1_SVT, 5),

/* matrixStats_methods.h.c */
	CALLMETHOD_DEF(C_colStats1_SVT, 6),
	CALLMETHOD_DEF(C_colStats2_SVT, 7),
	CALLMETHOD_DEF(C_colStats3_SVT, 8),

/* rowsum_methods.c */
	CALLMETHOD_DEF(C_rowsum_SVT, 6),
	CALLMETHOD_DEF(C_rowsum_dgCMatrix, 4),

/* randomSparseArray.c */
	CALLMETHOD_DEF(C_simple_rpois, 2),
	CALLMETHOD_DEF(C_poissonSparseArray, 2),

/* readSparseCSV.c */
	CALLMETHOD_DEF(C_readSparseCSV_as_SVT_SparseMatrix, 5),

/* test.c */
	CALLMETHOD_DEF(C_test, 0),

	{NULL, NULL, 0}
};

void R_init_SparseArray(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, 0);
	return;
}

