#include <R_ext/Rdynload.h>

#include "coerceVector2.h"
#include "OPBufTree.h"
#include "thread_control.h"
#include "leaf_utils.h"
#include "sparseMatrix_utils.h"
#include "SVT_SparseArray_class.h"
#include "SparseArray_dim_tuning.h"
#include "SparseArray_aperm.h"
#include "SparseArray_subsetting.h"
#include "SparseArray_subassignment.h"
#include "SparseArray_subassignment_OLD.h"
#include "SparseArray_abind.h"
#include "SparseArray_summarization.h"
#include "SparseArray_Arith_methods.h"
#include "SparseArray_Compare_methods.h"
#include "SparseArray_Logic_methods.h"
#include "SparseArray_Math_methods.h"
#include "SparseArray_Complex_methods.h"
#include "SparseArray_misc_methods.h"
#include "SparseArray_matrixStats.h"
#include "rowsum_methods.h"
#include "SparseMatrix_mult.h"
#include "randomSparseArray.h"
#include "readSparseCSV.h"
#include "test.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* coerceVector2.c */
	CALLMETHOD_DEF(C_coercion_can_introduce_zeros, 2),
	CALLMETHOD_DEF(C_coercion_can_introduce_NAs, 2),

/* OPBufTree.c */
	CALLMETHOD_DEF(C_free_global_OPBufTree, 0),

/* thread_control.c */
	CALLMETHOD_DEF(C_get_num_procs, 0),
	CALLMETHOD_DEF(C_get_max_threads, 0),
	CALLMETHOD_DEF(C_set_max_threads, 1),

/* leaf_utils.c */
	CALLMETHOD_DEF(C_lacunar_mode_is_on, 0),

/* sparseMatrix_utils.c */
	CALLMETHOD_DEF(C_colMins_dgCMatrix, 2),
	CALLMETHOD_DEF(C_colMaxs_dgCMatrix, 2),
	CALLMETHOD_DEF(C_colRanges_dgCMatrix, 2),
	CALLMETHOD_DEF(C_colVars_dgCMatrix, 2),

/* SVT_SparseArray_class.c */
	CALLMETHOD_DEF(C_set_SVT_SparseArray_type, 5),
	CALLMETHOD_DEF(C_nzcount_SVT_SparseArray, 2),
	CALLMETHOD_DEF(C_nzwhich_SVT_SparseArray, 3),
	CALLMETHOD_DEF(C_from_SVT_SparseArray_to_Rarray, 5),
	CALLMETHOD_DEF(C_build_SVT_from_Rarray, 3),
	CALLMETHOD_DEF(C_from_SVT_SparseMatrix_to_CsparseMatrix, 4),
	CALLMETHOD_DEF(C_build_SVT_from_CSC, 5),
	CALLMETHOD_DEF(C_build_SVT_from_CsparseMatrix, 2),
	CALLMETHOD_DEF(C_from_SVT_SparseArray_to_COO_SparseArray, 3),

/* SparseArray_dim_tuning.c */
	CALLMETHOD_DEF(C_tune_SVT_dims, 4),

/* SparseArray_aperm.c */
	CALLMETHOD_DEF(C_transpose_2D_SVT, 3),
	CALLMETHOD_DEF(C_aperm0_SVT, 4),
	CALLMETHOD_DEF(C_aperm_SVT, 4),

/* SparseArray_subsetting.c */
	CALLMETHOD_DEF(C_subset_SVT_by_Lindex, 4),
	CALLMETHOD_DEF(C_subset_SVT_by_Mindex, 4),
	CALLMETHOD_DEF(C_subset_SVT_by_Nindex, 4),

/* SparseArray_subassignment.c */
	CALLMETHOD_DEF(C_subassign_SVT_by_Lindex, 6),
	CALLMETHOD_DEF(C_subassign_SVT_by_Mindex, 5),
	CALLMETHOD_DEF(C_subassign_SVT_with_short_Rvector, 5),
	CALLMETHOD_DEF(C_subassign_SVT_with_Rarray, 5),
	CALLMETHOD_DEF(C_subassign_SVT_with_SVT, 7),

/* SparseArray_subassignment_OLD.c */
	CALLMETHOD_DEF(C_subassign_SVT_by_Lindex_OLD, 5),
	CALLMETHOD_DEF(C_subassign_SVT_by_Mindex_OLD, 5),

/* SparseArray_abind.c */
	CALLMETHOD_DEF(C_abind_SVT_SparseArray_objects, 4),

/* SparseArray_summarization.c */
	CALLMETHOD_DEF(C_summarize_SVT, 7),

/* SparseArray_Arith_methods.c */
	CALLMETHOD_DEF(C_unary_minus_SVT, 3),
	CALLMETHOD_DEF(C_Arith_SVT1_v2, 6),
	CALLMETHOD_DEF(C_Arith_SVT1_SVT2, 8),

/* SparseArray_Compare_methods.c */
	CALLMETHOD_DEF(C_Compare_SVT1_v2, 6),
	CALLMETHOD_DEF(C_Compare_SVT1_SVT2, 7),

/* SparseArray_Logic_methods.c */
	CALLMETHOD_DEF(C_Logic_SVT1_SVT2, 7),

/* SparseArray_Math_methods.c */
	CALLMETHOD_DEF(C_Math_SVT, 5),

/* SparseArray_Complex_methods.c */
	CALLMETHOD_DEF(C_Complex_SVT, 4),

/* SparseArray_misc_methods.c */
	CALLMETHOD_DEF(C_SVT_apply_isFUN, 4),

/* SparseArray_matrixStats.c */
	CALLMETHOD_DEF(C_colStats_SVT, 9),
	CALLMETHOD_DEF(C_rowStats_SVT, 9),

/* rowsum_methods.c */
	CALLMETHOD_DEF(C_rowsum_SVT, 6),
	CALLMETHOD_DEF(C_rowsum_dgCMatrix, 4),
	CALLMETHOD_DEF(C_colsum_SVT, 6),
	CALLMETHOD_DEF(C_colsum_dgCMatrix, 4),

/* SparseMatrix_mult.c */
	CALLMETHOD_DEF(C_crossprod2_SVT_mat, 7),
	CALLMETHOD_DEF(C_crossprod2_mat_SVT, 7),
	CALLMETHOD_DEF(C_crossprod2_SVT_SVT, 8),
	CALLMETHOD_DEF(C_crossprod1_SVT, 5),

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

