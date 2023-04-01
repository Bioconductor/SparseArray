#include <R_ext/Rdynload.h>

#include "sparseMatrix_utils.h"
#include "SVT_SparseArray_class.h"
#include "SparseArray_subassignment.h"
#include "SparseArray_subsetting.h"
#include "SparseArray_combine.h"
#include "SparseArray_summarization.h"
#include "SparseMatrix_mult.h"
#include "randomSparseArray.h"
#include "readSparseCSV.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* sparseMatrix_utils.c */
	CALLMETHOD_DEF(C_rowsum_dgCMatrix, 4),
	CALLMETHOD_DEF(C_colMins_dgCMatrix, 2),
	CALLMETHOD_DEF(C_colMaxs_dgCMatrix, 2),
	CALLMETHOD_DEF(C_colRanges_dgCMatrix, 2),
	CALLMETHOD_DEF(C_colVars_dgCMatrix, 2),

/* SVT_SparseArray_class.c */
	CALLMETHOD_DEF(C_get_SVT_SparseArray_nzcount, 2),
	CALLMETHOD_DEF(C_set_SVT_SparseArray_type, 4),
	CALLMETHOD_DEF(C_from_SVT_SparseArray_to_Rarray, 4),
	CALLMETHOD_DEF(C_build_SVT_from_Rarray, 2),
	CALLMETHOD_DEF(C_from_SVT_SparseMatrix_to_CsparseMatrix, 3),
	CALLMETHOD_DEF(C_build_SVT_from_CsparseMatrix, 2),
	CALLMETHOD_DEF(C_from_SVT_SparseArray_to_COO_SparseArray, 3),
	CALLMETHOD_DEF(C_transpose_SVT_SparseMatrix, 3),

/* SparseArray_subassignment.c */
	CALLMETHOD_DEF(C_subassign_SVT_by_Mindex, 5),
	CALLMETHOD_DEF(C_subassign_SVT_by_Lindex, 5),
	CALLMETHOD_DEF(C_subassign_SVT_with_short_Rvector, 5),
	CALLMETHOD_DEF(C_subassign_SVT_with_Rarray, 5),
	CALLMETHOD_DEF(C_subassign_SVT_with_SVT, 7),

/* SparseArray_subsetting.c */
	CALLMETHOD_DEF(C_drop_SVT_SparseArray_ineffective_dims, 4),
	CALLMETHOD_DEF(C_subset_SVT_SparseArray, 4),

/* SparseArray_combine.c */
	CALLMETHOD_DEF(C_abind_SVT_SparseArray_objects, 3),

/* SparseArray_summarization.c */
	CALLMETHOD_DEF(C_summarize_SVT_SparseArray, 6),
	CALLMETHOD_DEF(C_count_SVT_SparseArray_NAs, 3),
	CALLMETHOD_DEF(C_anyNA_SVT_SparseArray, 3),

/* SparseMatrix_mult.c */
	CALLMETHOD_DEF(C_SVT_SparseMatrix_crossprod, 5),

/* randomSparseArray.c */
	CALLMETHOD_DEF(C_simple_rpois, 2),
	CALLMETHOD_DEF(C_poissonSparseArray, 2),

/* readSparseCSV.c */
	CALLMETHOD_DEF(C_readSparseCSV_as_SVT_SparseMatrix, 5),

	{NULL, NULL, 0}
};

void R_init_SparseArray(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	return;
}

