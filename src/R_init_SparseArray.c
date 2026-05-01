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


/* Initialize global variables character0 and character1 declared in
   src/Rvector_utils.h */
static SEXP C_init_character0_character1(SEXP strings01)
{
	const char *errmsg = "SparseArray internal error in "
			     "C_init_character0_character1():\n"
			     "    'strings01' must be 'c(\"\", \"1\")'";
	if (!(IS_CHARACTER(strings01) && LENGTH(strings01) == 2))
		error("%s", errmsg);
	character0 = STRING_ELT(strings01, 0);  /* CHARSXP */
	character1 = STRING_ELT(strings01, 1);  /* CHARSXP */
	/* Some sanity checks.
	   Note that we're comparing the CHARSXPs' addresses, not their values.
	   However, this is much faster, but also, and most importantly, it's
	   equivalent to comparing their values. That's because CHARSXPs with
	   the same value are expected to have the same address, thanks to R's
	   global CHARSXP cache.
	   Also, because of this caching, there's no need to PROTECT() the
	   address returned by mkChar(). */
	if (character0 != mkChar("") || character1 != mkChar("1"))
		error("%s", errmsg);
	return R_NilValue;
}

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
	CALLMETHOD_DEF(C_get_initial_device, 0),
	CALLMETHOD_DEF(C_pause_resource, 2),

/* sparseMatrix_utils.c */
	CALLMETHOD_DEF(C_colMins_dgCMatrix, 2),
	CALLMETHOD_DEF(C_colMaxs_dgCMatrix, 2),
	CALLMETHOD_DEF(C_colRanges_dgCMatrix, 2),
	CALLMETHOD_DEF(C_colVars_dgCMatrix, 2),

/* SVT_SparseArray_class.c */
	CALLMETHOD_DEF(C_set_SVT_type, 5),
	CALLMETHOD_DEF(C_is_nonzero_SVT, 2),
	CALLMETHOD_DEF(C_nzcount_SVT, 2),
	CALLMETHOD_DEF(C_nzwhich_SVT, 3),
	CALLMETHOD_DEF(C_nzvals_SVT, 3),
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
	CALLMETHOD_DEF(C_subset_SVT_by_Lindex, 5),
	CALLMETHOD_DEF(C_subset_SVT_by_Mindex, 5),
	CALLMETHOD_DEF(C_subset_SVT_as_Rarray, 5),
	CALLMETHOD_DEF(C_subset_SVT_as_SVT, 4),

/* SparseArray_subassignment.c */
	CALLMETHOD_DEF(C_subassign_SVT_by_Lindex, 6),
	CALLMETHOD_DEF(C_subassign_SVT_by_Mindex, 6),
	CALLMETHOD_DEF(C_subassign_SVT_with_short_Rvector, 6),
	CALLMETHOD_DEF(C_subassign_SVT_with_Rarray, 6),
	CALLMETHOD_DEF(C_subassign_SVT_with_SVT, 9),

/* SparseArray_abind.c */
	CALLMETHOD_DEF(C_abind_SVT_SparseArray_objects, 4),

/* SparseArray_summarization.c */
	CALLMETHOD_DEF(C_summarize_SVT, 7),

/* SparseArray_Arith_methods.c */
	CALLMETHOD_DEF(C_unary_minus_SVT, 3),
	CALLMETHOD_DEF(C_Arith_SVT1_v2, 8),
	CALLMETHOD_DEF(C_Arith_v1_SVT2, 7),
	CALLMETHOD_DEF(C_Arith_SVT1_SVT2, 10),

/* SparseArray_Compare_methods.c */
	CALLMETHOD_DEF(C_Compare_SVT1_v2, 6),
	CALLMETHOD_DEF(C_Compare_SVT1_SVT2, 9),

/* SparseArray_Logic_methods.c */
	CALLMETHOD_DEF(C_logical_neg_NaSVT, 3),
	CALLMETHOD_DEF(C_Logic_NaSVT1_na, 4),
	CALLMETHOD_DEF(C_Logic_SVT1_SVT2, 9),

/* SparseArray_Math_methods.c */
	CALLMETHOD_DEF(C_Math_SVT, 6),

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
	CALLMETHOD_DEF(C_simple_omp_parallel_for_loop, 1),

/* this file */
	CALLMETHOD_DEF(C_init_character0_character1, 1),

	{NULL, NULL, 0}
};

void R_init_SparseArray(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, 0);

	/* Initialize global variables declared in src/Rvector_utils.h */
	intNA = NA_INTEGER;
	doubleNA = RcomplexNA.r = RcomplexNA.i = NA_REAL;
	characterNA = NA_STRING;            /* CHARSXP */
	list0 = R_NilValue;
	return;
}

