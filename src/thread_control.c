/****************************************************************************
 *               Wrappers for OpenMP thread control functions               *
 ****************************************************************************/
#include "thread_control.h"

#ifdef _OPENMP
/* <Rinternals.h> defines macro match that seems to break <omp.h> on
   some versions of Clang.
   See https://github.com/Bioconductor/SparseArray/issues/9 */
#undef match
#include <omp.h>
#endif


static int get_num_procs(void)
{
#ifdef _OPENMP
	return omp_get_num_procs();
#else
	return 0;
#endif
}

static int get_max_threads(void)
{
#ifdef _OPENMP
	return omp_get_max_threads();
#else
	return 0;
#endif
}

static void set_max_threads(int nthread)
{
#ifdef _OPENMP
	omp_set_num_threads(nthread);
#endif
	return;
}


/****************************************************************************
 * .Call ENTRY POINTS
 */

/* --- .Call ENTRY POINT --- */
SEXP C_get_num_procs(void)
{
	return ScalarInteger(get_num_procs());
}

/* --- .Call ENTRY POINT --- */
SEXP C_get_max_threads(void)
{
	return ScalarInteger(get_max_threads());
}

/* --- .Call ENTRY POINT --- */
SEXP C_set_max_threads(SEXP nthread)
{
	int prev_max_threads = get_max_threads();
	set_max_threads(INTEGER(nthread)[0]);
	return ScalarInteger(prev_max_threads);
}

