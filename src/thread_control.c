/****************************************************************************
 *               Wrappers for OpenMP thread control functions               *
 ****************************************************************************/
#include "thread_control.h"

#ifdef _OPENMP
/* <Rinternals.h> defines macro match that seems to break <omp.h> with
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

int _get_max_threads(void)
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

static int get_initial_device(void)
{
#ifdef _OPENMP
	return omp_get_initial_device();
#else
	return 0;
#endif
}

static int pause_resource(int hard_pause, int device_num)
{
#ifdef _OPENMP
	omp_pause_resource_t kind = hard_pause ? omp_pause_hard :
						 omp_pause_soft;
	return omp_pause_resource(kind, device_num);
#else
	return 0;
#endif
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
	return ScalarInteger(_get_max_threads());
}

/* --- .Call ENTRY POINT --- */
SEXP C_set_max_threads(SEXP nthread)
{
	int prev_max_threads = _get_max_threads();
	set_max_threads(INTEGER(nthread)[0]);
	return ScalarInteger(prev_max_threads);
}

/* --- .Call ENTRY POINT --- */
SEXP C_get_initial_device(void)
{
	return ScalarInteger(get_initial_device());
}

/* --- .Call ENTRY POINT --- */
SEXP C_pause_resource(SEXP hard_pause, SEXP device_num)
{
	int h = LOGICAL(hard_pause)[0], d = INTEGER(device_num)[0];
	return ScalarInteger(pause_resource(h, d));
}

