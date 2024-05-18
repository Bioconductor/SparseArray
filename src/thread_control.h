#ifndef _THREAD_CONTROL_H_
#define _THREAD_CONTROL_H_

#include <Rdefines.h>


/* Used in C_colStats_SVT(), C_subset_SVT_by_Lindex(), etc.. to find
   the dimension along which to perform parallel execution.
   Returns the 0-based position of the max value in 'x' (greatest position
   if max is found at more than one position in 'x'). */
static inline int which_max(const int *x, int x_len)
{
	if (x_len == 0)
		return -1;
	int max_idx = x_len - 1;
	for (int i = max_idx - 1; i >= 0; i--) {
		if (x[i] > x[max_idx])
			max_idx = i;
	}
	return max_idx;
}

SEXP C_get_num_procs(void);

SEXP C_get_max_threads(void);

SEXP C_set_max_threads(SEXP nthread);

#endif  /* _THREAD_CONTROL_H_ */

