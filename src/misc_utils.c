/****************************************************************************
 *                     Miscellaneous utility functions                      *
 ****************************************************************************/
#include "misc_utils.h"


R_xlen_t *_alloc_and_compute_cumprod(const int *x, int x_len)
{
	R_xlen_t *cumprod = (R_xlen_t *) R_alloc(x_len, sizeof(R_xlen_t));
	R_xlen_t prod = 1;
	for (int i = 0; i < x_len; i++) {
		prod *= x[i];
		cumprod[i] = prod;
	}
	return cumprod;
}

/* Used in C_colStats_SVT(), C_subset_SVT_by_Lindex(), etc.. to find
   the dimension along which to perform parallel execution.
   Returns the 0-based position of the max value in 'x' (greatest position
   if max is found at more than one position in 'x'). */
int _which_max(const int *x, int x_len)
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

