#ifndef _MISC_UTILS_H_
#define _MISC_UTILS_H_

#include <Rdefines.h>  /* for R_xlen_t */

R_xlen_t *_alloc_and_compute_cumprod(
	const int *x,
	int x_len
);

int _which_max(
	const int *x,
	int x_len
);

#endif  /* _MISC_UTILS_H_ */

