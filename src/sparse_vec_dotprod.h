#ifndef _SPARSE_VEC_DOTPROD_H_
#define _SPARSE_VEC_DOTPROD_H_

#include "sparse_vec.h"

double _dotprod_sparse_vecs(
	const struct sparse_vec *sv1,
	const struct sparse_vec *sv2
);

double _dotprod_sparse_vec_and_finite_col(
	const struct sparse_vec *sv1,
	const double *x2
);

double _dotprod_sparse_vec_and_double_col(
	const struct sparse_vec *sv1,
	const double *x2
);

double _dotprod_sparse_vec_and_noNA_int_col(
	const struct sparse_vec *sv1,
	const int *x2
);

double _dotprod_sparse_vec_and_int_col(
	const struct sparse_vec *sv1,
	const int *x2
);

double _dotprod0_int_col(
	const int *x,
	int x_len
);

double _dotprod0_double_col(
	const double *x,
	int x_len
);

double _dotprod0_sparse_vec(const struct sparse_vec *sv);

#endif  /* _SPARSE_VEC_DOTPROD_H_ */

