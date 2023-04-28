#ifndef _LEAF_VECTOR_DOTPROD_H_
#define	_LEAF_VECTOR_DOTPROD_H_ 

#include <Rdefines.h>

double _dotprod_leaf_vectors(
	SEXP lv1,
	SEXP lv2
);

double _dotprod_leaf_vector_and_finite_col(
	SEXP lv1,
	const double *x2
);

double _dotprod_leaf_vector_and_double_col(
	SEXP lv1,
	const double *x2,
	int x2_len
);

double _dotprod_leaf_vector_and_noNA_int_col(
	SEXP lv1,
	const int *x2
);

double _dotprod_leaf_vector_and_int_col(
	SEXP lv1,
	const int *x2,
	int x2_len
);

double _dotprod0_int_col(
	const int *x,
	int x_len
);

double _dotprod0_double_col(
	const double *x,
	int x_len
);

double _dotprod0_leaf_vector(SEXP lv);

#endif  /* _LEAF_VECTOR_DOTPROD_H_ */

