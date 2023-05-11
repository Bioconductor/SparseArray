#ifndef _SPARSEARRAY_DIM_TUNING_H_
#define	_SPARSEARRAY_DIM_TUNING_H_

#include <Rdefines.h>

SEXP C_tune_SVT_dims(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP dim_tuner
);

#endif  /* _SPARSEARRAY_DIM_TUNING_H_ */

