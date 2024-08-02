#ifndef _SPARSEARRAY_ABIND_H_
#define _SPARSEARRAY_ABIND_H_

#include <Rdefines.h>

SEXP C_abind_SVT_SparseArray_objects(
	SEXP objects,
	SEXP SVTslotname,
	SEXP along,
	SEXP ans_type
);

#endif  /* _SPARSEARRAY_ABIND_H_ */

