#ifndef _SPARSEMATRIX_UTILS_H_
#define _SPARSEMATRIX_UTILS_H_

#include <Rdefines.h>

SEXP C_rowsum_dgCMatrix(SEXP x, SEXP group, SEXP ngroup, SEXP na_rm);
SEXP C_colMins_dgCMatrix(SEXP x, SEXP na_rm);
SEXP C_colMaxs_dgCMatrix(SEXP x, SEXP na_rm);
SEXP C_colRanges_dgCMatrix(SEXP x, SEXP na_rm);
SEXP C_colVars_dgCMatrix(SEXP x, SEXP na_rm);

#endif  /* _SPARSEMATRIX_UTILS_H_ */

