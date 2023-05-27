#ifndef _SBT_UTILS_H_
#define	_SBT_UTILS_H_

#include <Rdefines.h>

/* A "Sparse Buffer Tree" is a "Sparse Vector Tree" where the leaves
   are "leaf buffers" instead of "leaf vectors". */

void _push_int_to_SBT(
	SEXP SBT,
	const int *dim,
	int ndim,
	const int *coords0,
	int val
);

void _push_double_to_SBT(
	SEXP SBT,
	const int *dim,
	int ndim,
	const int *coords0,
	double val
);

void _push_Rcomplex_to_SBT(
	SEXP SBT,
	const int *dim,
	int ndim,
	const int *coords0,
	Rcomplex val
);

void _push_Rbyte_to_SBT(
	SEXP SBT,
	const int *dim,
	int ndim,
	const int *coords0,
	Rbyte val
);

void _push_SEXP_to_SBT(
	SEXP SBT,
	const int *dim,
	int ndim,
	const int *coords0,
	SEXP val
);

void _SBT2SVT(
	SEXP SBT,
	const int *dim,
	int ndim,
	SEXPTYPE Rtype
);

#endif  /* _SBT_UTILS_H_ */

