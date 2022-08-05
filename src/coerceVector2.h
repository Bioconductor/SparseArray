#ifndef _COERCEVECTOR2_H_
#define _COERCEVECTOR2_H_

#include <Rdefines.h>

void _CoercionWarning(int warn);

SEXP _coerceVector2(
	SEXP v,
	SEXPTYPE type,
	int *warn
);

#endif  /* _COERCEVECTOR2_H_ */

