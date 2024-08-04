#ifndef _COERCEVECTOR2_H_
#define _COERCEVECTOR2_H_

#include <Rdefines.h>

void _CoercionWarning(int warn);

int _coercion_can_introduce_zeros(
	SEXPTYPE from_Rtype,
	SEXPTYPE to_Rtype
);

SEXP C_coercion_can_introduce_zeros(
	SEXP from_type,
	SEXP to_type
);

int _coercion_can_introduce_NAs(
	SEXPTYPE from_Rtype,
	SEXPTYPE to_Rtype
);

SEXP C_coercion_can_introduce_NAs(
	SEXP from_type,
	SEXP to_type
);

SEXP _coerceVector2(
	SEXP v,
	SEXPTYPE type,
	int *warn
);

#endif  /* _COERCEVECTOR2_H_ */

