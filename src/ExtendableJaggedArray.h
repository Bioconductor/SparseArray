#ifndef _EXTENDABLE_JAGGED_ARRAY_H_
#define _EXTENDABLE_JAGGED_ARRAY_H_

#include <Rdefines.h>

typedef struct extendable_jagged_array {
	size_t _ncol;
	int **_cols;
	size_t *_buflengths;
	size_t *_nelts;
} ExtendableJaggedArray;

ExtendableJaggedArray _new_ExtendableJaggedArray(size_t ncol);

void _free_ExtendableJaggedArray(ExtendableJaggedArray *x);

void _add_ExtendableJaggedArray_elt(
	ExtendableJaggedArray *x,
	int j,
	int val
);

SEXP _move_ExtendableJaggedArrays_to_SVT(
	ExtendableJaggedArray *offss,
	ExtendableJaggedArray *valss
);

#endif  /* _EXTENDABLE_JAGGED_ARRAY_H_ */
