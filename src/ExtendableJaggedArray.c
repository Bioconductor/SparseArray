/****************************************************************************
 *          A simple implementation of an Extendable Jagged Array           *
 ****************************************************************************/
#include "ExtendableJaggedArray.h"

#include "S4Vectors_interface.h"

#include "leaf_utils.h"

#include <stdlib.h>  /* for malloc(), free(), realloc() */
#include <string.h>  /* for memcpy() */

ExtendableJaggedArray _new_ExtendableJaggedArray(size_t ncol)
{
	ExtendableJaggedArray x;
	size_t j;

	x._ncol = ncol;
	x._cols = (int **) malloc(sizeof(int *) * ncol);
	if (x._cols == NULL)
		goto on_error;
	x._buflengths = (size_t *) malloc(sizeof(size_t) * ncol);
	if (x._buflengths == NULL) {
		free(x._cols);
		goto on_error;
	}
	x._nelts = (size_t *) malloc(sizeof(size_t) * ncol);
	if (x._nelts == NULL) {
		free(x._buflengths);
		free(x._cols);
		goto on_error;
	}
	for (j = 0; j < ncol; j++)
		x._buflengths[j] = x._nelts[j] = 0;
	return x;

    on_error:
	error("SparseArray internal error in "
	      "_new_ExtendableJaggedArray():\n"
	      "    memory allocation failed");
}

void _free_ExtendableJaggedArray(ExtendableJaggedArray *x)
{
	size_t j;

	for (j = 0; j < x->_ncol; j++) {
		if (x->_buflengths[j] != 0)
			free(x->_cols[j]);
	}
	free(x->_nelts);
	free(x->_buflengths);
	free(x->_cols);
	return;
}

static void extend_ExtendableJaggedArray_col(ExtendableJaggedArray *x, int j)
{
	size_t current_buflength, new_buflength, new_size;
	int *col;

	current_buflength = x->_buflengths[j];
	new_buflength = increase_buflength(current_buflength);
	new_size = sizeof(int) * new_buflength;
	if (current_buflength == 0) {
		col = (int *) malloc(new_size);
		if (col == NULL) {
			_free_ExtendableJaggedArray(x);
			error("SparseArray internal error in "
			      "extend_ExtendableJaggedArray_col():\n"
			      "    memory allocation failed");
		}
	} else {
		col = (int *) realloc(x->_cols[j], new_size);
		if (col == NULL) {
			_free_ExtendableJaggedArray(x);
			error("SparseArray internal error in "
			      "extend_ExtendableJaggedArray_col():\n"
			      "    memory reallocation failed");
		}
	}
	x->_cols[j] = col;
	x->_buflengths[j] = new_buflength;
	return;
}

void _add_ExtendableJaggedArray_elt(ExtendableJaggedArray *x,
				    int j, int val)
{
	if (x->_nelts[j] == x->_buflengths[j])
		extend_ExtendableJaggedArray_col(x, j);
	x->_cols[j][(x->_nelts[j])++] = val;
	return;
}

static SEXP make_leaf_vector(const int *offs, const int *vals, int lv_len)
{
	SEXP ans_offs, ans_vals, ans;

	ans_offs = PROTECT(NEW_INTEGER(lv_len));
	memcpy(INTEGER(ans_offs), offs, sizeof(int) * lv_len);
	ans_vals = PROTECT(NEW_INTEGER(lv_len));
	memcpy(INTEGER(ans_vals), vals, sizeof(int) * lv_len);
	ans = zip_leaf(ans_offs, ans_vals);  // unprotected!
	UNPROTECT(2);
	return ans;
}

/* 'offss' and 'valss' are **asumed** to have the same shape but we don't
   check this!
   The function frees the columns in 'offss' and 'valss' as it walks over
   them and copies their content to the SVT. Note that it's still the
   responsibility of the caller to call _free_ExtendableJaggedArray() on
   'offss' and 'valss' on return. */
SEXP _move_ExtendableJaggedArrays_to_SVT(ExtendableJaggedArray *offss,
					 ExtendableJaggedArray *valss)
{
	int SVT_len, is_empty, i, lv_len;
	SEXP ans, ans_elt;
	int *offs, *vals;

	SVT_len = offss->_ncol;
	ans = PROTECT(NEW_LIST(SVT_len));
	is_empty = 1;
	for (i = 0; i < SVT_len; i++) {
		lv_len = offss->_nelts[i];  // assumed to be the same
					    // as 'valss->_nelts[i]'
		if (lv_len != 0) {
			offs = offss->_cols[i];
			vals = valss->_cols[i];
			ans_elt = make_leaf_vector(offs, vals, lv_len);
			PROTECT(ans_elt);
			SET_VECTOR_ELT(ans, i, ans_elt);
			UNPROTECT(1);
			is_empty = 0;
		}
		if (offss->_buflengths[i] != 0) {
			free(offs);
			offss->_buflengths[i] = offss->_nelts[i] = 0;
		}
		if (valss->_buflengths[i] != 0) {
			free(vals);
			valss->_buflengths[i] = valss->_nelts[i] = 0;
		}
	}
	UNPROTECT(1);
	return is_empty ? R_NilValue : ans;
}

