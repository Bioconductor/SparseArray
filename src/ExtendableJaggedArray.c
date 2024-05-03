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

static SEXP make_leaf(const int *nzoffs, const int *nzvals, int nzcount)
{
	SEXP ans_nzoffs, ans_nzvals;
	SEXP ans = PROTECT(_alloc_and_unzip_leaf(INTSXP, nzcount,
						 &ans_nzoffs, &ans_nzvals));
	memcpy(INTEGER(ans_nzoffs), nzoffs, sizeof(int) * nzcount);
	memcpy(INTEGER(ans_nzvals), nzvals, sizeof(int) * nzcount);
	UNPROTECT(1);
	return ans;
}

/* 'nzoffss' and 'nzvalss' are **asumed** to have the same shape but we
   don't check this!
   The function frees the columns in 'nzoffss' and 'nzvalss' as it walks
   over them and copies their content to the SVT. Note that it's still the
   responsibility of the caller to call _free_ExtendableJaggedArray() on
   'nzoffss' and 'nzvalss' on return. */
SEXP _move_ExtendableJaggedArrays_to_SVT(ExtendableJaggedArray *nzoffss,
					 ExtendableJaggedArray *nzvalss)
{
	int SVT_len = nzoffss->_ncol;
	SEXP ans = PROTECT(NEW_LIST(SVT_len));
	int is_empty = 1;
	for (int i = 0; i < SVT_len; i++) {
		int nzcount = nzoffss->_nelts[i];  // assumed to be the same
						   // as 'nzvalss->_nelts[i]'
		int *nzoffs, *nzvals;
		if (nzcount != 0) {
			nzoffs = nzoffss->_cols[i];
			nzvals = nzvalss->_cols[i];
			SEXP ans_elt =
				PROTECT(make_leaf(nzoffs, nzvals, nzcount));
			SET_VECTOR_ELT(ans, i, ans_elt);
			UNPROTECT(1);
			is_empty = 0;
		}
		if (nzoffss->_buflengths[i] != 0) {
			free(nzoffs);
			nzoffss->_buflengths[i] = nzoffss->_nelts[i] = 0;
		}
		if (nzvalss->_buflengths[i] != 0) {
			free(nzvals);
			nzvalss->_buflengths[i] = nzvalss->_nelts[i] = 0;
		}
	}
	UNPROTECT(1);
	return is_empty ? R_NilValue : ans;
}

