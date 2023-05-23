/****************************************************************************
 *                            Sparse Buffer Tree                            *
 ****************************************************************************/
#include "SBT_utils.h"

#include "Rvector_utils.h"
#include "leaf_vector_utils.h"

#include <stdlib.h>  /* for malloc(), free(), realloc() */
#include <limits.h>  /* for INT_MAX */
#include <string.h>  /* for memcpy() */


/****************************************************************************
 * Manipulation of SparseBuf structs
 */

typedef union vals_buf_t {
	int *ints;
	double *doubles;
	Rcomplex *Rcomplexes;
	Rbyte *Rbytes;
	SEXP *SEXPs;
} ValsBuf;

typedef struct sparse_buf_t {
	int buflength;
	int nelt;
	int *offs;
	ValsBuf vals;
} SparseBuf;

/* Guaranteed to return a value > 'buflength', or to raise an error. */
#define	MAX_BUFLENGTH_INC 16777216  /* 2^24 */
static int increase_buflength(int buflength)
{
	unsigned int new_len;

	if (buflength == INT_MAX)
		error("SparseArray internal error in increase_buflength(): "
		      "max buflength reached");
	if (buflength <= 4)
		return 8;
	if (buflength <= 8)
		return 32;
	if (buflength <= 32)
		return 128;
	if (buflength <= MAX_BUFLENGTH_INC)
		return 2 * buflength;
	new_len = buflength + MAX_BUFLENGTH_INC;
	if (new_len > INT_MAX)
		new_len = INT_MAX;
	return new_len;
}

#define	FUNDEF_new_SparseBuf(type, union_member)			     \
	(int buflength)							     \
{									     \
	SparseBuf *buf;							     \
									     \
	buf = (SparseBuf *) malloc(sizeof(SparseBuf));			     \
	if (buf == NULL)						     \
		error("new_" #type "_SparseBuf: malloc() error");	     \
	buf->offs = (int *) malloc(sizeof(int) * buflength);		     \
	if (buf->offs == NULL) {					     \
		free(buf);						     \
		error("new_" #type "_SparseBuf: malloc() error");	     \
	}								     \
	buf->vals.union_member = (type *) malloc(sizeof(type) * buflength);  \
	if (buf->vals.union_member == NULL) {				     \
		free(buf->offs);					     \
		free(buf);						     \
		error("new_" #type "_SparseBuf: malloc() error");	     \
	}								     \
	buf->buflength = buflength;					     \
	buf->nelt = 0;							     \
	return buf;							     \
}

static SparseBuf *new_int_SparseBuf
	FUNDEF_new_SparseBuf(int, ints)
static SparseBuf *new_double_SparseBuf
	FUNDEF_new_SparseBuf(double, doubles)
static SparseBuf *new_Rcomplex_SparseBuf
	FUNDEF_new_SparseBuf(Rcomplex, Rcomplexes)
static SparseBuf *new_Rbyte_SparseBuf
	FUNDEF_new_SparseBuf(Rbyte, Rbytes)
static SparseBuf *new_SEXP_SparseBuf
	FUNDEF_new_SparseBuf(SEXP, SEXPs)

#define	FUNDEF_extend_SparseBuf(type, union_member)			     \
	(SparseBuf *buf, int new_buflength)				     \
{									     \
	int *new_offs;							     \
	type *new_vals;							     \
									     \
	new_offs = (int *) realloc(buf->offs, sizeof(int) * new_buflength);  \
	if (new_offs == NULL)						     \
		error("extend_" #type "_SparseBuf: realloc() error");	     \
	buf->offs = new_offs;						     \
	new_vals = (type *) realloc(buf->vals.union_member,		     \
				    sizeof(type) * new_buflength);	     \
	if (new_vals == NULL)						     \
		error("extend_" #type "_SparseBuf: realloc() error");	     \
	buf->vals.union_member = new_vals;				     \
	buf->buflength = new_buflength;					     \
	return;								     \
}

static void extend_int_SparseBuf
	FUNDEF_extend_SparseBuf(int, ints)
static void extend_double_SparseBuf
	FUNDEF_extend_SparseBuf(double, doubles)
static void extend_Rcomplex_SparseBuf
	FUNDEF_extend_SparseBuf(Rcomplex, Rcomplexes)
static void extend_Rbyte_SparseBuf
	FUNDEF_extend_SparseBuf(Rbyte, Rbytes)
static void extend_SEXP_SparseBuf
	FUNDEF_extend_SparseBuf(SEXP, SEXPs)

#define	FUNDEF_free_SparseBuf(union_member)				     \
	(SparseBuf *buf)						     \
{									     \
	/*printf("free SparseBuf -- ");*/				     \
	/*printf("nelt/buflength = %d/%d\n", buf->nelt, buf->buflength);*/   \
	free(buf->vals.union_member);					     \
	free(buf->offs);						     \
	free(buf);							     \
	return;								     \
}

static void free_int_SparseBuf FUNDEF_free_SparseBuf(ints)
static void free_double_SparseBuf FUNDEF_free_SparseBuf(doubles)
static void free_Rcomplex_SparseBuf FUNDEF_free_SparseBuf(Rcomplexes)
static void free_Rbyte_SparseBuf FUNDEF_free_SparseBuf(Rbytes)
static void free_SEXP_SparseBuf FUNDEF_free_SparseBuf(SEXPs)

/* 'buf' can NOT be empty! Returns a "leaf vector". */
static SEXP make_leaf_vector_from_SparseBuf(SEXPTYPE Rtype,
					    const SparseBuf *buf)
{
	int ans_len, k;
	SEXP ans_offs, ans_vals, ans;

	ans_len = buf->nelt;
	ans_offs = PROTECT(NEW_INTEGER(ans_len));
	memcpy(INTEGER(ans_offs), buf->offs, sizeof(int) * ans_len);
	ans_vals = PROTECT(allocVector(Rtype, ans_len));
	switch (Rtype) {
	    case LGLSXP: case INTSXP:
		memcpy(INTEGER(ans_vals), buf->vals.ints,
		       sizeof(int) * ans_len);
		break;
	    case REALSXP:
		memcpy(REAL(ans_vals), buf->vals.doubles,
		       sizeof(double) * ans_len);
		break;
	    case CPLXSXP:
		memcpy(COMPLEX(ans_vals), buf->vals.Rcomplexes,
		       sizeof(Rcomplex) * ans_len);
		break;
	    case RAWSXP:
		memcpy(RAW(ans_vals), buf->vals.Rbytes,
		       sizeof(Rbyte) * ans_len);
		break;
	    case STRSXP:
		for (k = 0; k < ans_len; k++)
			SET_STRING_ELT(ans_vals, k, buf->vals.SEXPs[k]);
		break;
	    case VECSXP:
		for (k = 0; k < ans_len; k++)
			SET_VECTOR_ELT(ans_vals, k, buf->vals.SEXPs[k]);
		break;
	    default:
		error("SparseArray internal error in "
		      "make_leaf_vector_from_SparseBuf():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));
	}
	ans = _new_leaf_vector(ans_offs, ans_vals);
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * Low-level manipulation of "leaf buffers"
 */

#define	FUNDEF_finalize_leaf_buffer(type)				     \
	(SEXP lb)							     \
{									     \
	SparseBuf *buf;							     \
									     \
	buf = R_ExternalPtrAddr(lb);					     \
	/*printf("finalize_" #type "_leaf_buffer(): ");*/		     \
	if (buf == NULL) {						     \
		/*printf("leaf buffer already finalized\n");*/		     \
		return;							     \
	}								     \
	free_ ## type ## _SparseBuf(buf);				     \
	R_SetExternalPtrAddr(lb, NULL);					     \
	return;								     \
}

static void finalize_int_leaf_buffer FUNDEF_finalize_leaf_buffer(int)
static void finalize_double_leaf_buffer FUNDEF_finalize_leaf_buffer(double)
static void finalize_Rcomplex_leaf_buffer FUNDEF_finalize_leaf_buffer(Rcomplex)
static void finalize_Rbyte_leaf_buffer FUNDEF_finalize_leaf_buffer(Rbyte)
static void finalize_SEXP_leaf_buffer FUNDEF_finalize_leaf_buffer(SEXP)

#define	FUNDEF_new_leaf_buffer(type, union_member)			     \
	(int buflength)							     \
{									     \
	SparseBuf *buf;							     \
	SEXP ans;							     \
									     \
	buf = new_ ## type ## _SparseBuf(buflength);			     \
	buf = (SparseBuf *) malloc(sizeof(SparseBuf));			     \
	buf->offs = (int *) malloc(sizeof(int) * buflength);		     \
	buf->vals.union_member = (type *) malloc(sizeof(type) * buflength);  \
	buf->buflength = buflength;					     \
	buf->nelt = 0;							     \
	ans = PROTECT(R_MakeExternalPtr(buf, R_NilValue, R_NilValue));	     \
	R_RegisterCFinalizer(ans, finalize_ ## type ## _leaf_buffer);	     \
	UNPROTECT(1);							     \
	return ans;							     \
}

static SEXP new_int_leaf_buffer
	FUNDEF_new_leaf_buffer(int, ints)
static SEXP new_double_leaf_buffer
	FUNDEF_new_leaf_buffer(double, doubles)
static SEXP new_Rcomplex_leaf_buffer
	FUNDEF_new_leaf_buffer(Rcomplex, Rcomplexes)
static SEXP new_Rbyte_leaf_buffer
	FUNDEF_new_leaf_buffer(Rbyte, Rbytes)
static SEXP new_SEXP_leaf_buffer
	FUNDEF_new_leaf_buffer(SEXP, SEXPs)

#define	FUNDEF_push_val_to_leaf_buffer(type, union_member)		     \
	(SEXP lb, int off, type val)					     \
{									     \
	SparseBuf *buf;							     \
	int new_buflength;						     \
									     \
	buf = R_ExternalPtrAddr(lb);					     \
	if (buf->nelt == buf->buflength) {				     \
		new_buflength = increase_buflength(buf->buflength);	     \
		extend_ ## type ## _SparseBuf(buf, new_buflength);	     \
	}								     \
	buf->offs[buf->nelt] = off;					     \
	buf->vals.union_member[buf->nelt] = val;			     \
	return ++buf->nelt;						     \
}

static inline int push_int_to_leaf_buffer
	FUNDEF_push_val_to_leaf_buffer(int, ints)
static inline int push_double_to_leaf_buffer
	FUNDEF_push_val_to_leaf_buffer(double, doubles)
static inline int push_Rcomplex_to_leaf_buffer
	FUNDEF_push_val_to_leaf_buffer(Rcomplex, Rcomplexes)
static inline int push_Rbyte_to_leaf_buffer
	FUNDEF_push_val_to_leaf_buffer(Rbyte, Rbytes)
static inline int push_SEXP_to_leaf_buffer
	FUNDEF_push_val_to_leaf_buffer(SEXP, SEXPs)

/* Returns a "leaf vector". */
static SEXP lb2lv(SEXP lb, SEXPTYPE Rtype)
{
	return make_leaf_vector_from_SparseBuf(Rtype, R_ExternalPtrAddr(lb));
}


/****************************************************************************
 * Push values to a Sparse Buffer Tree
 */

/* Must be called with 'ndim' >= 2. */
static inline SEXP descend_SBT(SEXP SBT, const int *dim, int ndim,
		const int *coords0,
		SEXP (*new_leaf_buffer_FUN)(int buflength))
{
	SEXP subSBT;
	int along, i;

	along = ndim - 1;
	do {
		i = coords0[along];
		subSBT = VECTOR_ELT(SBT, i);
		if (along == 1)
			break;
		if (subSBT == R_NilValue) {
			subSBT = PROTECT(NEW_LIST(dim[along - 1]));
			SET_VECTOR_ELT(SBT, i, subSBT);
			UNPROTECT(1);
		}
		SBT = subSBT;
		along--;
	} while (1);
	if (subSBT == R_NilValue) {
		subSBT = PROTECT(new_leaf_buffer_FUN(1));
		SET_VECTOR_ELT(SBT, i, subSBT);
		UNPROTECT(1);
	}
	return subSBT;
}

#define	FUNDEF_push_val_to_SBT(type)					     \
	(SEXP SBT, const int *dim, int ndim, const int *coords0, type val)   \
{									     \
	if (ndim >= 2)							     \
		SBT = descend_SBT(SBT, dim, ndim, coords0,		     \
				  new_ ## type ##_leaf_buffer);		     \
	push_ ## type ## _to_leaf_buffer(SBT, coords0[0], val);		     \
	return;								     \
}


/* For all the _push* functions below, 'SBT' can NOT be R_NilValue! It must
   be a "leaf buffer" (if 'ndim' == 1) or a list (if 'ndim' >= 2). */

/* For SBT of type "logical" or "integer". */
void _push_int_to_SBT		FUNDEF_push_val_to_SBT(int)
/* For SBT of type "double". */
void _push_double_to_SBT	FUNDEF_push_val_to_SBT(double)
/* For SBT of type "complex". */
void _push_Rcomplex_to_SBT	FUNDEF_push_val_to_SBT(Rcomplex)
/* For SBT of type "raw". */
void _push_Rbyte_to_SBT		FUNDEF_push_val_to_SBT(Rbyte)
/* For SBT of type "character" or "list". */
void _push_SEXP_to_SBT		FUNDEF_push_val_to_SBT(SEXP)


/****************************************************************************
 * _SBT2SVT()
 */

/* Recursive. 'SBT' can NOT be R_NilValue! */
static void REC_SBT2SVT(SEXP SBT, const int *dim, int ndim, SEXPTYPE Rtype)
{
	int SBT_len, i;
	SEXP subSBT, lv;

	SBT_len = LENGTH(SBT);
	for (i = 0; i < SBT_len; i++) {
		subSBT = VECTOR_ELT(SBT, i);
		if (subSBT == R_NilValue)
			continue;
		if (ndim >= 3) {
			REC_SBT2SVT(subSBT, dim, ndim - 1, Rtype);
			continue;
		}
		/* 'subSBT' is a "leaf buffer". */
		lv = PROTECT(lb2lv(subSBT, Rtype));
		SET_VECTOR_ELT(SBT, i, lv);
		finalize_int_leaf_buffer(subSBT);
		UNPROTECT(1);
	}
	return;
}

/* Replace all "leaf buffers" with "leaf vectors" in Sparse Buffer Tree 'SBT'.
   'SBT' can NOT be R_NilValue! Must be called with 'ndim' >= 2. */
void _SBT2SVT(SEXP SBT, const int *dim, int ndim, SEXPTYPE Rtype)
{
	REC_SBT2SVT(SBT, dim, ndim, Rtype);
	return;
}

