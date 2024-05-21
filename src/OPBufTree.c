/****************************************************************************
 *                      OPBuf and OPBufTree structures                      *
 ****************************************************************************/
#include "OPBufTree.h"

#include <stdlib.h>  /* for malloc(), free(), calloc(), realloc() */
#include <limits.h>  /* for INT_MAX */
#include <errno.h>


static void alloc_error(int errnum)
{
	error("SparseArray internal error: %s", strerror(errnum));
}

static OPBuf *alloc_empty_OPBuf(void)
{
	OPBuf *opbuf = (OPBuf *) malloc(sizeof(OPBuf));
	if (opbuf == NULL)
		alloc_error(errno);
	opbuf->buflen = opbuf->nelt = 0;
	opbuf->idx0s = NULL;
	opbuf->Loffs = NULL;
	opbuf->xLoffs = NULL;
	return opbuf;
}

static void free_OPBuf(OPBuf *opbuf)
{
	if (opbuf->idx0s != NULL)
		free(opbuf->idx0s);
	if (opbuf->Loffs != NULL)
		free(opbuf->Loffs);
	if (opbuf->xLoffs != NULL)
		free(opbuf->xLoffs);
	free(opbuf);
}

static int increase_buflen(int buflen)
{
	if (buflen == INT_MAX)
		return -1;
	if (buflen == 0)
		return 1;
	if (buflen <= 2)
		return 4;
	if (buflen < 32768)			/* 2^15 */
		return buflen * 4;
	if (buflen < 33554432)			/* 2^25 */
		return buflen * 2;
	if (buflen < 268435456)			/* 2^28 */
		return buflen + 33554432;	/* 2^25 increments */
	if (buflen < INT_MAX - 536870912)	/* INT_MAX - 2^29 */
		return buflen + 268435456;	/* 2^28 increments */
	return INT_MAX;
}

static R_xlen_t *alloc_xLoffs_and_init_with_Loffs(int buflen,
						  int *Loffs, int nelt)
{
	R_xlen_t *xLoffs = (R_xlen_t *) malloc(sizeof(R_xlen_t) * buflen);
	if (xLoffs == NULL)
		alloc_error(errno);
	if (Loffs != NULL) {
		for (int k = 0; k < nelt; k++)
			xLoffs[k] = (R_xlen_t) Loffs[k];
		free(Loffs);
	}
	return xLoffs;
}

static int extend_OPBuf(OPBuf *opbuf, int extend_xLoffs)
{
	int new_buflen = increase_buflen(opbuf->buflen);
	if (new_buflen < 0)
		return MAX_OPBUF_LEN_REACHED;
	if (opbuf->buflen == 0) {
		opbuf->idx0s = (int *)
			malloc(sizeof(int) * new_buflen);
		if (opbuf->idx0s == NULL)
			alloc_error(errno);
		if (extend_xLoffs) {
			opbuf->xLoffs = (R_xlen_t *)
				malloc(sizeof(R_xlen_t) * new_buflen);
			if (opbuf->xLoffs == NULL)
				alloc_error(errno);
		} else {
			opbuf->Loffs = (int *)
				malloc(sizeof(int) * new_buflen);
			if (opbuf->Loffs == NULL)
				alloc_error(errno);
		}
	} else {
		int *new_idx0s = (int *)
			realloc(opbuf->idx0s, sizeof(int) * new_buflen);
		if (new_idx0s == NULL)
			alloc_error(errno);
		opbuf->idx0s = new_idx0s;
		if (opbuf->xLoffs != NULL) {
			R_xlen_t *new_xLoffs = (R_xlen_t *)
				realloc(opbuf->xLoffs,
					sizeof(R_xlen_t) * new_buflen);
			if (new_xLoffs == NULL)
				alloc_error(errno);
			opbuf->xLoffs = new_xLoffs;
		} else if (extend_xLoffs) {
			opbuf->xLoffs = alloc_xLoffs_and_init_with_Loffs(
						new_buflen,
						opbuf->Loffs, opbuf->nelt);
			opbuf->Loffs = NULL;
		} else {
			int *new_Loffs = (int *)
				realloc(opbuf->Loffs, sizeof(int) * new_buflen);
			if (new_Loffs == NULL)
				alloc_error(errno);
			opbuf->Loffs = new_Loffs;
		}
	}
	return opbuf->buflen = new_buflen;
}

/* This is the **unsafe** version of append_idx0xLoff_to_OPBuf() below.
   It assumes that 'opbuf' has not been switched from using 'opbuf->Loffs'
   to using 'opbuf->xLoffs' yet (i.e. that 'opbuf->xLoffs' is NULL). This
   is NOT checked! Slightly faster than append_idx0xLoff_to_OPBuf(). */
static int append_idx0Loff_to_OPBuf(OPBuf *opbuf, int idx0, int Loff)
{
	if (opbuf->nelt >= opbuf->buflen) {
		int ret = extend_OPBuf(opbuf, 0);
		if (ret < 0)
			return ret;
	}
	opbuf->idx0s[opbuf->nelt] = idx0;
	opbuf->Loffs[opbuf->nelt] = Loff;
	return ++(opbuf->nelt);
}

/* Safe version of append_idx0Loff_to_OPBuf(). Will switch 'opbuf' from
   using 'opbuf->Loffs' to using 'opbuf->xLoffs' if necessary. Note that
   the switch is irreversible. */
static int append_idx0xLoff_to_OPBuf(OPBuf *opbuf, int idx0, R_xlen_t Loff)
{
	if (opbuf->xLoffs == NULL && Loff <= INT_MAX)
		return append_idx0Loff_to_OPBuf(opbuf, idx0, (int) Loff);
	if (opbuf->nelt >= opbuf->buflen) {
		int ret = extend_OPBuf(opbuf, 1);
		if (ret < 0)
			return ret;
	} else if (opbuf->xLoffs == NULL) {
		opbuf->xLoffs = alloc_xLoffs_and_init_with_Loffs(
					opbuf->buflen,
					opbuf->Loffs, opbuf->nelt);
		opbuf->Loffs = NULL;
	}
	opbuf->idx0s[opbuf->nelt] = idx0;
	opbuf->xLoffs[opbuf->nelt] = Loff;
	return ++(opbuf->nelt);
}

static InnerNode *alloc_InnerNode(int n)
{
	InnerNode *inner_node = (InnerNode *) malloc(sizeof(InnerNode));
	if (inner_node == NULL)
		alloc_error(errno);
	inner_node->n = n;
	inner_node->children = (OPBufTree *) calloc(n, sizeof(OPBufTree));
	if (inner_node->children == NULL) {
		free(inner_node);
		alloc_error(errno);
	}
	return inner_node;
}

static void free_InnerNode(InnerNode *inner_node)
{
	for (int i = 0; i < inner_node->n; i++)
		_free_OPBufTree(inner_node->children + i);
	free(inner_node->children);
	free(inner_node);
	return;
}

void _alloc_OPBufTree_children(OPBufTree *opbuf_tree, int n)
{
	if (opbuf_tree->node_type != NULL_NODE)
		error("SparseArray internal error in "
		      "_alloc_OPBufTree_children():\n"
		      "    opbuf_tree->node_type != NULL_NODE");
	opbuf_tree->node.inner_node_p = alloc_InnerNode(n);
	opbuf_tree->node_type = INNER_NODE;  /* must be last */
	return;
}

static void alloc_OPBufTree_leaf(OPBufTree *opbuf_tree)
{
	if (opbuf_tree->node_type != NULL_NODE)
		error("SparseArray internal error in "
		      "alloc_OPBufTree_leaf():\n"
		      "    opbuf_tree->node_type != NULL_NODE");
	opbuf_tree->node.opbuf_p = alloc_empty_OPBuf();
	opbuf_tree->node_type = LEAF_NODE;  /* must be last */
	return;
}

/* 'opbuf_tree' must be a node of type NULL_NODE or LEAF_NODE.
   Based on **unsafe** append_idx0Loff_to_OPBuf().
   See append_idx0Loff_to_OPBuf() above for the details. */
int _append_idx0Loff_to_OPBufTree_leaf(OPBufTree *opbuf_tree,
		int idx0, int Loff)
{
	if (opbuf_tree->node_type == NULL_NODE)
		alloc_OPBufTree_leaf(opbuf_tree);
	OPBuf *opbuf = get_OPBufTree_leaf(opbuf_tree);
	return append_idx0Loff_to_OPBuf(opbuf, idx0, Loff);
}

/* 'opbuf_tree' must be a node of type NULL_NODE or LEAF_NODE. */
int _append_idx0xLoff_to_OPBufTree_leaf(OPBufTree *opbuf_tree,
		int idx0, R_xlen_t Loff)
{
	if (opbuf_tree->node_type == NULL_NODE)
		alloc_OPBufTree_leaf(opbuf_tree);
	OPBuf *opbuf = get_OPBufTree_leaf(opbuf_tree);
	return append_idx0xLoff_to_OPBuf(opbuf, idx0, Loff);
}

void _free_OPBufTree(OPBufTree *opbuf_tree)
{
	if (opbuf_tree->node_type == NULL_NODE)
		return;
	if (opbuf_tree->node_type == INNER_NODE) {
		free_InnerNode(opbuf_tree->node.inner_node_p);
	} else {
		free_OPBuf(opbuf_tree->node.opbuf_p);
	}
	*opbuf_tree = OPBufTree0;
	return;
}

static void print_OPBuf(OPBuf *opbuf, const char *margin)
{
	/* Print 'opbuf->idx0s'. */
	Rprintf("%sidx0s : ", margin);
	if (opbuf->idx0s == NULL) {
		Rprintf("NULL");
	} else {
		for (int k = 0; k < opbuf->nelt; k++)
			Rprintf("%4d", opbuf->idx0s[k]);
	}
	Rprintf("\n");

	/* Print 'opbuf->Loffs'. */
	Rprintf("%sLoffs : ", margin);
	if (opbuf->Loffs == NULL) {
		Rprintf("NULL");
	} else {
		for (int k = 0; k < opbuf->nelt; k++)
			Rprintf("%4d", opbuf->Loffs[k]);
	}
	Rprintf("\n");

	/* Print 'opbuf->xLoffs'. */
	Rprintf("%sxLoffs: ", margin);
	if (opbuf->xLoffs == NULL) {
		Rprintf("NULL");
	} else {
		for (int k = 0; k < opbuf->nelt; k++)
			Rprintf("%4lu", opbuf->xLoffs[k]);
	}
	Rprintf("\n");
	return;
}

void _print_OPBufTree(const OPBufTree *opbuf_tree, int depth)
{
	if (opbuf_tree->node_type == NULL_NODE) {
		Rprintf("NULL\n");
		return;
	}
	char format[14], margin[100];
	if (opbuf_tree->node_type == LEAF_NODE) {
		OPBuf *opbuf = get_OPBufTree_leaf(opbuf_tree);
		Rprintf("OPBuf (buflen=%d)\n", opbuf->buflen);
		snprintf(format, sizeof(format), "%%%ds", 2 * (depth + 1));
		snprintf(margin, sizeof(margin), format, "");
		print_OPBuf(opbuf, margin);
		return;
	}
	InnerNode *inner_node = opbuf_tree->node.inner_node_p;
	Rprintf("InnerNode\n");
	snprintf(format, sizeof(format), "%%%ds", 2 * depth);
	snprintf(margin, sizeof(margin), format, "");
	for (int i = 0; i < inner_node->n; i++) {
		Rprintf("%so child %d/%d: ", margin, i + 1, inner_node->n);
		_print_OPBufTree(inner_node->children + i, depth + 1);
	}
	return;
}


/****************************************************************************
 * Manipulation of the global OPBufTree
 */

static OPBufTree global_opbuf_tree = OPBufTree0;

OPBufTree *_get_global_opbuf_tree(void)
{
	return &global_opbuf_tree;
}

/* --- .Call ENTRY POINT --- */
SEXP C_free_global_OPBufTree(void)
{
	_free_OPBufTree(&global_opbuf_tree);
	return R_NilValue;
}

