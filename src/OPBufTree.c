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

OPBuf *new_empty_OPBuf(void)
{
	OPBuf *opbuf = (OPBuf *) malloc(sizeof(OPBuf));
	if (opbuf == NULL)
		alloc_error(errno);
	opbuf->buflen = opbuf->nelt = 0;
	opbuf->loffs = NULL;
	opbuf->soffs = NULL;
	return opbuf;
}

static void free_OPBuf(OPBuf *opbuf)
{
	if (opbuf->loffs != NULL)
		free(opbuf->loffs);
	if (opbuf->soffs != NULL)
		free(opbuf->soffs);
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

static int extend_OPBuf(OPBuf *opbuf)
{
	int new_buflen = increase_buflen(opbuf->buflen);
	if (new_buflen < 0)
		return MAX_OPBUF_LEN_REACHED;
	if (opbuf->buflen == 0) {
		opbuf->loffs = (R_xlen_t *)
			malloc(sizeof(R_xlen_t) * new_buflen);
		if (opbuf->loffs == NULL)
			alloc_error(errno);
		opbuf->soffs = (int *)
			malloc(sizeof(int) * new_buflen);
		if (opbuf->soffs == NULL)
			alloc_error(errno);
	} else {
		R_xlen_t *new_loffs = (R_xlen_t *)
			realloc(opbuf->loffs, sizeof(R_xlen_t) * new_buflen);
		if (new_loffs == NULL)
			alloc_error(errno);
		opbuf->loffs = new_loffs;
		int *new_soffs = (int *)
			realloc(opbuf->soffs, sizeof(int) * new_buflen);
		if (new_soffs == NULL)
			alloc_error(errno);
		opbuf->soffs = new_soffs;
	}
	return opbuf->buflen = new_buflen;
}

int _append_to_OPBuf(OPBuf *opbuf, R_xlen_t loff, int soff)
{
	if (opbuf->nelt >= opbuf->buflen) {
		int ret = extend_OPBuf(opbuf);
		if (ret < 0)
			return ret;
	}
	opbuf->loffs[opbuf->nelt] = loff;
	opbuf->soffs[opbuf->nelt] = soff;
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
	InnerNode *inner_node = alloc_InnerNode(n);
	opbuf_tree->node.inner_node_p = inner_node;
	opbuf_tree->node_type = INNER_NODE;
	return;
}

void _alloc_OPBufTree_leaf(OPBufTree *opbuf_tree)
{
	if (opbuf_tree->node_type != NULL_NODE)
		error("SparseArray internal error in "
		      "_alloc_OPBufTree_leaf():\n"
		      "    opbuf_tree->node_type != NULL_NODE");
	OPBuf *opbuf = new_empty_OPBuf();
	opbuf_tree->node.opbuf_p = opbuf;
	opbuf_tree->node_type = LEAF_NODE;
	return;
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
	Rprintf("%sloffs: ", margin);
	for (int k = 0; k < opbuf->nelt; k++)
		Rprintf("%4lu", opbuf->loffs[k]);
	Rprintf("\n");
	Rprintf("%ssoffs: ", margin);
	for (int k = 0; k < opbuf->nelt; k++)
		Rprintf("%4d", opbuf->soffs[k]);
	Rprintf("\n");
	return;
}

void _print_OPBufTree(const OPBufTree *opbuf_tree, int depth)
{
	if (opbuf_tree->node_type == NULL_NODE) {
		Rprintf("NULL\n");
		return;
	}
	char format[10], margin[100];
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
	for (int i = 0; i < inner_node->n; i++) {
		snprintf(format, sizeof(format), "%%%ds", 2 * depth);
		snprintf(margin, sizeof(margin), format, "");
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

