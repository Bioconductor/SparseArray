#ifndef _OPBUF_TREE_H_
#define _OPBUF_TREE_H_

#include <Rdefines.h>


/* Buffer of Offset Pairs. Each pair is made of a "long offset" (R_xlen_t)
   and a "short offset" (int). */
typedef struct opbuf_t {
	int buflen;
	R_xlen_t *loffs;  /* array of long offsets */
	int *soffs;       /* array of short offsets */
	int nelt;
} OPBuf;

enum node_type { NULL_NODE, INNER_NODE, LEAF_NODE };

/* Tree of OPBuf structures (as leaves). */
typedef struct opbuf_tree_t {
	enum node_type node_type;
	union {
		struct inner_node_t {
			int n;
			struct opbuf_tree_t *children;  /* array */
		} *inner_node_p;  /* inner node: array of OPBufTree's */
		OPBuf *opbuf_p;   /* leaf node: pointer to an OPBuf */
	} node;
} OPBufTree;

typedef struct inner_node_t InnerNode;

static const OPBufTree OPBufTree0 = { NULL_NODE, {NULL}};

static inline OPBufTree *get_OPBufTree_child(const OPBufTree *opbuf_tree, int i)
{
	if (opbuf_tree->node_type != INNER_NODE)
		error("SparseArray internal error in "
		      "get_OPBufTree_child():\n"
		      "    opbuf_tree->node_type != INNER_NODE");
	return opbuf_tree->node.inner_node_p->children + i;
}

static inline OPBuf *get_OPBufTree_leaf(const OPBufTree *opbuf_tree)
{
	if (opbuf_tree->node_type != LEAF_NODE)
		error("SparseArray internal error in "
		      "get_OPBufTree_leaf():\n"
		      "    opbuf_tree->node_type != LEAF_NODE");
	return opbuf_tree->node.opbuf_p;
}

int _append_to_OPBuf(
	OPBuf *opbuf,
	R_xlen_t loff,
	int soff
);

void _alloc_OPBufTree_children(
	OPBufTree *opbuf_tree,
	int n
);

void _alloc_OPBufTree_leaf(OPBufTree *opbuf_tree);

void _free_OPBufTree(OPBufTree *opbuf_tree);

void _print_OPBufTree(
	const OPBufTree *opbuf_tree,
	int depth
);

OPBufTree *_get_global_opbuf_tree(void);

SEXP C_free_global_OPBufTree(void);

#endif  /* _OPBUF_TREE_H_ */

