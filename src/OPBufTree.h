#ifndef _OPBUF_TREE_H_
#define _OPBUF_TREE_H_

#include <Rdefines.h>


#define	MAX_OPBUF_LEN_REACHED -1

/* Buffer of Offset Pairs. Each (idx0,Loff) pair is made of:
   - idx0: an offset along the first dimension of the array to subset or
           subassign;
   - Loff: an offset along the L-index used for the subsetting or
           subassignment. */
typedef struct opbuf_t {
	int buflen;
	int *idx0s;        /* Array of offsets < dim(<array>)[0]  */
	int *Loffs;        /* Array of offsets < length(<L-index>).
			      Note that we could use unsigned int instead of
			      int to delay the switch to 'xLoffs' as much
			      as possible. */
	R_xlen_t *xLoffs;  /* Use instead of 'Loffs' to store "long offsets" */
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

static inline int get_OPBufTree_nchildren(const OPBufTree *opbuf_tree)
{
	if (opbuf_tree->node_type != INNER_NODE)
		error("SparseArray internal error in "
		      "get_OPBufTree_nchildren():\n"
		      "    opbuf_tree->node_type != INNER_NODE");
	return opbuf_tree->node.inner_node_p->n;
}

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

void _alloc_OPBufTree_children(
	OPBufTree *opbuf_tree,
	int n
);

int _append_idx0Loff_to_host_node(
	OPBufTree *host_node,
	int idx0,
	int Loff
);

int _append_idx0xLoff_to_host_node(
	OPBufTree *host_node,
	int idx0,
	R_xlen_t Loff
);

void _free_OPBufTree(OPBufTree *opbuf_tree);

void _print_OPBufTree(
	const OPBufTree *opbuf_tree,
	int depth
);

OPBufTree *_get_global_opbuf_tree(void);

SEXP C_free_global_OPBufTree(void);

#endif  /* _OPBUF_TREE_H_ */

