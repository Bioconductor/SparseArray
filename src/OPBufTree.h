#ifndef _OPBUF_TREE_H_
#define _OPBUF_TREE_H_

#include <Rdefines.h>


/****************************************************************************
 * OPBuf and OPBufTree structures and API
 *
 * The OPBuf/OPBufTree structs are used by the subsetting and subassignment
 * code in SparseArray_subsetting.c and SparseArray_subassignment.c.
 */

#define	MAX_OPBUF_LEN_REACHED -1

/* Buffer of Offset Pairs. Each (idx0,Loff) pair is made of:
   - idx0: an offset along the first dimension of the array to subset or
           subassign;
   - Loff: an offset along the L-index used for the subsetting or
           subassignment.
   Note that 'buflen' and 'nelt' are both of type int below. This means
   that the max length of an OPBuf is INT_MAX. */
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


/****************************************************************************
 * extract_long_idx0() and family
 *
 * Additional helper functions also used by the subsetting and subassignment
 * code in SparseArray_subsetting.c and SparseArray_subassignment.c.
 */

/* In addition to 0, extract_long_idx0() and extract_idx0() can return
   one of the five negative values defined below. Note that these are
   typically considered errors but not always. For example, some functions
   in SparseArray_subsetting.c treat SUBSCRIPT_ELT_IS_NA as an error but
   others don't.
   Be aware that MAX_OPBUF_LEN_REACHED is set to -1 (see above in this file)
   so do **not** use that value. */
#define BAD_SUBSCRIPT_TYPE             -2
#define SUBSCRIPT_IS_TOO_LONG          -3  /* returned by extract_idx0() only */
#define SUBSCRIPT_ELT_IS_LESS_THAN_ONE -4
#define SUBSCRIPT_ELT_IS_BEYOND_MAX    -5
#define SUBSCRIPT_ELT_IS_NA            -6

/* 'subscript' must be a numeric vector, possibly a long one. It is
   expected to contain 1-based indices that are >= 1 and <= 'max'.
   Returns 0 or one of the negative values defined above. Will set 'idx0'
   only if returning 0. In other words, caller **must** ignore 'idx0' if
   a non-zero value is returned. */
static inline int extract_long_idx0(SEXP subscript, R_xlen_t i,
				    R_xlen_t max, R_xlen_t *idx0)
{
	if (IS_INTEGER(subscript)) {
		int idx = INTEGER(subscript)[i];
		if (idx == NA_INTEGER)
			return SUBSCRIPT_ELT_IS_NA;
		idx--;  /* from 1-based to 0-based */
		if (idx < 0)
			return SUBSCRIPT_ELT_IS_LESS_THAN_ONE;
		if ((R_xlen_t) idx >= max)  /* >= yes! */
			return SUBSCRIPT_ELT_IS_BEYOND_MAX;
		*idx0 = (R_xlen_t) idx;
		return 0;
	}
	if (IS_NUMERIC(subscript)) {
		double idx = REAL(subscript)[i];
		/* ISNAN(): True for *both* NA and NaN. See <R_ext/Arith.h> */
		if (ISNAN(idx))
			return SUBSCRIPT_ELT_IS_NA;
		idx -= 1.0;  /* from 1-based to 0-based */
		if (idx < 0.0)
			return SUBSCRIPT_ELT_IS_LESS_THAN_ONE;
		if (idx >= (double) max)  /* >= yes! */
			return SUBSCRIPT_ELT_IS_BEYOND_MAX;
		*idx0 = (R_xlen_t) idx;
		return 0;
	}
	return BAD_SUBSCRIPT_TYPE;
}

/* Like extract_long_idx0(), except that 'subscript' cannot be a long vector
   nor can it contain values > INT_MAX ('max' must be supplied as an 'int'). */
static inline int extract_idx0(SEXP subscript, int i, int max, int *idx0)
{
	if (XLENGTH(subscript) > (R_xlen_t) INT_MAX)
		return SUBSCRIPT_IS_TOO_LONG;
	R_xlen_t lidx0;
	int ret = extract_long_idx0(subscript, (R_xlen_t) i,
				    (R_xlen_t) max, &lidx0);
	if (ret < 0)
		return ret;
	*idx0 = (int) lidx0;
	return 0;
}

void _bad_Lindex_error(int ret_code);
void _bad_Mindex_error(int ret_code);
void _bad_Nindex_error(int ret_code, int along1);

#endif  /* _OPBUF_TREE_H_ */

