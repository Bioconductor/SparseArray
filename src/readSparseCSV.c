/****************************************************************************
 *                     Workhorse behind readSparseCSV()                     *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "readSparseCSV.h"

#include "S4Vectors_interface.h"
#include "XVector_interface.h"

#include "leaf_vector_utils.h"
#include "ExtendableJaggedArray.h"

#include <R_ext/Connections.h>  /* for R_ReadConnection() */

#include <string.h>  /* for memcpy() */

#define	IOBUF_SIZE 8000002


/****************************************************************************
 * filexp_gets2(): A version of filexp_gets() that also works on connections
 * Copied from rtracklayer/src/readGFF.c
 */

Rconnection getConnection(int n);  /* not in <R_ext/Connections.h>, why? */

static char con_buf[250000];
static int con_buf_len, con_buf_offset;

static void init_con_buf()
{
	con_buf_len = con_buf_offset = 0;
	return;
}

static int filexp_gets2(SEXP filexp, char *buf, int buf_size, int *EOL_in_buf)
{
	Rconnection con;
	int buf_offset;
	char c;

	if (TYPEOF(filexp) == EXTPTRSXP)
		return filexp_gets(filexp, buf, buf_size, EOL_in_buf);
	buf_offset = *EOL_in_buf = 0;
	while (buf_offset < buf_size - 1) {
		if (con_buf_offset == con_buf_len) {
			con = getConnection(asInteger(filexp));
			con_buf_len = (int) R_ReadConnection(con,
					con_buf,
					sizeof(con_buf) / sizeof(char));
			if (con_buf_len == 0)
				break;
			con_buf_offset = 0;
		}
		c = con_buf[con_buf_offset++];
		buf[buf_offset++] = c;
		if (c == '\n') {
			*EOL_in_buf = 1;
			break;
		}
	}
	buf[buf_offset] = '\0';
	if (buf_offset == 0)
		return 0;
	if (con_buf_len == 0 || *EOL_in_buf)
		return 2;
	return 1;
}

/****************************************************************************
 * Using an environment as a growable list.
 */

static inline SEXP idx0_to_key(int idx0)
{
	char keybuf[11];

	snprintf(keybuf, sizeof(keybuf), "%010d", idx0);
	return mkChar(keybuf);  // unprotected!
}

static void set_env_elt(SEXP env, int i, SEXP val)
{
	SEXP key;

	key = PROTECT(idx0_to_key(i));
	defineVar(install(translateChar(key)), val, env);
	UNPROTECT(1);
	return;
}

static SEXP dump_env_as_list_or_R_NilValue(SEXP env, int ans_len)
{
	SEXP ans, key, ans_elt;
	int is_empty, i;

	ans = PROTECT(NEW_LIST(ans_len));
	is_empty = 1;
	for (i = 0; i < ans_len; i++) {
		key = PROTECT(idx0_to_key(i));
		ans_elt = findVar(install(translateChar(key)), env);
		UNPROTECT(1);
		if (ans_elt == R_UnboundValue)
			continue;
		SET_VECTOR_ELT(ans, i, ans_elt);
		is_empty = 0;
	}
	UNPROTECT(1);
	return is_empty ? R_NilValue : ans;
}


/****************************************************************************
 * Low-level helpers used by C_readSparseCSV_as_SVT_SparseMatrix()
 */

static char errmsg_buf[200];

static char get_sep_char(SEXP sep)
{
	SEXP sep0;

	if (!(IS_CHARACTER(sep) && LENGTH(sep) == 1))
		error("'sep' must be a single character");
	sep0 = STRING_ELT(sep, 0);
	if (sep0 == NA_STRING || length(sep0) != 1)
		error("'sep' must be a single character");
	return CHAR(sep0)[0];
}

static inline void IntAE_fast_append(IntAE *ae, int val)
{
	/* We don't use IntAE_get_nelt() for maximum speed. */
	if (ae->_nelt == ae->_buflength)
		IntAE_extend(ae, increase_buflength(ae->_buflength));
	ae->elts[ae->_nelt++] = val;
	return;
}

static inline void CharAEAE_fast_append(CharAEAE *aeae, CharAE *ae)
{
	/* We don't use CharAEAE_get_nelt() for maximum speed. */
	if (aeae->_nelt == aeae->_buflength)
		CharAEAE_extend(aeae, increase_buflength(aeae->_buflength));
	aeae->elts[aeae->_nelt++] = ae;
	return;
}

static void load_csv_rowname(const char *data, int data_len,
			     CharAEAE *csv_rownames_buf)
{
	CharAE *ae = new_CharAE(data_len);
	memcpy(ae->elts, data, data_len);
	/* We don't use CharAE_set_nelt() for maximum speed. */
	ae->_nelt = data_len;
	CharAEAE_fast_append(csv_rownames_buf, ae);
}


/****************************************************************************
 * C_readSparseCSV_as_SVT_SparseMatrix()
 */

/* 'offs_buf' and 'vals_buf' are **asumed** to have the same length and this
   length is **assumed** to be != 0. We don't check this! */
static SEXP make_leaf_vector_from_AEbufs(const IntAE *offs_buf,
					 const IntAE *vals_buf)
{
	SEXP ans_offs, ans_vals, ans;

	ans_offs = PROTECT(new_INTEGER_from_IntAE(offs_buf));
	ans_vals = PROTECT(new_INTEGER_from_IntAE(vals_buf));
	ans = _new_leaf_vector(ans_offs, ans_vals);  // unprotected!
	UNPROTECT(2);
	return ans;
}

static void store_AEbufs_in_env_as_leaf_vector(
		const IntAE *offs_buf, const IntAE *vals_buf,
		int idx0, SEXP env)
{
	int lv_len;
	SEXP lv;

	lv_len = IntAE_get_nelt(offs_buf);
	if (lv_len == 0)
		return;
	lv = PROTECT(make_leaf_vector_from_AEbufs(offs_buf, vals_buf));
	set_env_elt(env, idx0, lv);
	UNPROTECT(1);
	return;
}

static void load_csv_data_to_AEbufs(const char *data, int data_len,
		int off, IntAE *offs_buf, IntAE *vals_buf)
{
	int val;

	data_len = delete_trailing_LF_or_CRLF(data, data_len);
	if (data_len == 0 || (val = as_int(data, data_len)) == 0)
		return;
	IntAE_fast_append(offs_buf, off);
	IntAE_fast_append(vals_buf, val);
	return;
}

static void load_csv_data_to_offss_and_valss(const char *data, int data_len,
		int row_idx0, int col_idx0,
		ExtendableJaggedArray *offss, ExtendableJaggedArray *valss)
{
	int val;

	data_len = delete_trailing_LF_or_CRLF(data, data_len);
	if (data_len == 0 || (val = as_int(data, data_len)) == 0)
		return;
	_add_ExtendableJaggedArray_elt(offss, col_idx0, row_idx0);
	_add_ExtendableJaggedArray_elt(valss, col_idx0, val);
	return;
}

/* Used to load the sparse data when 'transpose' is TRUE. */
static void load_csv_row_to_AEbufs(const char *line, char sep,
		CharAEAE *csv_rownames_buf, IntAE *offs_buf, IntAE *vals_buf)
{
	int col_idx, i, data_len;
	const char *data;
	char c;

	IntAE_set_nelt(offs_buf, 0);
	IntAE_set_nelt(vals_buf, 0);
	col_idx = i = 0;
	data = line;
	data_len = 0;
	while ((c = line[i++])) {
		if (c != sep) {
			data_len++;
			continue;
		}
		if (col_idx == 0) {
			load_csv_rowname(data, data_len, csv_rownames_buf);
		} else {
			load_csv_data_to_AEbufs(data, data_len, col_idx - 1,
						offs_buf, vals_buf);
		}
		col_idx++;
		data = line + i;
		data_len = 0;
	}
	load_csv_data_to_AEbufs(data, data_len, col_idx - 1,
				offs_buf, vals_buf);
	return;
}

/* Used to load the sparse data when 'transpose' is FALSE. */
static void load_csv_row_to_offss_and_valss(const char *line, char sep,
		int row_idx0, CharAEAE *csv_rownames_buf,
		ExtendableJaggedArray *offss, ExtendableJaggedArray *valss)
{
	int col_idx, i, data_len;
	const char *data;
	char c;

	col_idx = i = 0;
	data = line;
	data_len = 0;
	while ((c = line[i++])) {
		if (c != sep) {
			data_len++;
			continue;
		}
		if (col_idx == 0) {
			load_csv_rowname(data, data_len, csv_rownames_buf);
		} else {
			load_csv_data_to_offss_and_valss(data, data_len,
						  row_idx0, col_idx - 1,
						  offss, valss);
		}
		col_idx++;
		data = line + i;
		data_len = 0;
	}
	load_csv_data_to_offss_and_valss(data, data_len,
				  row_idx0, col_idx - 1,
				  offss, valss);
	return;
}

static const char *read_sparse_csv(
		SEXP filexp, char sep, int transpose,
		CharAEAE *csv_rownames_buf,
		ExtendableJaggedArray *offss, ExtendableJaggedArray *valss,
		SEXP tmpenv)
{
	IntAE *offs_buf, *vals_buf;
	int row_idx0, lineno, ret_code, EOL_in_buf;
	char buf[IOBUF_SIZE];

	if (transpose) {
		offs_buf = new_IntAE(0, 0, 0);
		vals_buf = new_IntAE(0, 0, 0);
	}
	if (TYPEOF(filexp) != EXTPTRSXP)
		init_con_buf();
	row_idx0 = 0;
	for (lineno = 1;
	     (ret_code = filexp_gets2(filexp, buf, IOBUF_SIZE, &EOL_in_buf));
	     lineno += EOL_in_buf)
	{
		if (ret_code == -1) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "read error while reading characters "
				 "from line %d", lineno);
			return errmsg_buf;
		}
		if (!EOL_in_buf) {
			snprintf(errmsg_buf, sizeof(errmsg_buf),
				 "cannot read line %d, "
				 "line is too long", lineno);
			return errmsg_buf;
		}
		if (lineno == 1)
			continue;
		if (transpose) {
			/* Turn the CSV rows into leaf vectors as we go and
			   store them in 'tmpenv'. */
			load_csv_row_to_AEbufs(buf, sep,
					       csv_rownames_buf,
					       offs_buf, vals_buf);
			store_AEbufs_in_env_as_leaf_vector(
					       offs_buf, vals_buf,
					       row_idx0, tmpenv);
		} else {
			load_csv_row_to_offss_and_valss(buf, sep, row_idx0,
						csv_rownames_buf,
						offss, valss);
		}
		row_idx0++;
	}
	return NULL;
}

/* --- .Call ENTRY POINT ---
 * Args:
 *   filexp:    A "file external pointer" (see src/io_utils.c in the XVector
 *              package). TODO: Support connections (see src/readGFF.c in the
 *              rtracklayer package for how to do that).
 *   sep:       A single string (i.e. character vector of length 1) made of a
 *              single character.
 *   transpose: A single logical (TRUE or FALSE).
 *   csv_ncol:  Number of columns of data in the CSV file (1st column
 *              containing the rownames doesn't count).
 *   tmpenv:    An environment that will be used to grow the list
 *              of "leaf vectors" when 'transpose' is TRUE. Unused when
 *              'transpose' is FALSE.
 * Returns 'list(csv_rownames, SVT)'.
 */
SEXP C_readSparseCSV_as_SVT_SparseMatrix(SEXP filexp, SEXP sep,
					 SEXP transpose, SEXP csv_ncol,
					 SEXP tmpenv)
{
	int transpose0, nrow0, ncol0;
	CharAEAE *csv_rownames_buf;
	ExtendableJaggedArray offss, valss;
	const char *errmsg;
	SEXP ans, ans_elt;

	transpose0 = LOGICAL(transpose)[0];
	ncol0 = INTEGER(csv_ncol)[0];
	csv_rownames_buf = new_CharAEAE(0, 0);
	if (!transpose0) {
		offss = _new_ExtendableJaggedArray(ncol0);
		valss = _new_ExtendableJaggedArray(ncol0);
	}

	errmsg = read_sparse_csv(filexp, get_sep_char(sep), transpose0,
				 csv_rownames_buf, &offss, &valss,
				 tmpenv);
	if (errmsg != NULL)
		error("reading file: %s", errmsg);

	ans = PROTECT(NEW_LIST(2));

	ans_elt = PROTECT(new_CHARACTER_from_CharAEAE(csv_rownames_buf));
	SET_VECTOR_ELT(ans, 0, ans_elt);
	UNPROTECT(1);

	if (transpose0) {
		nrow0 = CharAEAE_get_nelt(csv_rownames_buf);
		ans_elt = dump_env_as_list_or_R_NilValue(tmpenv, nrow0);
	} else {
		ans_elt = _move_ExtendableJaggedArrays_to_SVT(&offss, &valss);
		_free_ExtendableJaggedArray(&offss);
		_free_ExtendableJaggedArray(&valss);
	}
	PROTECT(ans_elt);
	SET_VECTOR_ELT(ans, 1, ans_elt);
	UNPROTECT(1);

	UNPROTECT(1);
	return ans;
}

