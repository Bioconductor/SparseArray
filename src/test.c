/****************************************************************************
 ****************************************************************************/
#include "test.h"

#define	IOBUF_SIZE 8000002


static char readcon_buf[250000];
static int con_buf_len, con_buf_offset;

static void init_con_buf()
{
	con_buf_len = con_buf_offset = 0;
	return;
}

static FILE *open_log(const char *logfile, const char *mode)
{
	char filepath[256];
	int ret;
	FILE *f;

	ret = snprintf(filepath, sizeof(filepath),
		       "/home/hpages/github/Bioconductor/%s", logfile);
	if (ret >= sizeof(filepath))
		error("output of snprintf() got truncated");
	f = fopen(filepath, mode);
	if (f == NULL)
		error("failed to open %s", filepath);
	return f;
}

static void start_log(const char *logfile)
{
	FILE *f = open_log(logfile, "w");
	Rprintf("--- start %s ---\n", logfile);
	fprintf(f, "--- start %s ---\n", logfile);
	fclose(f);
	return;
}

static void write_log(const char *logfile, const char *text)
{
	FILE *f = open_log(logfile, "a");
	Rprintf("%s: %s", logfile, text);
	fprintf(f, "%s", text);
	fclose(f);
	return;
}

/* Work-in-progress: filexp_gets3() tries to reproduce the mysterious "segfault
   from C stack overflow" error that we get when readSparseCSV() is called
   in the context of creating the vignette with 'R CMD build', but with
   C code that is as simplified as possible. */
static int filexp_gets3(char *buf, int buf_size, int *EOL_in_buf)
{
	int buf_offset;
	char c;

	buf_offset = *EOL_in_buf = 0;
	start_log("debug-filexp_gets3.log");
	while (buf_offset < buf_size - 1) {
		Rprintf("buf_offset = %d -- buf_size = %d\n",
			buf_offset, buf_size);
		Rprintf("con_buf_offset = %d -- con_buf_len = %d\n",
			con_buf_offset, con_buf_len);
		//write_log("debug-filexp_gets3.log", "ok1\n");
		if (con_buf_offset == con_buf_len) {
			//write_log("debug-filexp_gets3.log", "ok2\n");
			con_buf_len = snprintf(readcon_buf, sizeof(readcon_buf),
					       "%s\n", ",A,B,C,D,E,F");
			Rprintf("con_buf_len = %d\n", con_buf_len);
			if (con_buf_len == 0)
				break;
			con_buf_offset = 0;
		}
		c = readcon_buf[con_buf_offset++];
		Rprintf("%d\n", (int) c);
		buf[buf_offset++] = c;
		if (c == '\n') {
			*EOL_in_buf = 1;
			break;
		}
	}
	write_log("debug-filexp_gets3.log", "DONE\n");
	buf[buf_offset] = '\0';
	if (buf_offset == 0)
		return 0;
	if (con_buf_len == 0 || *EOL_in_buf)
		return 2;
	return 1;
}

/* --- .Call ENTRY POINT --- */
SEXP C_test()
{

	int EOL_in_buf;
	char buf[IOBUF_SIZE];

	start_log("debug-Call.log");
	write_log("debug-Call.log", "ok0\n");
	init_con_buf();
	int ret_code = filexp_gets3(buf, IOBUF_SIZE, &EOL_in_buf);
	write_log("debug-Call.log", "DONE\n");
	return NEW_LIST(2);
}

