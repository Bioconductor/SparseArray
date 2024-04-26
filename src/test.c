/****************************************************************************
 ****************************************************************************/
#include "test.h"

/*
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
*/

#define	IOBUF_SIZE 8000002

static const char *string1 = "ABCDEF";

static void write_to_buf(char *buf, int buf_size)
{
	int offset;
	char c;

	offset = 0;
	while ((c = string1[offset])) {
		Rprintf("%c\n", c);
		buf[offset] = c;
		offset++;
	}
	return;
}

/* --- .Call ENTRY POINT ---
   Reproduces the mysterious "segfault from C stack overflow" error that
   we use to get when readSparseCSV() was called in the context of creating
   the vignette with 'R CMD build', but with C code that is much simpler.
   We observe this segfault if 'buf' is declared as:

       char buf[IOBUF_SIZE];

   but not if it's declared as:

       static char buf[IOBUF_SIZE];

   So some strange memory corruption seems to happen if 'buf' is not declared
   as static. Note that we can also observe this memory corruption problem
   when running the following code interactively:

       library(SparseArray)
       ## More that one call to SparseArray:::test() might be needed!
       SparseArray:::test()
       SparseArray:::test()
       for (i in 1:10) SparseArray:::test()

   In this case, the error message is slightly different: "C stack usage
   8033324 is too close to the limit".

   readSparseCSV() was fixed in 0.99.4 by declaring 'buf' as static in
   read_sparse_csv() (see src/readSparseCSV.c).
*/
SEXP C_test(void)
{
	/* Not using the 'static' keyword will produce the memory corruption
	   described above. */
	char buf[IOBUF_SIZE];
	//static char buf[IOBUF_SIZE];

	//start_log("debug-Call.log");
	write_to_buf(buf, IOBUF_SIZE);
	//write_log("debug-Call.log", "DONE\n");
	return R_NilValue;
}

