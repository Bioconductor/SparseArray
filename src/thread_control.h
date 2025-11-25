#ifndef _THREAD_CONTROL_H_
#define _THREAD_CONTROL_H_

#include <Rdefines.h>

SEXP C_get_num_procs(void);

SEXP C_get_max_threads(void);

SEXP C_set_max_threads(SEXP nthread);

#endif  /* _THREAD_CONTROL_H_ */

