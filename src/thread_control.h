#ifndef _THREAD_CONTROL_H_
#define _THREAD_CONTROL_H_

#include <Rdefines.h>

int _get_max_threads(void);

SEXP C_get_num_procs(void);

SEXP C_get_max_threads(void);

SEXP C_set_max_threads(SEXP nthread);

SEXP C_get_initial_device(void);

SEXP C_pause_resource(SEXP hard_pause, SEXP device_num);

#endif  /* _THREAD_CONTROL_H_ */

