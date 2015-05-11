/*
** $Id: t_error.c,v 1.1 2001/04/19 00:03:06 rosinski Exp $
*/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

/*
** t_error: error return routine to print a message and return a failure
** value.
**
** Input arguments:
**   fmt: format string
**   variable list of additional arguments for vfprintf
**
** Return value: -1 (failure)
*/

int t_error (const char *fmt, ...)
{
  va_list args;

#if ( ! defined DISABLE_TIMERS )
  va_start (args, fmt);

  if (fmt != NULL)
    (void) vfprintf (stderr, fmt, args);

  va_end (args);
#endif
  
  return (-1);
}
